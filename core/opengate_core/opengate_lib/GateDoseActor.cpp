/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   ------------------------------------ -------------- */

#include "GateDoseActor.h"
#include "G4Navigator.hh"
#include "G4RandomTools.hh"
#include "G4RunManager.hh"

#include "G4Deuteron.hh"
#include "G4Electron.hh"
#include "G4EmCalculator.hh"
#include "G4Gamma.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"

#include "GateHelpers.h"
#include "GateHelpersDict.h"
#include "GateHelpersImage.h"

#include <iostream>

// Mutex that will be used by thread to write in the edep/dose image
G4Mutex SetPixelMutex = G4MUTEX_INITIALIZER;

GateDoseActor::GateDoseActor(py::dict &user_info)
    : GateVActor(user_info, true) {
  // Create the image pointer
  // (the size and allocation will be performed on the py side)
  cpp_edep_image = ImageType::New();
  // Action for this actor: during stepping
  fActions.insert("SteppingAction");
  fActions.insert("BeginOfRunAction");
  fActions.insert("EndSimulationAction");
  // Option: compute uncertainty
  fUncertaintyFlag = DictGetBool(user_info, "uncertainty");
  // Option: compute dose in Gray
  fGrayFlag = DictGetBool(user_info, "gray");
  // Option: compute dose in Gray
  fDoseToWaterFlag = DictGetBool(user_info, "dose_to_water");
  // translation
  fInitialTranslation = DictGetG4ThreeVector(user_info, "translation");
  // Hit type (random, pre, post etc)
  fHitType = DictGetStr(user_info, "hit_type");
}

void GateDoseActor::ActorInitialize() {
  if (fUncertaintyFlag) {
    cpp_square_image = ImageType::New();
    cpp_temp_image = ImageType::New();
    cpp_last_id_image = ImageType::New();
  }
  if (fDoseToWaterFlag) {
    // emcalc = new G4EmCalculator; // removed now
    fGrayFlag = true;
  }
  if (fGrayFlag) {
    cpp_dose_image = ImageType::New();
  }
}

void GateDoseActor::BeginOfRunAction(const G4Run *) {
  // Important ! The volume may have moved, so we re-attach each run
  AttachImageToVolume<ImageType>(cpp_edep_image, fPhysicalVolumeName,
                                 fInitialTranslation);
  // compute volume of a dose voxel
  auto sp = cpp_edep_image->GetSpacing();
  fVoxelVolume = sp[0] * sp[1] * sp[2];
}

void GateDoseActor::SteppingAction(G4Step *step) {
  auto preGlobal = step->GetPreStepPoint()->GetPosition();
  auto postGlobal = step->GetPostStepPoint()->GetPosition();
  auto touchable = step->GetPreStepPoint()->GetTouchable();

  // FIXME If the volume has multiple copy, touchable->GetCopyNumber(0) ?

  // consider random position between pre and post
  auto position = postGlobal;
  if (fHitType == "pre") {
    position = preGlobal;
  }
  if (fHitType == "random") {
    auto x = G4UniformRand();
    auto direction = postGlobal - preGlobal;
    position = preGlobal + x * direction;
  }
  if (fHitType == "middle") {
    auto direction = postGlobal - preGlobal;
    position = preGlobal + 0.5 * direction;
  }
  auto localPosition =
      touchable->GetHistory()->GetTransform(0).TransformPoint(position);

  // convert G4ThreeVector to itk PointType
  ImageType::PointType point;
  point[0] = localPosition[0];
  point[1] = localPosition[1];
  point[2] = localPosition[2];

  // get edep in MeV (take weight into account)
  auto w = step->GetTrack()->GetWeight();
  auto edep = step->GetTotalEnergyDeposit() / CLHEP::MeV * w;

  // get pixel index
  ImageType::IndexType index;
  bool isInside = cpp_edep_image->TransformPhysicalPointToIndex(point, index);

  // set value
  if (isInside) {
    // With mutex (thread)
    G4AutoLock mutex(&SetPixelMutex);

    // If uncertainty: consider edep per event
    if (fUncertaintyFlag) {
      auto event_id =
          G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
      auto previous_id = cpp_last_id_image->GetPixel(index);
      cpp_last_id_image->SetPixel(index, event_id);
      if (event_id == previous_id) {
        // Same event : continue temporary edep
        ImageAddValue<ImageType>(cpp_temp_image, index, edep);
      } else {
        // Different event : update previous and start new event
        auto e = cpp_temp_image->GetPixel(index);
        ImageAddValue<ImageType>(cpp_edep_image, index, e);
        ImageAddValue<ImageType>(cpp_square_image, index, e * e);
        // new temp value
        cpp_temp_image->SetPixel(index, edep);
      }
    } else {
      ImageAddValue<ImageType>(cpp_edep_image, index, edep);
    }

    // Compute the dose in Gray ?
    if (fGrayFlag) {
      auto *material_currstep = step->GetPreStepPoint()->GetMaterial();
      auto density_currstep = material_currstep->GetDensity();
      auto dose = edep / density_currstep / fVoxelVolume / CLHEP::gray;

      if (fDoseToWaterFlag) {
        double dedx_cut = DBL_MAX;
        // dedx
        double dedx_currstep = 0., dedx_water = 0.;
        double density_water = 1.0;
        // other material
        const G4ParticleDefinition *p =
            step->GetTrack()->GetParticleDefinition();
        static G4Material *water =
            G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER");
        auto energy1 = step->GetPreStepPoint()->GetKineticEnergy();
        auto energy2 = step->GetPostStepPoint()->GetKineticEnergy();
        auto energy = (energy1 + energy2) / 2;
        // Accounting for particles with dedx=0; i.e. gamma and neutrons
        // For gamma we consider the dedx of electrons instead - testing
        // with 1.3 MeV photon beam or 150 MeV protons or 1500 MeV carbon ion
        // beam showed that the error induced is 0 		when comparing
        // dose and dosetowater in the material G4_WATER For neutrons the dose
        // is neglected - testing with 1.3 MeV photon beam or 150 MeV protons or
        // 1500 MeV carbon ion beam showed that the error induced is < 0.01%
        //		when comparing dose and dosetowater in the material
        // G4_WATER (we are systematically missing a little bit of dose of
        // course with this solution)

        if (p == G4Gamma::Gamma())
          p = G4Electron::Electron();
        auto &l = fThreadLocalData.Get().emcalc;

        dedx_currstep =
            l.ComputeTotalDEDX(energy, p, material_currstep, dedx_cut);
        dedx_water = l.ComputeTotalDEDX(energy, p, water, dedx_cut);
        // dedx_currstep =
        //     emcalc->ComputeTotalDEDX(energy, p, material_currstep, dedx_cut);
        // dedx_water = emcalc->ComputeTotalDEDX(energy, p, water, dedx_cut);
        density_water = water->GetDensity();
        double spr = dedx_currstep / dedx_water;
        double mspr =
            (density_currstep / density_water) * (dedx_water / dedx_currstep);
        // std::cout <<"density_currstep: " << density_currstep
        // *(CLHEP::g/CLHEP::cm3)<< spr<< std::endl;
        /*
        std::cout <<"water name: " <<water->GetName() << std::endl;
        std::cout <<"water getDensity: " <<water->GetDensity() << std::endl;
        std::cout <<"mat name: " << material_currstep->GetName() << std::endl;
        std::cout <<"mat temperature: " << material_currstep->GetTemperature()
        << std::endl; std::cout <<"mat MassMolecu: " <<
        material_currstep->GetMassOfMolecule() << std::endl; std::cout <<"mat
        getdensity: " << material_currstep->GetDensity() << std::endl; std::cout
        <<"density_currstep: " << density_currstep << std::endl; std::cout
        <<"SPR: " << spr<< std::endl; std::cout <<"mSPR: " << 1/mspr<<
        std::endl<< std::endl;
        */

        // In current implementation, dose deposited directly by neutrons is
        // neglected - the below lines prevent "inf or NaN"
        if (dedx_currstep == 0 || dedx_water == 0) {
          dose = 0.;
        } else {
          // std::cout << "Overwrite dose: "<< std::endl;
          dose *=
              (density_currstep / density_water) * (dedx_water / dedx_currstep);
        }
      }

      // std::cout << "Add to image: "<< std::endl;
      ImageAddValue<ImageType>(cpp_dose_image, index, dose);
    }

  } // else : outside the image
}

void GateDoseActor::EndSimulationAction() {}
