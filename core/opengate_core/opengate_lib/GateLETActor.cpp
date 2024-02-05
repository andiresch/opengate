/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   ------------------------------------ -------------- */

#include "GateLETActor.h"
#include "G4Navigator.hh"
#include "G4RandomTools.hh"
#include "G4RunManager.hh"
#include "GateHelpers.h"
#include "GateHelpersDict.h"
#include "GateHelpersImage.h"

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

#include "G4LinInterpolation.hh"

// Mutex that will be used by thread to write in the edep/dose image
G4Mutex SetLETPixelMutex = G4MUTEX_INITIALIZER;

GateLETActor::GateLETActor(py::dict &user_info) : GateVActor(user_info, true) {
  // Create the image pointer
  // (the size and allocation will be performed on the py side)
  cpp_numerator_image = ImageType::New();
  cpp_denominator_image = ImageType::New();
  // Action for this actor: during stepping
  fActions.insert("SteppingAction");
  fActions.insert("BeginOfRunAction");
  fActions.insert("EndSimulationAction");
  // Option: compute uncertainty
  fmMKM = DictGetBool(user_info, "enable_rbe");
  fdoseAverage = DictGetBool(user_info, "dose_average");
  ftrackAverage = DictGetBool(user_info, "track_average");
  fLETtoOtherMaterial = DictGetBool(user_info, "let_to_other_material");
  fotherMaterial = DictGetStr(user_info, "other_material");
  // fQAverage = DictGetBool(user_info, "qAverage");
  fInitialTranslation = DictGetG4ThreeVector(user_info, "translation");
  // Hit type (random, pre, post etc)
  fHitType = DictGetStr(user_info, "hit_type");
  // Option: RBE model type (mkm, etc)
  if (fmMKM) {
    // create lookuptable
    // energies = new G4DataVector;
    table = new std::vector<G4DataVector *>;
    CreateLookupTable(user_info);
    fRBEmodel = DictGetStr(user_info, "rbe_model");
    if (fRBEmodel == "mkm") {
      fAlpha0 = DictGetDouble(user_info, "alpha_0");
      fBeta = DictGetDouble(user_info, "beta");
    }
  }
}

void GateLETActor::ActorInitialize() {}

void GateLETActor::BeginOfRunAction(const G4Run *) {
  // Important ! The volume may have moved, so we re-attach each run
  AttachImageToVolume<ImageType>(cpp_numerator_image, fPhysicalVolumeName,
                                 fInitialTranslation);
  AttachImageToVolume<ImageType>(cpp_denominator_image, fPhysicalVolumeName,
                                 fInitialTranslation);
  // compute volume of a dose voxel
  auto sp = cpp_numerator_image->GetSpacing();
  fVoxelVolume = sp[0] * sp[1] * sp[2];
  static G4Material *water =
      G4NistManager::Instance()->FindOrBuildMaterial(fotherMaterial);
}

void GateLETActor::SteppingAction(G4Step *step) {
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

  // get pixel index
  ImageType::IndexType index;
  bool isInside =
      cpp_numerator_image->TransformPhysicalPointToIndex(point, index);

  // set value
  if (isInside) {
    // With mutex (thread)
    G4AutoLock mutex(&SetLETPixelMutex);

    // get edep in MeV (take weight into account)
    auto w = step->GetTrack()->GetWeight();
    auto edep = step->GetTotalEnergyDeposit() / CLHEP::MeV * w;
    double dedx_cut = DBL_MAX;
    // dedx
    auto *current_material = step->GetPreStepPoint()->GetMaterial();
    auto density = current_material->GetDensity() / CLHEP::g * CLHEP::cm3;
    // double dedx_currstep = 0., dedx_water = 0.;
    // double density_water = 1.0;
    //  other material
    const G4ParticleDefinition *p = step->GetTrack()->GetParticleDefinition();

    auto energy1 = step->GetPreStepPoint()->GetKineticEnergy() / CLHEP::MeV;
    auto energy2 = step->GetPostStepPoint()->GetKineticEnergy() / CLHEP::MeV;
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
    auto dedx_currstep =
        l.ComputeElectronicDEDX(energy, p, current_material, dedx_cut) /
        CLHEP::MeV * CLHEP::mm;

    auto steplength = step->GetStepLength() / CLHEP::mm;
    double scor_val_num = 0.;
    double scor_val_den = 0.;

    if (fLETtoOtherMaterial) {
      auto density_water = water->GetDensity() / CLHEP::g * CLHEP::cm3;
      auto dedx_water = l.ComputeElectronicDEDX(energy, p, water, dedx_cut) /
                        CLHEP::MeV * CLHEP::mm;
      auto SPR_otherMaterial = dedx_water / dedx_currstep;
      edep *= SPR_otherMaterial;
      dedx_currstep *= SPR_otherMaterial;
    }

    if (fdoseAverage) {
      scor_val_num = edep * dedx_currstep / CLHEP::MeV / CLHEP::MeV * CLHEP::mm;
      scor_val_den = edep / CLHEP::MeV;
    } else if (ftrackAverage) {
      scor_val_num = steplength * dedx_currstep * w / CLHEP::MeV;
      scor_val_den = steplength * w / CLHEP::mm;
    } else if (fmMKM) {
      const G4ParticleDefinition *p = step->GetTrack()->GetParticleDefinition();
      if (p == G4Gamma::Gamma())
        p = G4Electron::Electron();
      /*auto dedx_currstep =
          emcalc->ComputeElectronicDEDX(energy, p, current_material, dedx_cut) /
          CLHEP::MeV * CLHEP::mm;*/
      auto charge = int(p->GetAtomicNumber());
      auto mass = p->GetAtomicMass();
      auto table_value = GetValue(charge, energy); // energy has unit?
      //       auto alpha_currstep = fAlpha0 + fBeta * table_value;
      auto alpha_currstep = table_value;
      /*
      ========== GATE 9.3
          if(mDoseEnergyByZ[i][k]>energy){
                                            efficiency=mDoseEfficiencyByZ[i][k-1]+(mDoseEfficiencyByZ[i][k]-mDoseEfficiencyByZ[i][k-1])/(mDoseEnergyByZ[i][k]-mDoseEnergyByZ[i][k-1])*(energy-mDoseEnergyByZ[i][k-1]);
                                            k=mDoseEnergyByZ[i].size();
                                          }
                                          else
      if(k==mDoseEnergyByZ[i].size()-1){ GateMessage("Actor", 0, "WARNING
      particle energy larger than energies available in the file:
      "<<mDoseEfficiencyFileByZ.at(i)<<" Efficiency = 1 instead"<<Gateendl);
                                          }
                                          ======= GATE 9.3

                                  */
      //            std::cout<< "energy:" << energy << ", mass: " << mass <<
      //            std::endl; std::cout << "Charge: " << charge << ",
      //            energy/mass: " << energy/mass << std::endl; std::cout
      //            <<"z*_1D: " << table_value <<
      //            ", alpha_step: " << alpha_currstep<< std::endl;

      // auto steplength = step->GetStepLength() / CLHEP::mm;
      //       double scor_val_num = 0.;
      //       double scor_val_den = 0.;
      if (alpha_currstep > 0) {
        scor_val_num = edep * alpha_currstep;
        scor_val_den = edep;
      } else {
        scor_val_num = 0.;
        scor_val_den = 0.;
      }
    }
    ImageAddValue<ImageType>(cpp_numerator_image, index, scor_val_num);
    ImageAddValue<ImageType>(cpp_denominator_image, index, scor_val_den);
    //}

  } // else : outside the image
}

// ========================== copied from RBE Actor =====================

void GateLETActor::CreateLookupTable(py::dict &user_info) {
  // get lookup table
  std::vector<std::vector<double>> lookupTab =
      DictGetVecofVecDouble(user_info, "lookup_table");
  // energies = VectorToG4DataVector(lookupTab[0]);

  for (int i = 0; i < lookupTab.size(); i++) {
    table->push_back(VectorToG4DataVector(lookupTab[i]));
  }
}

double GateLETActor::GetValue(int Z, float energy) {
  // std::cout << "GetValue: Z: " << Z << ", energy[MeV/u]: " << energy <<
  // std::endl;
  // initalize value
  G4double y = 0;
  // get table values for the given Z
  //   if (Z > 6 || Z < 1 ){
  // 	  return 0;}
  //   G4DataVector *data = (*table)[Z - 1];
  G4DataVector *Z_vec = new G4DataVector();
  Z_vec->insertAt(0, Z);
  int bin_table = -1;
  G4DataVector *energies;
  G4DataVector *data;

  //   std::cout<<"Get data for Z : "<<Z <<std::endl;
  for (int i = 0; i < table->size(); i += 3) {
    //   std::cout<<"*(*table)["<< i << "][0]: "<<(*(*table)[i])[0] <<std::endl;
    if (*(*table)[i] == *Z_vec) {
      bin_table = i;
      energies = (*table)[i + 1];
      data = (*table)[i + 2];
      //       std::cout<<"Z_vec: "<<Z_vec <<std::endl;
      //       std::cout<<"energies: "<< energies<<std::endl;
      //       std::cout<<"data: "<< data <<std::endl;
    }
  }
  if (bin_table == -1) {
    return 0;
  }
  // find the index of the lower bound energy to the given energy
  size_t bin = FindLowerBound(energy, energies);
  // std::cout << "interpolation bin: " << bin << std::endl;
  G4LinInterpolation linearAlgo;
  // get table value for the given energy
  y = linearAlgo.Calculate(energy, bin, *energies, *data);
  // std::cout<<"interpolation output:" << y << std::endl;

  return y;
}

size_t GateLETActor::FindLowerBound(G4double x, G4DataVector *values) const {
  size_t lowerBound = 0;
  size_t upperBound(values->size() - 1);
  if (x < (*values)[0]) {
    return 0;
  }
  if (x > (*values).back()) {
    return values->size() - 1;
  }
  while (lowerBound <= upperBound) {
    size_t midBin((lowerBound + upperBound) / 2);
    // std::cout<<"upper: "<<upperBound<<" lower: "<<lowerBound<<std::endl;
    // std::cout<<(*values)[midBin]<<std::endl;
    if (x < (*values)[midBin])
      upperBound = midBin - 1;
    else
      lowerBound = midBin + 1;
  }
  return upperBound;
}

void GateLETActor::EndSimulationAction() {}
