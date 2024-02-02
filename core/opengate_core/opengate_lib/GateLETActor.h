/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#ifndef GateLETActor_h
#define GateLETActor_h

#include "G4Cache.hh"
#include "G4EmCalculator.hh"
#include "G4NistManager.hh"
#include "G4VPrimitiveScorer.hh"
#include "GateVActor.h"
#include "itkImage.h"
#include <pybind11/stl.h>

#include "G4DataVector.hh"

namespace py = pybind11;

class GateLETActor : public GateVActor {

public:
  // Constructor
  GateLETActor(py::dict &user_info);

  virtual void ActorInitialize();

  // Main function called every step in attached volume
  virtual void SteppingAction(G4Step *);

  // Called every time a Run starts (all threads)
  virtual void BeginOfRunAction(const G4Run *run);

  virtual void EndSimulationAction();

  // Image type is 3D float by default
  // TODO double precision required
  typedef itk::Image<double, 3> ImageType;

  // The image is accessible on py side (shared by all threads)
  ImageType::Pointer cpp_numerator_image;
  ImageType::Pointer cpp_denominator_image;

  // Option: indicate if we must compute dose in Gray also
  std::string fPhysicalVolumeName;
  bool fdoseAverage;
  bool ftrackAverage;
  bool fLETtoOtherMaterial;
  std::string fotherMaterial;
  // RBE specific
  std::string fRBEmodel;

  double fAlpha0;
  double fBeta;
  bool fmMKM;

private:
  // RBE specific
  // G4DataVector *energies;
  std::vector<G4DataVector *> *table;

  void CreateLookupTable(py::dict &user_info);
  double GetValue(int Z, float energy);
  size_t FindLowerBound(G4double x, G4DataVector *values) const;

  double fVoxelVolume;

  G4ThreeVector fInitialTranslation;
  std::string fHitType;

  G4Material *water;

  struct threadLocalT {
    G4EmCalculator emcalc;
  };
  G4Cache<threadLocalT> fThreadLocalData;
};

#endif // GateLETActor_h
