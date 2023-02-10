/* --------------------------------------------------
   Copyright (C): OpenGATE Collaboration
   This software is distributed under the terms
   of the GNU Lesser General  Public Licence (LGPL)
   See LICENSE.md for further details
   -------------------------------------------------- */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include "GateDosimetryActor.h"

void init_GateDosimetryActor(py::module &m) {
  py::class_<GateDosimetryActor,
             std::unique_ptr<GateDosimetryActor, py::nodelete>, GateVActor>(
      m, "GateDosimetryActor")
      .def(py::init<py::dict &>())
      .def_readwrite("cpp_edep_image", &GateDosimetryActor::cpp_edep_image)
      .def_readwrite("cpp_square_image", &GateDosimetryActor::cpp_square_image)
      .def_readwrite("cpp_temp_image", &GateDosimetryActor::cpp_temp_image)
      .def_readwrite("cpp_dose_image", &GateDosimetryActor::cpp_dose_image)
      .def_readwrite("cpp_last_id_image",
                     &GateDosimetryActor::cpp_last_id_image)
      .def_readwrite("fPhysicalVolumeName",
                     &GateDosimetryActor::fPhysicalVolumeName);
}
