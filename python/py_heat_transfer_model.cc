/**
 * Copyright (©) 2019-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * This file is part of Akantu
 *
 * Akantu is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 */

/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
#include "py_constitutive_law.hh"
#include "py_constitutive_laws_handler.hh"
/* -------------------------------------------------------------------------- */
#include <diffusion_law.hh>
#include <diffusion_model.hh>
#include <heat_diffusion.hh>
#include <heat_transfer_model.hh>
#include <model_options.hh>
#include <non_linear_solver.hh>
/* -------------------------------------------------------------------------- */
// #include <pybind11/operators.h>
#include <pybind11/pybind11.h>
// #include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace {
template <typename _DiffusionLaw> class PyDiffusionLaw : public _DiffusionLaw {
public:
  /* Inherit the constructors */
  using _DiffusionLaw::_DiffusionLaw;

  void computeDiffusivityGradUOnQuadPoints(ElementType type,
                                           GhostType ghost_type) override {
    PYBIND11_OVERRIDE(void, _DiffusionLaw, computeDiffusivityGradUOnQuadPoints,
                      type, ghost_type);
  }

  void computeDiffusivityOnQuadPoints(ElementType type,
                                      GhostType ghost_type) override {
    PYBIND11_OVERRIDE(void, _DiffusionLaw, computeDiffusivityOnQuadPoints, type,
                      ghost_type);
  }

  Real getStableTimeStep(Real element_size) override {
    PYBIND11_OVERRIDE(Real, _DiffusionLaw, getStableTimeStep, element_size);
  }
};

/* ------------------------------------------------------------------------ */
template <typename _DiffusionLaw>
void register_diffusion_law_classes(py::module & mod, const ID & name) {
  py::class_<_DiffusionLaw, DiffusionLaw, PyDiffusionLaw<_DiffusionLaw>>(
      mod, name.c_str(), py::multiple_inheritance())
      .def(py::init<DiffusionModel &, const ID &>());
}
} // namespace
/* -------------------------------------------------------------------------- */
void register_heat_transfer_model(py::module & mod) {
  register_constitutive_law<DiffusionModel>(mod);

  /* ------------------------------------------------------------------------ */
  py::class_<DiffusionLaw, ConstitutiveLaw<DiffusionModel>>(
      mod, "DiffusionLaw", py::multiple_inheritance());

  register_diffusion_law_classes<HeatDiffusion<1>>(mod, "HeatDiffusion1D");
  register_diffusion_law_classes<HeatDiffusion<2>>(mod, "HeatDiffusion2D");
  register_diffusion_law_classes<HeatDiffusion<3>>(mod, "HeatDiffusion3D");

  /* ------------------------------------------------------------------------ */
  register_constitutive_laws_handler<DiffusionLaw, Model>(mod);

  /* ------------------------------------------------------------------------ */
  py::class_<DiffusionModel, ConstitutiveLawsHandler<DiffusionLaw, Model>>(
      mod, "DiffusionModel", py::multiple_inheritance())
      .def(
          "initFull",
          [](DiffusionModel & self, const AnalysisMethod & _analysis_method) {
            self.initFull(HeatTransferModelOptions(_analysis_method));
          },
          py::arg("_analysis_method"))
      .def("setTimeStep", &DiffusionModel::setTimeStep, py::arg("time_step"),
           py::arg("solver_id") = "")
      .def("getStableTimeStep", &DiffusionModel::getStableTimeStep)
      .def("getBlockedDOFs", &DiffusionModel::getBlockedDOFs);

  /* ------------------------------------------------------------------------ */
  py::class_<HeatTransferModel, DiffusionModel>(mod, "HeatTransferModel",
                                                py::multiple_inheritance())
      .def(py::init<Mesh &, Int, const ID &>(), py::arg("mesh"),
           py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "diffusion_model", py::keep_alive<1, 2>())
      .def(
          "getTemperature",
          [](HeatTransferModel & self) -> decltype(auto) {
            return self.getTemperature();
          },
          py::return_value_policy::reference);
}

} // namespace akantu
