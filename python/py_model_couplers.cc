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
/* -------------------------------------------------------------------------- */
#include <aka_error.hh>
#include <cohesive_contact_solvercallback.hh>
#include <coupler_solid_cohesive_contact.hh>
#include <coupler_solid_contact.hh>
#include <element_synchronizer.hh>
#include <node_synchronizer.hh>
#include <non_linear_solver.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/attr.h>
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace {
template <class CouplerSolidContact_>
auto register_coupler_solid_contact(py::module & mod, const std::string & name)
    -> py::class_<CouplerSolidContact_, Model> {
  return py::class_<CouplerSolidContact_, Model>(mod, name.c_str(),
                                                 py::multiple_inheritance())
      .def(py::init<Mesh &, UInt, const ID &, std::shared_ptr<DOFManager>>(),
           py::arg("mesh"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("id") = "coupler_solid_contact",
           py::arg("dof_manager") = nullptr, py::keep_alive<1, 2>())
      .def("applyBC",
           [](CouplerSolidContact_ & self,
              BC::Dirichlet::DirichletFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("applyBC",
           [](CouplerSolidContact_ & self, BC::Neumann::NeumannFunctor & func,
              const std::string & element_group) {
             self.applyBC(func, element_group);
           })
      .def("getContactMechanicsModel",
           &CouplerSolidContact_::getContactMechanicsModel,
           py::return_value_policy::reference)
      .def("getDisplacement", &CouplerSolidContact_::getDisplacement,
           py::return_value_policy::reference)
      .def("getVelocity", &CouplerSolidContact_::getVelocity,
           py::return_value_policy::reference)
      .def("getAcceleration", &CouplerSolidContact_::getAcceleration,
           py::return_value_policy::reference)
      .def("getMass", &CouplerSolidContact_::getMass,
           py::return_value_policy::reference)
      .def("getContactForce", &CouplerSolidContact_::getContactForce,
           py::return_value_policy::reference)
      .def("getContactDetector", &CouplerSolidContact_::getContactDetector,
           py::return_value_policy::reference)
      .def("setTimeStep", &CouplerSolidContact_::setTimeStep,
           py::arg("time_step"), py::arg("solver_id") = "")
      .def("getStableTimeStep", &CouplerSolidContact_::getStableTimeStep);
}
} // namespace

/* -------------------------------------------------------------------------- */
void register_model_couplers(py::module & mod) {
  register_coupler_solid_contact<CouplerSolidContact>(mod,
                                                      "CouplerSolidContact")
      .def(
          "getSolidMechanicsModel",
          [](CouplerSolidContact & self) -> decltype(auto) {
            return self.getSolidMechanicsModel();
          },
          py::return_value_policy::reference)
      .def(
          "initFull",
          [](CouplerSolidContact & self,
             const AnalysisMethod & analysis_method) {
            self.initFull(_analysis_method = analysis_method);
          },
          py::arg("_analysis_method") = _explicit_lumped_mass);

  register_coupler_solid_contact<CouplerSolidCohesiveContact>(
      mod, "CouplerSolidCohesiveContact")
      .def(
          "initFull",
          [](CouplerSolidCohesiveContact & self,
             const AnalysisMethod & analysis_method, bool is_extrinsic) {
            self.initFull(_analysis_method = analysis_method,
                          _is_extrinsic = is_extrinsic);
          },
          py::arg("_analysis_method") = _explicit_lumped_mass,
          py::arg("_is_extrinsic") = false)
      .def("checkCohesiveStress",
           [](CouplerSolidCohesiveContact & self) {
             return self.checkCohesiveStress();
           })
      .def(
          "getSolidMechanicsModelCohesive",
          [](CouplerSolidCohesiveContact & self) -> decltype(auto) {
            return self.getSolidMechanicsModelCohesive();
          },
          py::return_value_policy::reference);
}

} // namespace akantu
