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
/* -------------------------------------------------------------------------- */
#include <constitutive_law.hh>
#include <material_selector.hh>
#include <solid_mechanics_model.hh>
#if defined(AKANTU_COHESIVE_ELEMENT)
#include <material_cohesive_linear.hh>
#include <material_cohesive_linear_friction.hh>
#include <solid_mechanics_model_cohesive.hh>
#endif
#include <material_elastic.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace {
template <typename _Material> class PyMaterial : public _Material {
public:
  /* Inherit the constructors */
  using _Material::_Material;

  void initMaterial() override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(void, _Material, initMaterial, );
  };
  void computeStress(ElementType el_type,
                     GhostType ghost_type = _not_ghost) override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(void, _Material, computeStress, el_type, ghost_type);
  }
  void computeTangentModuli(ElementType el_type, Array<Real> & tangent_matrix,
                            GhostType ghost_type = _not_ghost) override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(void, _Material, computeTangentModuli, el_type,
                      tangent_matrix, ghost_type);
  }

  void computePotentialEnergy(ElementType el_type) override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(void, _Material, computePotentialEnergy, el_type);
  }

  [[nodiscard]] auto getPushWaveSpeed(const Element & element) const
      -> Real override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(Real, _Material, getPushWaveSpeed, element);
  }

  [[nodiscard]] auto getShearWaveSpeed(const Element & element) const
      -> Real override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(Real, _Material, getShearWaveSpeed, element);
  }

  [[nodiscard]] auto hasStiffnessMatrixChanged() const -> bool override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(bool, _Material, hasStiffnessMatrixChanged, );
  }
};

/* ------------------------------------------------------------------------ */
template <typename _Material>
void register_material_classes(py::module & mod, const std::string & name) {
  py::class_<_Material, Material, PyMaterial<_Material>>(
      mod, name.c_str(), py::multiple_inheritance())
      .def(py::init<SolidMechanicsModel &, const ID &>());
}

#if defined(AKANTU_COHESIVE_ELEMENT)
// trampoline for the cohesive materials
template <typename _Material,
          std::enable_if_t<std::is_base_of_v<MaterialCohesive, _Material>> * =
              nullptr>
class PyMaterialCohesive : public PyMaterial<_Material> {
public:
  using Parent = PyMaterial<_Material>;
  /* Inherit the constructors */
  using Parent::Parent;

  void checkInsertion(bool check_only) override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(void, _Material, checkInsertion, check_only);
  }

  void computeTraction(ElementType el_type,
                       GhostType ghost_type = _not_ghost) override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(void, _Material, computeTraction, el_type, ghost_type);
  }

  void computeTangentTraction(ElementType el_type, Array<Real> & tangent_matrix,
                              GhostType ghost_type = _not_ghost) override {
    // NOLINTNEXTLINE
    PYBIND11_OVERRIDE(void, _Material, computeTangentTraction, el_type,
                      tangent_matrix, ghost_type);
  }

  void computeStress(ElementType /*el_type*/,
                     GhostType /*ghost_type*/ = _not_ghost) final {}

  void computeTangentModuli(ElementType /*el_type*/,
                            Array<Real> & /*tangent_matrix*/,
                            GhostType /*ghost_type*/ = _not_ghost) final {}
};

template <typename _Material>
void register_material_cohesive_classes(py::module & mod,
                                        const std::string & name) {
  py::class_<_Material, MaterialCohesive, PyMaterialCohesive<_Material>>(
      mod, name.c_str(), py::multiple_inheritance())
      .def(py::init<SolidMechanicsModelCohesive &, const ID &>());
}
#endif
} // namespace

/* -------------------------------------------------------------------------- */
void register_material(py::module & mod) {
  register_constitutive_law<SolidMechanicsModel>(mod);

  py::class_<Material, PyMaterial<Material>,
             ConstitutiveLaw<SolidMechanicsModel>>(mod, "Material",
                                                   py::multiple_inheritance())
      .def(py::init<SolidMechanicsModel &, const ID &>())
      .def(
          "getGradU",
          [](Material & self, ElementType el_type,
             GhostType ghost_type = _not_ghost) -> decltype(auto) {
            return self.getGradU(el_type, ghost_type);
          },
          py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "getStress",
          [](Material & self, ElementType el_type,
             GhostType ghost_type = _not_ghost) -> decltype(auto) {
            return self.getStress(el_type, ghost_type);
          },
          py::arg("el_type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "getPotentialEnergy",
          [](Material & self, ElementType el_type) -> decltype(auto) {
            return self.getPotentialEnergy(el_type);
          },
          py::arg("el_type"), py::return_value_policy::reference)
      .def(
          "getPotentialEnergy",
          [](Material & self, Element element) -> Real {
            return self.getPotentialEnergy(element);
          },
          py::arg("element"))
      .def("getPotentialEnergy",
           [](Material & self) -> Real { return self.getPotentialEnergy(); })
      .def("initMaterial", &Material::initMaterial)
      .def("getModel",
           [](Material & self) -> SolidMechanicsModel & {
             return self.getModel();
           })
      .def("getPushWaveSpeed", &Material::getPushWaveSpeed)
      .def("getShearWaveSpeed", &Material::getShearWaveSpeed)
      .def("hasMatrixChanged", &Material::hasMatrixChanged)
      .def("hasStiffnessMatrixChanged", &Material::hasStiffnessMatrixChanged);

  register_material_classes<MaterialElastic<1>>(mod, "MaterialLinearElastic1D");
  register_material_classes<MaterialElastic<2>>(mod, "MaterialLinearElastic2D");
  register_material_classes<MaterialElastic<3>>(mod, "MaterialLinearElastic3D");

#if defined(AKANTU_COHESIVE_ELEMENT)
  /* ------------------------------------------------------------------------ */
  py::class_<MaterialCohesive, Material, PyMaterialCohesive<MaterialCohesive>>(
      mod, "MaterialCohesive", py::multiple_inheritance())
      .def(py::init<SolidMechanicsModelCohesive &, const ID &>())
      .def(
          "getFacetFilter",
          [](MaterialCohesive & self) -> decltype(auto) {
            return self.getFacetFilter();
          },
          py::return_value_policy::reference)
      .def(
          "getFacetFilter",
          [](MaterialCohesive & self, ElementType type,
             GhostType ghost_type) -> decltype(auto) {
            return self.getFacetFilter(type, ghost_type);
          },
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "getOpening",
          [](MaterialCohesive & self, ElementType type, GhostType ghost_type)
              -> decltype(auto) { return self.getOpening(type, ghost_type); },
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "getTraction",
          [](MaterialCohesive & self, ElementType type, GhostType ghost_type)
              -> decltype(auto) { return self.getTraction(type, ghost_type); },
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference);

  register_material_cohesive_classes<MaterialCohesiveLinear<2>>(
      mod, "MaterialCohesiveLinear2D");
  register_material_cohesive_classes<MaterialCohesiveLinear<3>>(
      mod, "MaterialCohesiveLinear3D");
  register_material_cohesive_classes<MaterialCohesiveLinearFriction<2>>(
      mod, "MaterialCohesiveLinearFriction2D");
  register_material_cohesive_classes<MaterialCohesiveLinearFriction<3>>(
      mod, "MaterialCohesiveLinearFriction3D");
#endif

  py::class_<MaterialFactory>(mod, "MaterialFactory")
      .def_static(
          "getInstance",
          []() -> MaterialFactory & { return Material::getFactory(); },
          py::return_value_policy::reference)
      .def("registerAllocator",
           [](MaterialFactory & self, const std::string id, py::function func) {
             // Keep the Python function alive by storing it
             auto stored_func = py::reinterpret_borrow<py::function>(func);
             stored_func.inc_ref();

             self.registerAllocator(
                 id,
                 [stored_func](Int dim, const ID & option,
                               SolidMechanicsModel & model,
                               const ID & id) -> std::unique_ptr<Material> {
                   py::object obj = stored_func(dim, option, model, id);
                   auto & ptr = py::cast<Material &>(obj);

                   obj.release();
                   return std::unique_ptr<Material>(&ptr);
                 });
           })
      .def("getPossibleAllocators", &MaterialFactory::getPossibleAllocators);
}

} // namespace akantu
