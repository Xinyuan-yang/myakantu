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
#include "aka_config.hh"
/* -------------------------------------------------------------------------- */
#include "py_aka_array.hh"
/* -------------------------------------------------------------------------- */
#include <element_type_map.hh>
#include <fe_engine.hh>
#include <mesh.hh>
#include <mesh_accessor.hh>
#include <mesh_utils.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

namespace {
/* ------------------------------------------------------------------------ */
template <typename T>
auto register_element_type_map_array(py::module & mod,
                                     const std::string & name) {
  return py::class_<ElementTypeMapArray<T>,
                    std::shared_ptr<ElementTypeMapArray<T>>>(
             mod, ("ElementTypeMapArray" + name).c_str())
      .def(py::init<const ID &, const ID &>(),
           py::arg("id") = "by_element_type_array",
           py::arg("parent_id") = "no_parent")
      .def(
          "__call__",
          [](ElementTypeMapArray<T> & self, ElementType type,
             GhostType ghost_type) -> decltype(auto) {
            return self(type, ghost_type);
          },
          py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference, py::keep_alive<0, 1>())
      .def(
          "elementTypes",
          [](ElementTypeMapArray<T> & self, Int _dim, GhostType _ghost_type,
             ElementKind _kind) -> std::vector<ElementType> {
            auto types = self.elementTypes(_dim, _ghost_type, _kind);
            std::vector<ElementType> _types;
            for (auto && t : types) {
              _types.push_back(t);
            }
            return _types;
          },
          py::arg("dim") = _all_dimensions, py::arg("ghost_type") = _not_ghost,
          py::arg("kind") = _ek_regular)
      .def(
          "initialize",
          [](ElementTypeMapArray<T> & self, const Mesh & mesh,
             GhostType ghost_type, Int nb_component, Int spatial_dimension,
             ElementKind element_kind, bool with_nb_element,
             bool with_nb_nodes_per_element, T default_value,
             bool do_not_default) {
            self.initialize(
                mesh, _ghost_type = ghost_type, _nb_component = nb_component,
                _spatial_dimension =
                    (spatial_dimension == -2 ? mesh.getSpatialDimension()
                                             : spatial_dimension),
                _element_kind = element_kind,
                _with_nb_element = with_nb_element,
                _with_nb_nodes_per_element = with_nb_nodes_per_element,
                _default_value = default_value,
                _do_not_default = do_not_default);
          },
          py::arg("mesh"), py::kw_only(), py::arg("ghost_type") = _casper,
          py::arg("nb_component") = 1, py::arg("spatial_dimension") = -2,
          py::arg("element_kind") = _ek_not_defined,
          py::arg("with_nb_element") = false,
          py::arg("with_nb_nodes_per_element") = false,
          py::arg("default_value") = T{}, py::arg("do_not_default") = false)
      .def(
          "initialize",
          [](ElementTypeMapArray<T> & self, const FEEngine & fe_engine,
             GhostType ghost_type, Int nb_component, Int spatial_dimension,
             ElementKind element_kind, bool with_nb_element,
             bool with_nb_nodes_per_element, T default_value,
             bool do_not_default) {
            self.initialize(fe_engine, _ghost_type = ghost_type,
                            _nb_component = nb_component,
                            _spatial_dimension =
                                (spatial_dimension == -2
                                     ? fe_engine.getMesh().getSpatialDimension()
                                     : spatial_dimension),
                            _element_kind = element_kind,
                            _with_nb_element = with_nb_element,
                            _with_nb_nodes_per_element =
                                with_nb_nodes_per_element,
                            _default_value = default_value,
                            _do_not_default = do_not_default);
          },
          py::arg("mesh"), py::kw_only(), py::arg("ghost_type") = _casper,
          py::arg("nb_component") = 1, py::arg("spatial_dimension") = -2,
          py::arg("element_kind") = _ek_not_defined,
          py::arg("with_nb_element") = false,
          py::arg("with_nb_nodes_per_element") = false,
          py::arg("default_value") = T{}, py::arg("do_not_default") = false)
      .def(
          "initialize",
          [](ElementTypeMapArray<T> & self, const FEEngine & fe_engine,
             GhostType ghost_type, Int nb_component, Int spatial_dimension,
             ElementKind element_kind, bool with_nb_element,
             bool with_nb_nodes_per_element, T default_value,
             bool do_not_default) {
            self.initialize(fe_engine, _ghost_type = ghost_type,
                            _nb_component = nb_component,
                            _spatial_dimension =
                                (spatial_dimension == -2
                                     ? fe_engine.getMesh().getSpatialDimension()
                                     : spatial_dimension),
                            _element_kind = element_kind,
                            _with_nb_element = with_nb_element,
                            _with_nb_nodes_per_element =
                                with_nb_nodes_per_element,
                            _default_value = default_value,
                            _do_not_default = do_not_default);
          },
          py::arg("fe_engine"), py::arg("ghost_type") = _casper,
          py::arg("nb_component") = 1, py::arg("spatial_dimension") = -2,
          py::arg("element_kind") = _ek_not_defined,
          py::arg("with_nb_element") = false,
          py::arg("with_nb_nodes_per_element") = false,
          py::arg("default_value") = T{}, py::arg("do_not_default") = false)
      .def("getID", &ElementTypeMapArray<T>::getID);
}
} // namespace
/* -------------------------------------------------------------------------- */
void register_mesh(py::module & mod) {
  py::class_<Mesh::PeriodicSlaves>(mod, "PeriodicSlaves")
      .def(
          "__iter__",
          [](Mesh::PeriodicSlaves & _this) {
            return py::make_iterator(_this.begin(), _this.end());
          },
          py::keep_alive<0, 1>());

  py::class_<MeshData, std::shared_ptr<MeshData>>(mod, "MeshData")
      .def(
          "getElementalDataUInt",
          [](MeshData & _this, const ID & name) -> decltype(auto) {
            return _this.getElementalData<UInt>(name);
          },
          py::return_value_policy::reference)
      .def(
          "getElementalDataReal",
          [](MeshData & _this, const ID & name) -> decltype(auto) {
            return _this.getElementalData<Real>(name);
          },
          py::return_value_policy::reference);

  py::class_<Mesh, GroupManager, Dumpable, MeshData, std::shared_ptr<Mesh>>(
      mod, "Mesh", py::multiple_inheritance())
      .def(py::init<Int, const ID &>(), py::arg("spatial_dimension"),
           py::arg("id") = "mesh")
      .def("read", &Mesh::read, py::arg("filename"),
           py::arg("mesh_io_type") = _miot_auto, "read the mesh from a file")
      .def(
          "getNodes",
          [](Mesh & self) -> decltype(auto) { return self.getNodes(); },
          py::return_value_policy::reference)
      .def("getNbNodes", &Mesh::getNbNodes)
      .def(
          "getConnectivity",
          [](Mesh & self, ElementType type) -> decltype(auto) {
            return self.getConnectivity(type);
          },
          py::return_value_policy::reference)
      .def(
          "getConnectivities",
          [](Mesh & self) -> decltype(auto) {
            return self.getConnectivities();
          },
          py::return_value_policy::reference)
      .def(
          "addConnectivityType",
          [](Mesh & self, ElementType type, GhostType ghost_type) -> void {
            self.addConnectivityType(type, ghost_type);
          },
          py::arg("type"), py::arg("ghost_type") = _not_ghost)
      .def("distribute", [](Mesh & self) { self.distribute(); })
      .def("isDistributed",
           [](const Mesh & self) { return self.isDistributed(); })
      .def("fillNodesToElements", &Mesh::fillNodesToElements,
           py::arg("dimension") = _all_dimensions)
      .def("getAssociatedElements",
           [](Mesh & self, const Idx & node, py::list list) {
             Array<Element> elements;
             self.getAssociatedElements(node, elements);
             for (auto && element : elements) {
               list.append(element);
             }
           })
      .def("makePeriodic",
           [](Mesh & self, const SpatialDirection & direction) {
             self.makePeriodic(direction);
           })
      .def(
          "getNbElement",
          [](Mesh & self, const Int spatial_dimension, GhostType ghost_type,
             ElementKind kind) {
            return self.getNbElement(spatial_dimension, ghost_type, kind);
          },
          py::kw_only(), py::arg("spatial_dimension") = _all_dimensions,
          py::arg("ghost_type") = _not_ghost, py::arg("kind") = _ek_not_defined)
      .def(
          "getNbElement",
          [](Mesh & self, ElementType type, GhostType ghost_type) {
            return self.getNbElement(type, ghost_type);
          },
          py::kw_only(), py::arg("type"), py::arg("ghost_type") = _not_ghost)
      .def_static(
          "getSpatialDimension",
          [](ElementType & type) { return Mesh::getSpatialDimension(type); })
      .def(
          "getDataReal",
          [](Mesh & _this, const ID & name, ElementType type,
             GhostType ghost_type) -> decltype(auto) {
            return _this.getData<Real>(name, type, ghost_type);
          },
          py::arg("name"), py::arg("type"), py::arg("ghost_type") = _not_ghost,
          py::return_value_policy::reference)
      .def(
          "hasDataReal",
          [](Mesh & _this, const ID & name, ElementType type,
             GhostType ghost_type) -> bool {
            return _this.hasData<Real>(name, type, ghost_type);
          },
          py::arg("name"), py::arg("type"), py::arg("ghost_type") = _not_ghost)
      .def("isPeriodic", [](const Mesh & _this) { return _this.isPeriodic(); })
      .def("getPeriodicMaster", &Mesh::getPeriodicMaster)
      .def("getPeriodicSlaves", &Mesh::getPeriodicSlaves)
      .def("isPeriodicSlave", &Mesh::isPeriodicSlave)
      .def("isPeriodicMaster", &Mesh::isPeriodicMaster)
      .def(
          "getMeshFacets",
          [](const Mesh & self) -> const Mesh & {
            return self.getMeshFacets();
          },
          py::return_value_policy::reference)
      .def("getUpperBounds", &Mesh::getUpperBounds)
      .def("getLowerBounds", &Mesh::getLowerBounds)
      .def("initMeshFacets", &Mesh::initMeshFacets,
           py::arg("id") = "mesh_facets", py::return_value_policy::reference);

  /* ------------------------------------------------------------------------ */
  py::class_<MeshUtils>(mod, "MeshUtils")
      .def_static("buildFacets", &MeshUtils::buildFacets);

  py::class_<MeshAccessor>(mod, "MeshAccessor")
      .def(py::init<Mesh &>(), py::arg("mesh"))
      .def(
          "resizeConnectivity",
          [](MeshAccessor & self, Int new_size, ElementType type, GhostType gt)
              -> void { self.resizeConnectivity(new_size, type, gt); },
          py::arg("new_size"), py::arg("type"),
          py::arg("ghost_type") = _not_ghost)
      .def(
          "resizeNodes",
          [](MeshAccessor & self, Int new_size) -> void {
            self.resizeNodes(new_size);
          },
          py::arg("new_size"))
      .def("makeReady", &MeshAccessor::makeReady);

  register_element_type_map_array<Real>(mod, "Real");
  register_element_type_map_array<UInt>(mod, "UInt");
  register_element_type_map_array<Int>(mod, "Int");
  register_element_type_map_array<bool>(mod, "bool");
}
} // namespace akantu
