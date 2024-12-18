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
#include <dumper_iohelper_paraview.hh>
#include <mesh.hh>
/* -------------------------------------------------------------------------- */
#include <dumpable_inline_impl.hh>
/* -------------------------------------------------------------------------- */
#include <pybind11/pybind11.h>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {

std::vector<detail::ArrayProxy<Real>> tmp_array;

void register_dumpable(py::module & mod) {
  /* ------------------------------------------------------------------------ */
  py::class_<Dumpable>(mod, "Dumpable")
      .def("registerDumperParaview", &Dumpable::registerDumper<DumperParaview>,
           py::arg("dumper_name"), py::arg("file_name"),
           py::arg("is_default") = false)
      .def("addDumpMeshToDumper", &Dumpable::addDumpMeshToDumper,
           py::arg("dumper_name"), py::arg("mesh"), py::arg("dimension"),
           py::arg("ghost_type") = _not_ghost,
           py::arg("element_kind") = _ek_regular)
      .def("addDumpMesh", &Dumpable::addDumpMesh, py::arg("mesh"),
           py::arg("dimension"), py::arg("ghost_type") = _not_ghost,
           py::arg("element_kind") = _ek_regular)
      .def("addDumpField", &Dumpable::addDumpField, py::arg("field_id"))
      .def("addDumpFieldToDumper", &Dumpable::addDumpFieldToDumper,
           py::arg("dumper_name"), py::arg("field_id"))
      .def(
          "addDumpFieldExternal",
          [](Dumpable & _this, const std::string & field_id,
             std::shared_ptr<dumpers::Field> field) {
            return _this.addDumpFieldExternal(field_id, field);
          },
          py::arg("field_id"), py::arg("field"))
      .def(
          "addDumpFieldExternal",
          [](Dumpable & _this, const std::string & field_id,
             Array<Real> & field) {
            auto & tmp = dynamic_cast<detail::ArrayProxy<Real> &>(field);
            tmp_array.push_back(tmp);
            return _this.addDumpFieldExternal(field_id, tmp_array.back());
          },
          py::arg("field_id"), py::arg("field"))
      .def(
          "addDumpFieldExternalToDumper",
          [](Dumpable & _this, const std::string & dumper_name,
             const std::string & field_id,
             std::shared_ptr<dumpers::Field> field) {
            return _this.addDumpFieldExternalToDumper(dumper_name, field_id,
                                                      field);
          },
          py::arg("dumper_name"), py::arg("field_id"), py::arg("field"))

      .def("dump", [](Dumpable & self) { self.dump(); })
      .def(
          "dump", [](Dumpable & self, UInt step) { self.dump(step); },
          py::arg("step"))
      .def(
          "dump",
          [](Dumpable & self, Real time, UInt step) { self.dump(time, step); },
          py::arg("time"), py::arg("step"))
      .def(
          "dump",
          [](Dumpable & self, const std::string & dumper) {
            self.dump(dumper);
          },
          py::arg("dumper_name"))
      .def(
          "dump",
          [](Dumpable & self, const std::string & dumper, UInt step) {
            self.dump(dumper, step);
          },
          py::arg("dumper_name"), py::arg("step"))
      .def(
          "dump",
          [](Dumpable & self, const std::string & dumper, Real time,
             UInt step) { self.dump(dumper, time, step); },
          py::arg("dumper_name"), py::arg("time"), py::arg("step"));
}

/* -------------------------------------------------------------------------- */
PYBIND11_MODULE(dumper_module, mod) {
  mod.attr("__name__") = "py11_akantu.dumper";

  /* ------------------------------------------------------------------------ */
  py::class_<dumpers::Field, std::shared_ptr<dumpers::Field>>(mod, "Field");

  /* ------------------------------------------------------------------------ */
  py::class_<dumpers::ElementalField<UInt>, dumpers::Field,
             std::shared_ptr<dumpers::ElementalField<UInt>>>(
      mod, "ElementalFieldUInt", py::multiple_inheritance())
      .def(py::init<dumpers::ElementalField<UInt>::field_type &, UInt,
                    GhostType, ElementKind>(),
           py::arg("field"), py::arg("spatial_dimension") = _all_dimensions,
           py::arg("ghost_type") = _not_ghost,
           py::arg("element_kind") = _ek_not_defined);
}

} // namespace akantu
