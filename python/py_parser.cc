/**
 * Copyright (©) 2020-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include <aka_common.hh>
#include <aka_iterators.hh>
#include <parameter_registry.hh>
#include <parsable.hh>
#include <parser.hh>
/* -------------------------------------------------------------------------- */
#include <map>
#include <pybind11/attr.h>
#include <pybind11/detail/common.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <variant>
/* -------------------------------------------------------------------------- */
namespace py = pybind11;
/* -------------------------------------------------------------------------- */

namespace akantu {
std::map<void *, std::map<std::string, void *>> map_params;

void register_parser(py::module & mod) {
  py::enum_<ParameterAccessType>(mod, "ParameterAccessType", py::arithmetic())
      .value("_pat_internal", _pat_internal)
      .value("_pat_writable", _pat_writable)
      .value("_pat_readable", _pat_readable)
      .value("_pat_modifiable", _pat_modifiable)
      .value("_pat_parsable", _pat_parsable)
      .value("_pat_parsmod", _pat_parsmod)
      .export_values();

  using ParsableVariant =
      std::variant<Real, Int, bool, std::string, Matrix<Real>, Vector<Real>>;

  py::class_<ParameterRegistry>(mod, "ParameterRegistry",
                                py::multiple_inheritance())
      .def("registerParamReal",
           [](ParameterRegistry & self, const std::string & name, int type,
              const std::string & description) {
             Real * p = new Real; // NOLINT
             map_params[&self][name] = p;
             self.registerParam<Real>(name, *p, ParameterAccessType(type),
                                      description);
           })
      .def("registerParamReal",
           [](ParameterRegistry & self, const Real & _default,
              const std::string & name, int type,
              const std::string & description) {
             Real * p = new Real; // NOLINT
             map_params[&self][name] = p;
             self.registerParam<Real>(name, *p, _default,
                                      ParameterAccessType(type), description);
           })

      .def("set",
           [](ParameterRegistry & self, std::string param,
              ParsableVariant value) {
             std::visit(
                 [&param, &self](auto && arg) {
                   std::any val = arg;
                   self.set(param, val);
                 },
                 value);
           })

      .def("get",
           [](ParameterRegistry & self, const std::string & param) {
             ParsableVariant ret;
             auto && parameter = self.get(param);
             if (parameter.isOfType<Real>()) {
               ret = parameter.get<Real>();
             } else if (parameter.isOfType<Int>()) {
               ret = parameter.get<Int>();
             } else if (parameter.isOfType<bool>()) {
               ret = parameter.get<bool>();
             } else if (parameter.isOfType<std::string>()) {
               ret = parameter.get<std::string>();
             } else if (parameter.isOfType<Matrix<Real>>()) {
               ret = parameter.get<Matrix<Real>>();
             } else if (parameter.isOfType<Vector<Real>>()) {
               ret = parameter.get<Vector<Real>>();
             } else {
               AKANTU_EXCEPTION(
                   "Paramter not of a known type or with wrong access rights");
             }

             return ret;
           })

      .def("getReal",
           [](ParameterRegistry & self, const std::string & name) -> Real {
             return self.get(name);
           })

      .def("getBool",
           [](ParameterRegistry & self, const std::string & name) -> bool {
             return self.get(name);
           })
      .def("getString",
           [](ParameterRegistry & self,
              const std::string & name) -> std::string {
             std::string tmp = self.get(name);
             return tmp;
           })
      .def("getInt",
           [](ParameterRegistry & self, const std::string & name) -> Int {
             return self.get(name);
           })

      .def(
          "getMatrix",
          [](ParameterRegistry & self, const std::string & name) {
            const Matrix<Real> & res =
                static_cast<const Matrix<Real> &>(self.get(name));
            return res;
          },
          py::return_value_policy::copy)

      .def("setReal", [](ParameterRegistry & self, const std::string & name,
                         const Real value) { self.set(name, value); })
      .def("setBool", [](ParameterRegistry & self, const std::string & name,
                         const bool value) { self.set(name, value); })
      .def("setString",
           [](ParameterRegistry & self, const std::string & name,
              const std::string & value) { self.set(name, value); })

      .def("setInt", [](ParameterRegistry & self, const std::string & name,
                        const Int value) { self.set(name, value); })
      .def("setRandom",
           [](ParameterRegistry & self, const std::string & name, Real value,
              const RandomDistributionType & distribution, Real a, Real b) {
             auto param = make_random_parameter(value, distribution, a, b);
             self.set(name, value);
           });

  py::class_<Parsable, ParameterRegistry>(mod, "Parsable")
      .def(py::init<const ParserType &, const ID &>())
      .def("parseSection", &Parsable::parseSection);

  py::class_<ParserParameter>(mod, "ParserParameter")
      .def("getName", &ParserParameter::getName)
      .def("getValue", &ParserParameter::getValue)
      .def(
          "asInt",
          [](ParserParameter & self) {
            Int i = self;
            return i;
          },
          py::return_value_policy::copy)
      .def(
          "asReal",
          [](ParserParameter & self) {
            Real r = self;
            return r;
          },
          py::return_value_policy::copy)
      .def(
          "asBool",
          [](ParserParameter & self) {
            bool b = self;
            return b;
          },
          py::return_value_policy::copy)
      .def(
          "asVector",
          [](ParserParameter & self) {
            Vector<Real> v = self;
            return v;
          },
          py::return_value_policy::copy)
      .def(
          "asMatrix",
          [](ParserParameter & self) {
            Matrix<Real> m = self;
            return m;
          },
          py::return_value_policy::reference);

  py::class_<ParserSection>(mod, "ParserSection")
      .def(py::init<const std::string &, ParserType>())
      .def(py::init<const std::string &, ParserType, const std::string &,
                    const ParserSection &>())

      .def(
          "getParameter",
          [](ParserSection & self, const std::string & name) {
            return self.getParameter(name);
          },
          py::return_value_policy::reference_internal)
      .def(
          "getSubSections",
          [](ParserSection & self, ParserType type, int id = 0) {
            auto sections = self.getSubSections(type);
            int i = 0;
            for (auto && section : sections) {
              if (i == id) {
                return section;
              }
              ++i;
            }
            AKANTU_EXCEPTION("Sub section not found");
          },
          py::arg("type"), py::arg("id") = 0,
          py::return_value_policy::reference);

  py::class_<Parser, ParserSection>(mod, "Parser").def(py::init());

  mod.def(
      "parseInput",
      [](const std::string & input_file) {
        getStaticParser().parse(input_file);
      },
      "Parse an Akantu input file");

  mod.def(
      "getStaticParser", []() { return getStaticParser(); },
      "Get the global parser", py::return_value_policy::reference);
}
} // namespace akantu
