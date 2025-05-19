/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "dumper_text.hh"
#include "communicator.hh"
#include "dumper_compute.hh"
#include "dumper_homogenizing_field.hh"
#include "dumper_nodal_field.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include <io_helper.hh>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
DumperText::DumperText(const std::string & basename,
                       iohelper::TextDumpMode mode, bool parallel) {
  this->dumper = std::make_unique<iohelper::DumperText>(mode);
  this->setBaseName(basename);

  this->setParallelContext(parallel);
}

/* -------------------------------------------------------------------------- */
void DumperText::registerMesh(const Mesh & mesh, Int /*spatial_dimension*/,
                              GhostType /*ghost_type*/,
                              ElementKind /*element_kind*/) {
  registerField("position",
                std::make_shared<dumpers::NodalField<Real>>(mesh.getNodes()));
  // in parallel we need node type
  auto nb_proc = mesh.getCommunicator().getNbProc();
  if (nb_proc > 1) {
    auto func = std::make_unique<dumpers::ComputeIntFromEnum<ContactState>>();
    std::shared_ptr<dumpers::Field> field =
        std::make_shared<dumpers::NodalField<NodeFlag>>(mesh.getNodesFlags());
    field =
        dumpers::FieldComputeProxy::createFieldCompute(field, std::move(func));
    registerField("nodes_type", field);
  }
}

/* -------------------------------------------------------------------------- */
void DumperText::registerFilteredMesh(
    const Mesh & mesh, const ElementTypeMapArray<Idx> & /*elements_filter*/,
    const Array<Idx> & nodes_filter, Idx /*spatial_dimension*/,
    GhostType /*ghost_type*/, ElementKind /*element_kind*/) {
  registerField("position", std::make_shared<dumpers::NodalField<Real, true>>(
                                mesh.getNodes(), 0, 0, &nodes_filter));

  // in parallel we need node type
  auto nb_proc = mesh.getCommunicator().getNbProc();
  if (nb_proc > 1) {
    auto func = std::make_unique<dumpers::ComputeIntFromEnum<NodeFlag>>();
    std::shared_ptr<dumpers::Field> field =
        std::make_shared<dumpers::NodalField<NodeFlag, true>>(
            mesh.getNodesFlags(), 0, 0, &nodes_filter);
    field =
        dumpers::FieldComputeProxy::createFieldCompute(field, std::move(func));
    registerField("nodes_type", field);
  }
}

/* -------------------------------------------------------------------------- */
void DumperText::setBaseName(const std::string & basename) {
  DumperIOHelper::setBaseName(basename);
  static_cast<iohelper::DumperText *>(this->dumper.get())
      ->setDataSubDirectory(this->filename + "-DataFiles");
}

/* -------------------------------------------------------------------------- */
void DumperText::setPrecision(Int prec) {
  static_cast<iohelper::DumperText *>(this->dumper.get())->setPrecision(prec);
}

} // namespace akantu
