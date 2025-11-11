/**
 * Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
// simtools
#include "ntrf_contact.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
NTRFContact::NTRFContact(SolidMechanicsModel & model, const ID & id)
    : NTNBaseContact(model, id), reference_point(model.getSpatialDimension()),
      normal(model.getSpatialDimension()) {
  AKANTU_DEBUG_IN();

  is_ntn_contact = false;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::setReferencePoint(Real x, Real y, Real z) {
  AKANTU_DEBUG_IN();

  Real coord[3];
  coord[0] = x;
  coord[1] = y;
  coord[2] = z;

  Int dim = this->model.getSpatialDimension();
  for (Int d = 0; d < dim; ++d) {
    this->reference_point(d) = coord[d];
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::setNormal(Real x, Real y, Real z) {
  AKANTU_DEBUG_IN();

  Int dim = this->model.getSpatialDimension();

  Real coord[3];
  coord[0] = x;
  coord[1] = y;
  coord[2] = z;

  for (Int d = 0; d < dim; ++d) {
    this->normal(d) = coord[d];
  }

  this->normal.normalize();

  this->updateNormals();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::addSurface(const ID & surf) {
  AKANTU_DEBUG_IN();

  const Mesh & mesh_ref = this->model.getMesh();

  try {
    const ElementGroup & boundary = mesh_ref.getElementGroup(surf);
    this->contact_surfaces.insert(&boundary);

    // find slave nodes
    for (auto && node : boundary.getNodeGroup().getNodes()) {
      if (not mesh_ref.isPeriodicSlave(node)) {
        this->addSplitNode(node);
      }
    }
  } catch (debug::Exception & e) {
    AKANTU_DEBUG_INFO("NTRFContact addSurface did not found subboundary "
                      << surf
                      << " and ignored it. Other procs might have it :)");
  }

  // synchronize with depending nodes
  findBoundaryElements(this->slaves.getArray(), this->slave_elements);
  updateInternalData();
  syncArrays(_added);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::addNodes(Array<UInt> & nodes) {
  AKANTU_DEBUG_IN();

  UInt nb_nodes = nodes.size();
  UInt nb_compo = nodes.getNbComponent();
  for (UInt n = 0; n < nb_nodes; ++n) {
    for (UInt c = 0; c < nb_compo; ++c) {
      this->addSplitNode(nodes(n, c));
    }
  }

  // synchronize with depending nodes
  findBoundaryElements(this->slaves.getArray(), this->slave_elements);
  updateInternalData();
  syncArrays(_added);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::updateNormals() {
  AKANTU_DEBUG_IN();

  // normal is the same for all slaves
  this->normals.set(this->normal);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::updateImpedance() {
  AKANTU_DEBUG_IN();

  auto nb_contact_nodes = getNbContactNodes();
  Real delta_t = this->model.getTimeStep();
  AKANTU_DEBUG_ASSERT(delta_t != NAN,
                      "Time step is NAN. Have you set it already?");

  const Array<Real> & mass = this->model.getMass();

  for (Int n = 0; n < nb_contact_nodes; ++n) {
    auto slave = this->slaves(n);

    Real imp = this->lumped_boundary_slaves(n) / mass(slave);
    imp = 2 / delta_t / imp;
    this->impedance(n) = imp;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::computeRelativeTangentialField(
    const Array<Real> & field, Array<Real> & rel_tang_field) const {
  AKANTU_DEBUG_IN();

  // resize arrays to zero
  rel_tang_field.resize(0);

  auto dim = this->model.getSpatialDimension();

  auto it_field = field.begin(dim);
  auto it_normal = this->normals.getArray().begin(dim);

  Vector<Real> rfv(dim);
  Vector<Real> np_rfv(dim);

  auto nb_contact_nodes = this->slaves.size();
  for (Int n = 0; n < nb_contact_nodes; ++n) {
    // nodes
    auto node = this->slaves(n);

    // relative field vector
    rfv = it_field[node];

    // normal projection of relative field
    np_rfv = it_normal[n] * rfv.dot(it_normal[n]);

    // subtract normal projection from relative field to get the tangential
    // projection
    rfv -= np_rfv;
    rel_tang_field.push_back(rfv);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::computeNormalGap(Array<Real> & gap) const {
  AKANTU_DEBUG_IN();

  gap.resize(0);

  Int dim = this->model.getSpatialDimension();

  auto it_cur_pos = this->model.getCurrentPosition().begin(dim);
  auto it_normal = this->normals.getArray().begin(dim);

  Vector<Real> gap_v(dim);

  Int nb_contact_nodes = this->getNbContactNodes();
  for (Int n = 0; n < nb_contact_nodes; ++n) {
    // nodes
    auto node = this->slaves(n);

    // gap vector
    gap_v = it_cur_pos[node] - this->reference_point;

    // length of normal projection of gap vector
    gap.push_back(gap_v.dot(it_normal[n]));
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTRFContact::computeRelativeNormalField(
    const Array<Real> & field, Array<Real> & rel_normal_field) const {
  // resize arrays to zero
  rel_normal_field.resize(0);

  auto dim = this->model.getSpatialDimension();

  auto it_field = field.begin(dim);
  auto it_normal = this->normals.getArray().begin(dim);

  auto nb_contact_nodes = this->getNbContactNodes();
  for (Int n = 0; n < nb_contact_nodes; ++n) {
    // nodes
    auto node = this->slaves(n);

    rel_normal_field.push_back(it_field[node].dot(it_normal[n]));
  }
}

/* -------------------------------------------------------------------------- */
void NTRFContact::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "NTRFContact [" << std::endl;
  NTNBaseContact::printself(stream, indent);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void NTRFContact::addDumpFieldToDumper(const std::string & dumper_name,
                                       const std::string & field_id) {
  NTNBaseContact::addDumpFieldToDumper(dumper_name, field_id);
}

} // namespace akantu
