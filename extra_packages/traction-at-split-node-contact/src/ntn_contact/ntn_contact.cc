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
#include "ntn_contact.hh"
#include "dumper_nodal_field.hh"
#include "dumper_text.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
NTNContact::NTNContact(SolidMechanicsModel & model, const ID & id)
    : NTNBaseContact(model, id),
      masters(0, 1, 0, id + ":masters", std::numeric_limits<UInt>::quiet_NaN(),
              "masters"),
      lumped_boundary_masters(0, 1, 0, id + ":lumped_boundary_masters",
                              std::numeric_limits<Real>::quiet_NaN(),
                              "lumped_boundary_masters"),
      master_elements("master_elements", id) {
  AKANTU_DEBUG_IN();

  const Mesh & mesh = this->model.getMesh();
  Int spatial_dimension = this->model.getSpatialDimension();

  this->master_elements.initialize(mesh, _nb_component = 1,
                                   _spatial_dimension = spatial_dimension - 1);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::pairInterfaceNodes(const ElementGroup & slave_boundary,
                                    const ElementGroup & master_boundary,
                                    Int surface_normal_dir, const Mesh & mesh,
                                    Array<Idx> & pairs) {
  AKANTU_DEBUG_IN();

  pairs.resize(0);
  AKANTU_DEBUG_ASSERT(pairs.getNbComponent() == 2,
                      "Array of node pairs should have nb_component = 2,"
                          << " but has nb_component = "
                          << pairs.getNbComponent());

  Int dim = mesh.getSpatialDimension();
  AKANTU_DEBUG_ASSERT(surface_normal_dir < dim,
                      "Mesh is of " << dim << " dimensions"
                                    << " and cannot have direction "
                                    << surface_normal_dir
                                    << " for surface normal");

  // offset for projection computation
  Vector<Int> offset(dim - 1);
  for (Int i = 0, j = 0; i < dim; ++i) {
    if (surface_normal_dir != i) {
      offset(j) = i;
      ++j;
    }
  }

  // find projected node coordinates
  const auto & coordinates = mesh.getNodes();

  // find slave nodes
  Array<Real> proj_slave_coord(slave_boundary.getNbNodes(), dim - 1, 0.);
  Array<Idx> slave_nodes(slave_boundary.getNbNodes());
  Int n(0);
  for (auto && slave_node : slave_boundary.getNodeGroup().getNodes()) {
    for (Int d = 0; d < dim - 1; ++d) {
      proj_slave_coord(n, d) = coordinates(slave_node, offset[d]);
      slave_nodes(n) = slave_node;
    }
    ++n;
  }

  // find master nodes
  Array<Real> proj_master_coord(master_boundary.getNbNodes(), dim - 1, 0.);
  Array<Idx> master_nodes(master_boundary.getNbNodes());
  n = 0;
  for (auto && master_node : master_boundary.getNodeGroup().getNodes()) {
    for (Int d = 0; d < dim - 1; ++d) {
      proj_master_coord(n, d) = coordinates(master_node, offset[d]);
      master_nodes(n) = master_node;
    }
    ++n;
  }

  // find minimum distance between slave nodes to define tolerance
  Real min_dist = std::numeric_limits<Real>::max();
  for (Int i = 0; i < proj_slave_coord.size(); ++i) {
    for (Int j = i + 1; j < proj_slave_coord.size(); ++j) {
      Real dist = 0.;
      for (Int d = 0; d < dim - 1; ++d) {
        dist += (proj_slave_coord(i, d) - proj_slave_coord(j, d)) *
                (proj_slave_coord(i, d) - proj_slave_coord(j, d));
      }
      if (dist < min_dist) {
        min_dist = dist;
      }
    }
  }
  min_dist = std::sqrt(min_dist);
  Real local_tol = 0.1 * min_dist;

  // find master slave node pairs
  for (Int i = 0; i < proj_slave_coord.size(); ++i) {
    for (Int j = 0; j < proj_master_coord.size(); ++j) {
      Real dist = 0.;
      for (Int d = 0; d < dim - 1; ++d) {
        dist += (proj_slave_coord(i, d) - proj_master_coord(j, d)) *
                (proj_slave_coord(i, d) - proj_master_coord(j, d));
      }
      dist = std::sqrt(dist);
      if (dist < local_tol) { // it is a pair
        pairs.push_back(Vector<Idx>{slave_nodes(i), master_nodes(j)});
        continue; // found master do not need to search further for this slave
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::addSurfacePair(const ID & slave, const ID & master,
                                Int surface_normal_dir) {
  AKANTU_DEBUG_IN();

  const auto & mesh = this->model.getMesh();
  const auto & slave_boundary = mesh.getElementGroup(slave);
  const auto & master_boundary = mesh.getElementGroup(master);
  this->contact_surfaces.insert(&slave_boundary);
  this->contact_surfaces.insert(&master_boundary);

  Array<Idx> pairs(0, 2);
  NTNContact::pairInterfaceNodes(slave_boundary, master_boundary,
                                 surface_normal_dir, this->model.getMesh(),
                                 pairs);

  // eliminate pairs which contain a pbc slave node
  Array<Idx> pairs_no_PBC_slaves(0, 2);

  for (auto && pair : make_view<2>(pairs)) {
    if (not mesh.isPeriodicSlave(pair(0)) and
        not mesh.isPeriodicSlave(pair(1))) {
      pairs_no_PBC_slaves.push_back(pair);
    }
  }

  this->addNodePairs(pairs_no_PBC_slaves);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::addNodePairs(const Array<Idx> & pairs) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_ASSERT(pairs.getNbComponent() == 2,
                      "Array of node pairs should have nb_component = 2,"
                          << " but has nb_component = "
                          << pairs.getNbComponent());
  for (auto && pair : make_view<2>(pairs)) {
    this->addSplitNode(pair(0), pair(1));
  }

  // synchronize with depending nodes
  findBoundaryElements(this->slaves.getArray(), this->slave_elements);
  findBoundaryElements(this->masters.getArray(), this->master_elements);

  updateInternalData();
  syncArrays(_added);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::getNodePairs(Array<Idx> & pairs) const {
  AKANTU_DEBUG_IN();

  pairs.resize(0);
  AKANTU_DEBUG_ASSERT(pairs.getNbComponent() == 2,
                      "Array of node pairs should have nb_component = 2,"
                          << " but has nb_component = "
                          << pairs.getNbComponent());

  auto nb_pairs = this->getNbContactNodes();
  for (Int n = 0; n < nb_pairs; ++n) {
    pairs.push_back(Vector<Idx>{this->slaves(n), this->masters(n)});
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::addSplitNode(Idx slave, Idx master) {
  AKANTU_DEBUG_IN();

  NTNBaseContact::addSplitNode(slave);

  this->masters.push_back(master);
  this->lumped_boundary_masters.push_back(
      std::numeric_limits<Real>::quiet_NaN());

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/*
  This function only works for surface elements with one quad point. For
  surface elements with more quad points, it computes still, but the result
  might not be what you are looking for.
 */
void NTNContact::updateNormals() {
  AKANTU_DEBUG_IN();

  // set normals to zero
  this->normals.zero();

  // contact information
  auto dim = this->model.getSpatialDimension();
  this->synch_registry->synchronize(
      SynchronizationTag::_cf_nodal); // synchronize current pos
  const auto & cur_pos = this->model.getCurrentPosition();

  auto & boundary_fem = this->model.getFEEngineBoundary();
  const auto & mesh = this->model.getMesh();

  for (auto ghost_type : ghost_types) {
    for (const auto & type : mesh.elementTypes(dim - 1, ghost_type)) {
      // compute the normals
      Array<Real> quad_normals(0, dim);
      boundary_fem.computeNormalsOnIntegrationPoints(cur_pos, quad_normals,
                                                     type, ghost_type);

      auto nb_quad_points =
          boundary_fem.getNbIntegrationPoints(type, ghost_type);

      // new version: compute normals only based on master elements (and not all
      // boundary elements)
      // -------------------------------------------------------------------------------------

      auto nb_nodes_per_element = mesh.getNbNodesPerElement(type);
      const auto & connectivity = mesh.getConnectivity(type, ghost_type);

      // loop over contact nodes
      for (auto & element : (this->master_elements)(type, ghost_type)) {
        for (Int q = 0; q < nb_nodes_per_element; ++q) {
          auto node = connectivity(element, q);
          auto node_index = this->masters.find(node);
          AKANTU_DEBUG_ASSERT(node_index != -1, "Could not find node "
                                                    << node
                                                    << " in the array!");

          for (Int q = 0; q < nb_quad_points; ++q) {
            // add quad normal to master normal
            for (Int d = 0; d < dim; ++d) {
              this->normals(node_index, d) +=
                  quad_normals(element * nb_quad_points + q, d);
            }
          }
        }
      }
    }
  }

  for (auto && n : make_view(this->normals, dim)) {
    n.normalize();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();

  NTNBaseContact::dumpRestart(file_name);
  this->masters.dumpRestartFile(file_name);
  this->lumped_boundary_masters.dumpRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();

  NTNBaseContact::readRestart(file_name);
  this->masters.readRestartFile(file_name);
  this->lumped_boundary_masters.readRestartFile(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::updateImpedance() {
  AKANTU_DEBUG_IN();

  Real delta_t = this->model.getTimeStep();
  AKANTU_DEBUG_ASSERT(delta_t != NAN,
                      "Time step is NAN. Have you set it already?");

  const Array<Real> & mass = this->model.getMass();
  
  for(auto && [imp,m,s,lbm,lbs]:
	zip(this->impedance, this->masters, this->slaves, this->lumped_boundary_masters, this->lumped_boundary_slaves)){
    imp = 2 / delta_t / ( ( lbm/mass(m) + lbs/mass(s)));
  }
  
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::updateLumpedBoundary() {
  AKANTU_DEBUG_IN();

  internalUpdateLumpedBoundary(this->slaves.getArray(), this->slave_elements,
                               this->lumped_boundary_slaves);

  internalUpdateLumpedBoundary(this->masters.getArray(), this->master_elements,
                               this->lumped_boundary_masters);

  AKANTU_DEBUG_OUT();
}

/* --------------------------------------------------------------------------
void NTNContact::applyContactPressure() {
  auto nb_ntn_pairs = getNbContactNodes();
  auto dim = this->model.getSpatialDimension();

  auto & residual = this->model.getInternalForce();

  for (Int n = 0; n < nb_ntn_pairs; ++n) {
    auto master = this->masters(n);
    auto slave = this->slaves(n);

    for (Int d = 0; d < dim; ++d) {
      residual(master, d) +=
          this->lumped_boundary_masters(n) * this->contact_pressure(n, d);
      residual(slave, d) -=
          this->lumped_boundary_slaves(n) * this->contact_pressure(n, d);
    }
  }
}

*/

/* -------------------------------------------------------------------------- */
void NTNContact::assembleGlobalContactPressure() {
  auto dim = this->model.getSpatialDimension();
  auto & global_contact_pressure = this->getGlobalContactPressure();
  global_contact_pressure.zero();
  
  for (auto && [mgcp, sgcp, cp, lbm, lbs] :
       zip(filter(masters, make_view(global_contact_pressure, dim)),
           filter(slaves, make_view(global_contact_pressure, dim)),
           make_view(this->contact_pressure, dim),
           this->lumped_boundary_masters, this->lumped_boundary_slaves)) {
    mgcp = lbm * cp;
    sgcp = -lbs * cp;
  }
}

/* -------------------------------------------------------------------------- */
void NTNContact::computeRelativeTangentialField(
    const Array<Real> & field, Array<Real> & rel_tang_field) const {
  // resize arrays to zero
  rel_tang_field.resize(0);

  Int dim = this->model.getSpatialDimension();

  auto it_field = field.begin(dim);
  auto it_normal = this->normals.getArray().begin(dim);

  Vector<Real> rfv(dim);
  Vector<Real> np_rfv(dim);

  auto nb_contact_nodes = this->slaves.size();
  for (Int n = 0; n < nb_contact_nodes; ++n) {
    // nodes
    auto slave = this->slaves(n);
    auto master = this->masters(n);

    // relative field vector (slave - master)
    rfv = it_field[slave] - it_field[master];

    // normal projection of relative field
    const Vector<Real> normal_v = it_normal[n];
    np_rfv = normal_v * rfv.dot(normal_v);

    // subract normal projection from relative field to get the tangential
    // projection
    rfv -= np_rfv;
    rel_tang_field.push_back(rfv);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::computeRelativeNormalField(
    const Array<Real> & field, Array<Real> & rel_normal_field) const {
  AKANTU_DEBUG_IN();

  // resize arrays to zero
  rel_normal_field.resize(0);

  Int dim = this->model.getSpatialDimension();

  auto it_field = make_view(field, dim).begin();
  auto it_normal = make_view(this->normals.getArray(), dim).begin();

  Vector<Real> rfv(dim);
  
  auto nb_contact_nodes = this->getNbContactNodes();
  for (Int n = 0; n < nb_contact_nodes; ++n) {
    // nodes
    auto slave = this->slaves(n);
    auto master = this->masters(n);

    // relative field vector (slave - master)
    rfv = it_field[slave] - it_field[master];

    // length of normal projection of relative field
    rel_normal_field.push_back(rfv.dot(it_normal[n]));
  }
}

/* -------------------------------------------------------------------------- */
Int NTNContact::getNodeIndex(Idx node) const {
  auto slave_i = NTNBaseContact::getNodeIndex(node);
  auto master_i = this->masters.find(node);

  return std::max(slave_i, master_i);
}

/* -------------------------------------------------------------------------- */
void NTNContact::printself(std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space(AKANTU_INDENT, indent);

  stream << space << "NTNContact [" << std::endl;
  NTNBaseContact::printself(stream, indent);
  stream << space << " + masters       : " << std::endl;
  this->masters.printself(stream, indent + 2);
  stream << space << " + lumped_boundary_mastres : " << std::endl;
  this->lumped_boundary_masters.printself(stream, indent + 2);

  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::syncArrays(SyncChoice sync_choice) {
  AKANTU_DEBUG_IN();

  NTNBaseContact::syncArrays(sync_choice);

  this->masters.syncElements(sync_choice);
  this->lumped_boundary_masters.syncElements(sync_choice);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNContact::addDumpFieldToDumper(const std::string & dumper_name,
                                      const std::string & field_id) {
  AKANTU_DEBUG_IN();

  /*
  const Array<UInt> & nodal_filter = this->slaves.getArray();

#define ADD_FIELD(field_id, field, type)				\
  internalAddDumpFieldToDumper(dumper_name,				\
                   field_id,				\
                   new DumperIOHelper::NodalField< type, true, \
                                   Array<type>, \
                                   Array<UInt> >(field, 0, 0, &nodal_filter))
  */

  if (field_id == "lumped_boundary_master") {
    internalAddDumpFieldToDumper(dumper_name, field_id,
                                 std::make_unique<dumpers::NodalField<Real>>(
                                     this->lumped_boundary_masters.getArray()));
  } else {
    NTNBaseContact::addDumpFieldToDumper(dumper_name, field_id);
  }

  /*
#undef ADD_FIELD
  */

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
