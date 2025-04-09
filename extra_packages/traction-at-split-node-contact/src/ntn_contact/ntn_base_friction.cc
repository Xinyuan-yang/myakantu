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
#include "ntn_base_friction.hh"
#include "aka_common.hh"
#include "dof_manager_default.hh"
#include "dumpable.hh"
#include "dumper_elemental_field.hh"
#include "dumper_nodal_field.hh"
#include "dumper_text.hh"
#include "non_linear_solver_lumped.hh"
#include "ntn_base_contact.hh"
#include "synchronized_array.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
#include <ostream>
#include <string>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
NTNBaseFriction::NTNBaseFriction(NTNBaseContact & contact, const ID & id)
    : Parsable(ParserType::_friction, id), contact(contact),
      is_sticking(0, 1, true, id + ":is_sticking", true, "is_sticking"),
      frictional_strength(0, 1, 0., id + ":frictional_strength", 0.,
                          "frictional_strength"),
      friction_traction(0, contact.getModel().getSpatialDimension(), 0.,
                        id + ":friction_traction", 0., "friction_traction"),
      friction_traction_old(0, contact.getModel().getSpatialDimension(), 0.,
                            id + ":friction_traction_old", 0.,
                            "friction_traction_old"),
      slip(0, 1, 0., id + ":slip", 0., "slip"),
      cumulative_slip(0, 1, 0., id + ":cumulative_slip", 0., "cumulative_slip"),
      slip_velocity(0, contact.getModel().getSpatialDimension(), 0.,
                    id + ":slip_velocity", 0., "slip_velocity"),
      slip_velocity_norm(0, 1, 0., id + ":slip_velocity_norm", 0.,
                         "slip_velocity_norm"),
      friction_state(0, 1, 0., id + ":friction_state", 0., "friction_state"),
      global_friction_traction(std::make_unique<akantu::Array<Real>>(
          0, contact.getModel().getSpatialDimension(), 0.,
          id + ":global_friction_traction")),
      global_friction_traction_old(std::make_unique<akantu::Array<Real>>(
          0, contact.getModel().getSpatialDimension(),
          id + ":global_friction_traction_old")) {

  AKANTU_DEBUG_IN();

  this->contact.registerSynchronizedArray(this->is_sticking);
  this->contact.registerSynchronizedArray(this->frictional_strength);
  this->contact.registerSynchronizedArray(this->friction_traction);
  this->contact.registerSynchronizedArray(this->friction_traction_old);
  this->contact.registerSynchronizedArray(this->slip);
  this->contact.registerSynchronizedArray(this->cumulative_slip);
  this->contact.registerSynchronizedArray(this->slip_velocity);
  this->contact.registerSynchronizedArray(this->slip_velocity_norm);
  this->contact.registerSynchronizedArray(this->friction_state);

  this->registerExternalDumper(contact.getDumper().shared_from_this(),
                               contact.getDefaultDumperName(), true);

  const auto & mesh = this->contact.getModel().getMesh();
  auto nb_nodes = mesh.getNbNodes();

  // Array<Real>
  // global_friction_traction(nb_global_nodes,contact.getModel().getSpatialDimension());

  this->global_friction_traction->resize(nb_nodes, 0.);
  this->global_friction_traction_old->resize(nb_nodes, 0.);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::updateSlip() {
  AKANTU_DEBUG_IN();

  SolidMechanicsModel & model = this->contact.getModel();
  Int dim = model.getSpatialDimension();
  Real delta_t = model.getTimeStep();

  // synchronize increment
  this->contact.getSynchronizerRegistry().synchronize(
      SynchronizationTag::_cf_incr);

  Array<Real> rel_tan_incr(0, dim);
  this->contact.computeRelativeTangentialField(model.getIncrement(),
                                               rel_tan_incr);
  auto it = rel_tan_incr.begin(dim);

  for (auto && [rti, is_sticking, slip, slip_vel_norm, slip_vel,
                cumulative_slip] :
       zip(make_view(rel_tan_incr, dim), this->is_sticking, this->slip,
           this->slip_velocity_norm, make_view(this->slip_velocity, dim),
           this->cumulative_slip)) {

    if (is_sticking) {
      slip = 0.;
      slip_vel.zero();
      slip_vel_norm = 0.;
    } else {
      slip += rti.norm();
      cumulative_slip += rti.norm();

      slip_vel = rti / delta_t;
      slip_vel_norm = slip_vel.norm();
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::computeFrictionTraction() {
  AKANTU_DEBUG_IN();

  this->computeStickTraction();
  this->computeFrictionalStrength();

  auto & model = this->contact.getModel();
  auto dim = model.getSpatialDimension();

  // get contact arrays
  const SynchronizedArray<bool> & is_in_contact =
      this->contact.getIsInContact();

  auto & traction = this->friction_traction.getArray();

  this->is_sticking.zero(); // set to not sticking

  for (auto && [n, fric_trac] : enumerate(make_view(traction, dim))) {
    // node pair is in contact
    if (is_in_contact(n)) {
      // check if it is larger than frictional strength
      auto abs_fric = fric_trac.norm();

      if (abs_fric != 0.) {
        auto alpha = this->frictional_strength(n) / abs_fric;

        // larger -> sliding
        if (alpha < 1.) {
          fric_trac *= alpha;
        } else {
          this->is_sticking(n) = true;
        }
      } else {
        // frictional traction is already zero
        this->is_sticking(n) = true;
      }
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::computeStickTraction() {
  auto & model = this->contact.getModel();
  auto dim = model.getSpatialDimension();
  auto delta_t = model.getTimeStep();

  auto nb_contact_nodes = this->contact.getNbContactNodes();

  // get contact arrays
  const auto & impedance = this->contact.getImpedance();
  const auto & is_in_contact = this->contact.getIsInContact();

  Array<Real> acceleration(0, dim);
  this->contact.computeAcceleration(acceleration);

  // compute relative normal fields of velocity and acceleration
  Array<Real> r_velo(0, dim);
  Array<Real> r_acce(0, dim);
  Array<Real> r_old_acce(0, dim);
  this->contact.computeRelativeTangentialField(model.getVelocity(), r_velo);
  this->contact.computeRelativeTangentialField(acceleration, r_acce);
  this->contact.computeRelativeTangentialField(model.getAcceleration(),
                                               r_old_acce);

  AKANTU_DEBUG_ASSERT(r_velo.size() == nb_contact_nodes,
                      "computeRelativeNormalField does not give back arrays "
                          << "size == nb_contact_nodes. nb_contact_nodes = "
                          << nb_contact_nodes
                          << " | array size = " << r_velo.size());

  // compute tangential gap_dot array for all nodes
  Array<Real> gap_dot(nb_contact_nodes, dim);
  for (auto && [gap_dot, r_velo, r_acce, r_old_acce] :
       zip(make_view(gap_dot), make_view(r_velo), make_view(r_acce),
           make_view(r_old_acce))) {
    gap_dot = r_velo + delta_t * r_acce - 1. / 2. * delta_t * r_old_acce;
  }

  // compute friction traction to stop sliding
  auto & traction = this->friction_traction.getArray();
  for (auto && [fric_trac, in_contact, impedance, gap_dot] :
       zip(make_view(traction, dim), is_in_contact, impedance,
           make_view(gap_dot, dim))) {
    fric_trac = (in_contact ? 1. : 0.) * impedance * gap_dot / 2;
  }
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::assembleGlobalFrictionTraction() {
  auto & global_friction_traction = this->getGlobalFrictionTraction();
  global_friction_traction.zero();

  auto dim = global_friction_traction.getNbComponent();

  const auto & lumped_boundary_slaves = this->contact.getLumpedBoundarySlaves();
  const auto & slaves = this->contact.getSlaves();

  for (auto && [gft, ft, lbs] :
       zip(filter(slaves, make_view(global_friction_traction, dim)),
           make_view(friction_traction, dim), lumped_boundary_slaves)) {
    gft = -lbs * ft;
  }

  const auto & mesh = this->contact.getModel().getMesh();
  mesh.makeConsistentForPeriodicity(global_friction_traction,
                                    ReduceOperation::_id);
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::registerSynchronizedArray(SynchronizedArrayBase & array) {
  this->frictional_strength.registerDependingArray(array);
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::dumpRestart(const std::string & file_name) const {
  this->is_sticking.dumpRestartFile(file_name);
  this->frictional_strength.dumpRestartFile(file_name);
  this->friction_traction.dumpRestartFile(file_name);
  this->friction_traction_old.dumpRestartFile(file_name);
  this->slip.dumpRestartFile(file_name);
  this->cumulative_slip.dumpRestartFile(file_name);
  this->slip_velocity.dumpRestartFile(file_name);
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::readRestart(const std::string & file_name) {
  this->is_sticking.readRestartFile(file_name);
  this->frictional_strength.readRestartFile(file_name);
  this->friction_traction.readRestartFile(file_name);
  this->friction_traction_old.readRestartFile(file_name);
  this->cumulative_slip.readRestartFile(file_name);
  this->slip_velocity.readRestartFile(file_name);
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::setParam(const std::string & name, UInt node,
                               Real value) {
  auto & array = this->get(name).get<SynchronizedArray<Real>>();
  Int index = this->contact.getNodeIndex(node);
  if (index < 0) {
    AKANTU_DEBUG_WARNING("Node "
                         << node << " is not a contact node. "
                         << "Therefore, cannot set interface parameter!!");
  } else {
    array(index) = value; // put value
  }
}

/* -------------------------------------------------------------------------- */
UInt NTNBaseFriction::getNbStickingNodes() const {
  Int nb_stick = 0;

  auto nb_nodes = this->contact.getNbContactNodes();
  const auto & nodes = this->contact.getSlaves();
  const auto & is_in_contact = this->contact.getIsInContact();

  const auto & mesh = this->contact.getModel().getMesh();

  for (Int n = 0; n < nb_nodes; ++n) {
    auto is_local_node = mesh.isLocalOrMasterNode(nodes(n));
    auto is_pbc_slave_node = mesh.isPeriodicSlave(nodes(n));
    if (is_local_node && !is_pbc_slave_node && is_in_contact(n) &&
        this->is_sticking(n)) {
      nb_stick++;
    }
  }

  mesh.getCommunicator().allReduce(nb_stick, SynchronizerOperation::_sum);
  return nb_stick;
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "NTNBaseFriction [" << std::endl;
  Parsable::printself(stream, indent);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::addDumpFieldToDumper(const std::string & dumper_name,
                                           const std::string & field_id) {
  AKANTU_DEBUG_IN();

  if (field_id == "is_sticking") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<bool>>(
            this->is_sticking.getArray()));
  } else if (field_id == "frictional_strength") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(
            this->frictional_strength.getArray()));
  } else if (field_id == "friction_traction") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(
            this->friction_traction.getArray()));
  } else if (field_id == "friction_traction_old") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(
            this->friction_traction_old.getArray()));
  } else if (field_id == "global_friction_traction_old") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(
            *this->global_friction_traction_old));
  } else if (field_id == "global_friction_traction") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(
            *this->global_friction_traction));
  } else if (field_id == "slip") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(this->slip.getArray()));
  } else if (field_id == "cumulative_slip") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(
            this->cumulative_slip.getArray()));
  } else if (field_id == "slip_velocity") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(
            this->slip_velocity.getArray()));
  } else if (field_id == "slip_velocity_norm") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(
            this->slip_velocity_norm.getArray()));
  } else {
    this->contact.addDumpFieldToDumper(dumper_name, field_id);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::updateFrictionState(UInt node, Real value) {
  friction_state(node) = value;
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::addNoiseToStateVariable(std::map<Real, Real> & noise) {
  AKANTU_DEBUG_IN();

  this->addNoiseToStateVariable(noise);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::savePreviousState() {
  global_friction_traction_old->copy(*global_friction_traction);
}

/* -------------------------------------------------------------------------- */
void NTNBaseFriction::restorePreviousState() {
  global_friction_traction->copy(*global_friction_traction_old);
}

} // namespace akantu
