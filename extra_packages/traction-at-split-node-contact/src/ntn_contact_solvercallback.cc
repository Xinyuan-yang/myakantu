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
#include "ntn_contact_solvercallback.hh"
#include "aka_common.hh"
#include "boundary_condition.hh"
#include "mesh.hh"
#include "ntn_base_contact.hh"
#include "ntn_contact.hh"
#include "ntn_initiation_function.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */

namespace akantu {

// NTNContactSolverCallback::NTNContactSolverCallback(
//   SolidMechanicsModel & solid, NTNBaseContact & base_contact,
//   NTNBaseFriction & friction)
//   : InterceptSolverCallback(solid), solid(solid), base_contact(base_contact),
//     friction(friction) {}

/* -------------------------------------------------------------------------- */
NTNContactSolverCallback::NTNContactSolverCallback(SolidMechanicsModel & solid,
                                                   const ID & slave,
                                                   const ID & master,
                                                   Int surface_normal_dir,
                                                   Real time_step_factor)
    : InterceptSolverCallback(solid), solid(solid) {
  AKANTU_DEBUG_IN();

  // Create note to node contact
  contact = std::make_shared<NTNContact>(solid);

  // Init explicit solid mechanics model
  solid.initFull(SolidMechanicsModelOptions(_explicit_lumped_mass));

  // Necessary for node to node contact !
  solid.assembleMassLumped();

  contact->addSurfacePair(slave, master, surface_normal_dir);

  contact->updateNormals();
  contact->updateLumpedBoundary();

  Real stable_time_step = solid.getStableTimeStep();
  Real time_step = stable_time_step * time_step_factor;
  solid.setTimeStep(time_step);
  contact->updateImpedance();

  friction = initializeNTNFriction(*contact);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
auto NTNContactSolverCallback::getExternalWorkIncrement() const
    -> Vector<Real, 2> {
  Vector<Real, 2> work{0., 0.};
  const auto & contact_pressure = contact->getGlobalContactPressure();
  const auto & friction_traction = friction->getGlobalFrictionTraction();
  const auto & mesh = solid.getMesh();
  auto spatial_dimension = mesh.getSpatialDimension();

  for (auto && [n, inc, f_cont, f_fric] :
       enumerate(make_view(solid.getIncrement(), spatial_dimension),
                 make_view(contact_pressure, spatial_dimension),
                 make_view(friction_traction, spatial_dimension))) {
    auto is_local_node = mesh.isLocalOrMasterNode(n);
    auto is_not_pbc_slave_node = not mesh.isPeriodicSlave(n);
    Real count_node = is_local_node && is_not_pbc_slave_node ? 1. : 0.;

    work[0] -= count_node * f_cont.dot(inc);
    work[1] -= count_node * f_fric.dot(inc);
  }

  return work;
}

/* -------------------------------------------------------------------------- */
auto NTNContactSolverCallback::getExternalWork() const -> Vector<Real, 2> {
  Vector<Real, 2> work = external_work;
  solid.getMesh().getCommunicator().allReduce(work,
                                              SynchronizerOperation::_sum);
  return work;
}

/* -------------------------------------------------------------------------- */
void NTNContactSolverCallback::assembleResidual() {

  auto && dof_manager = solid.getDOFManager();

  auto & contact_pressure = contact->getGlobalContactPressure();
  auto & friction_traction = friction->getGlobalFrictionTraction();

  if (external_work_release != solid.getDisplacementRelease()) {
    external_work += this->getExternalWorkIncrement();
    external_work_release = solid.getDisplacementRelease();
  }

  // Assemble the residual of the solid mechanics problem
  // auto & internal_force = solid.getInternalForce();
  // auto & external_force = solid.getExternalForce();
  // solid.assembleInternalForces();
  // dof_manager.assembleToResidual("displacement", external_force, 1);
  // dof_manager.assembleToResidual("displacement", internal_force, 1);
  aka::as_type<SolverCallback>(solid).assembleResidual();

  // Compute the contact - friction part of the problem and add to the residual
  contact->computeContactPressure();
  friction->updateSlip();
  friction->computeFrictionTraction();

  // Assemble the global array for the contact pressure and the friction
  // traction
  contact->assembleGlobalContactPressure();
  friction->assembleGlobalFrictionTraction();

  // The two functions contact.applyContactPressure and
  // friction.applyFrictionTraction() are replaced by assembling the forces to
  // the residual directly
  dof_manager.assembleToResidual("displacement", contact_pressure, 1);
  dof_manager.assembleToResidual("displacement", friction_traction, 1);
}

/* -------------------------------------------------------------------------- */
void NTNContactSolverCallback::beforeSolveStep() {
  aka::as_type<SolverCallback>(solid).beforeSolveStep();

  this->friction->savePreviousState();
  this->contact->savePreviousState();
}

/* -------------------------------------------------------------------------- */
void NTNContactSolverCallback::afterSolveStep(bool converged) {
  aka::as_type<SolverCallback>(solid).afterSolveStep(converged);

  if (not converged) {
    this->friction->restorePreviousState();
    this->contact->restorePreviousState();
    return;
  }
}

} // namespace akantu
