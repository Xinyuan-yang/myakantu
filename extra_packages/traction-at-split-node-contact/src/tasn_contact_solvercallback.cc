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
#include "tasn_contact_solvercallback.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

TasnContactSolverCallback::TasnContactSolverCallback(
    SolidMechanicsModel & solid, NTNBaseContact & contact,
    NTNBaseFriction & friction)
    : InterceptSolverCallback(solid), solid(solid), contact(contact),
      friction(friction) {}

/* -------------------------------------------------------------------------- */
void TasnContactSolverCallback::assembleResidual() {

  // Assemble the residual of the solid mechanics problem
  solver_callback.assembleResidual();

  // Compute the contact - friction part of the problem and add to the residual
  contact.computeContactPressure();
  friction.updateSlip();
  friction.computeFrictionTraction();

  // Assemble the global array for the contact pressure and the friction
  // traction
  contact.assembleGlobalContactPressure();
  friction.assembleGlobalFrictionTraction();

  auto & contact_pressure = contact.getGlobalContactPressure();
  auto & friction_traction = friction.getGlobalFrictionTraction();

  // The two functions contact.applyContactPressure and
  // friction.applyFrictionTraction() are replaced by assembling the forces to
  // the residual directly
  auto && dof_manager = solid.getDOFManager();
  dof_manager.assembleToResidual("displacement", contact_pressure, 1);
  dof_manager.assembleToResidual("displacement", friction_traction, 1);
}

} // namespace akantu
