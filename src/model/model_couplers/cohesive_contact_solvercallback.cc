/**
 * Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "cohesive_contact_solvercallback.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

CohesiveContactSolverCallback::CohesiveContactSolverCallback(
    SolidMechanicsModelCohesive & solid, ContactMechanicsModel & contact,
    AnalysisMethod & method)
    : solid(solid), contact(contact), method(method) {}

void CohesiveContactSolverCallback::assembleMatrix(const ID & matrix_id) {

  if (matrix_id == "K") {
    solid.assembleStiffnessMatrix();

    switch (method) {
    case _static:
    case _implicit_dynamic: {
      contact.assembleStiffnessMatrix();
      break;
    }
    default:
      break;
    }

  } else if (matrix_id == "M") {
    solid.assembleMass();
  }
}

/* -------------------------------------------------------------------------- */
void CohesiveContactSolverCallback::assembleLumpedMatrix(const ID & matrix_id) {
  if (matrix_id == "M") {
    solid.assembleMassLumped();
  }
}

/* -------------------------------------------------------------------------- */
void CohesiveContactSolverCallback::assembleResidual() {
  // computes the internal forces

  switch (method) {
  case _explicit_lumped_mass: {
    auto & current_positions = contact.getContactDetector().getPositions();
    current_positions.copy(solid.getCurrentPosition());
    contact.search();
    break;
  }
  default:
    break;
  }

  solid.assembleInternalForces();
  contact.assembleInternalForces();

  auto & internal_force = solid.getInternalForce();
  auto & external_force = solid.getExternalForce();

  auto & contact_force = contact.getInternalForce();

  solid.getDOFManager().assembleToResidual("displacement", external_force, 1);
  solid.getDOFManager().assembleToResidual("displacement", internal_force, 1);
  solid.getDOFManager().assembleToResidual("displacement", contact_force, 1);
}

/* -------------------------------------------------------------------------- */
void CohesiveContactSolverCallback::predictor() {
  auto & solid_model_solver = aka::as_type<ModelSolver>(solid);
  solid_model_solver.predictor();
}

/* -------------------------------------------------------------------------- */
void CohesiveContactSolverCallback::beforeSolveStep() {
  auto & solid_model_solver = aka::as_type<ModelSolver>(solid);
  solid_model_solver.beforeSolveStep();
}

/* -------------------------------------------------------------------------- */
void CohesiveContactSolverCallback::afterSolveStep(bool converged) {
  auto & solid_model_solver = aka::as_type<ModelSolver>(solid);
  solid_model_solver.afterSolveStep(converged);
}

/* -------------------------------------------------------------------------- */
void CohesiveContactSolverCallback::corrector() {
  auto & solid_model_solver = aka::as_type<ModelSolver>(solid);
  solid_model_solver.corrector();

  switch (method) {
  case _static:
  case _implicit_dynamic: {
    auto & current_positions = contact.getContactDetector().getPositions();
    current_positions.copy(solid.getCurrentPosition());
    contact.search();
    break;
  }
  default:
    break;
  }
}

/* -------------------------------------------------------------------------- */
MatrixType
CohesiveContactSolverCallback::getMatrixType(const ID & /*matrix_id*/) const {
  return _symmetric;
}

} // namespace akantu
