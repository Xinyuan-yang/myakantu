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
#include "aka_common.hh"
#include "ntn_base_contact.hh"
#include "ntn_base_friction.hh"
#include "ntn_contact.hh"
#include "solid_mechanics_model.hh"
#include "synchronized_array.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NTN_CONTACT_SOLVERCALLBACK_HH_
#define AKANTU_NTN_CONTACT_SOLVERCALLBACK_HH_

/* -------------------------------------------------------------------------- */

namespace akantu {

class NTNContactSolverCallback : public InterceptSolverCallback {

public:
  NTNContactSolverCallback(SolidMechanicsModel & solid, const ID & slave,
                           const ID & master, Int surface_normal_dir,
                           Real time_step_factor);

  /// implementation of SolverCallback::assembleResidual
  void assembleResidual() override;

  auto getContactPressure() -> const SynchronizedArray<Real> & {
    return this->contact->getContactPressure();
  }

  auto getNbContactNodes() -> Int { return this->contact->getNbContactNodes(); }

  [[nodiscard]] auto getFriction() const -> std::shared_ptr<NTNBaseFriction> {
    return friction;
  }

  [[nodiscard]] auto getContact() const -> std::shared_ptr<NTNContact> {
    return contact;
  }

  [[nodiscard]] auto getExternalWorkIncrement() const -> Vector<Real, 2>;
  [[nodiscard]] auto getExternalWorkIncrementNew() const -> Vector<Real, 2>;
  [[nodiscard]] auto getExternalWork() const -> Vector<Real, 2>;

  void beforeSolveStep() override;
  void afterSolveStep(bool converged = true) override;

private:
  /// model for the solid mechanics part of the coupling
  SolidMechanicsModel & solid;

  /// model for the contact part of the node to node algorithm
  std::shared_ptr<NTNContact> contact;

  /// model for the contact part of the node to node algorithm
  std::shared_ptr<NTNBaseFriction> friction;

  Vector<Real, 2> external_work{0., 0.};
  Int external_work_release{-1};
};

} // namespace akantu

#endif /* AKANTU_NTN_CONTACT_SOLVERCALLBACK_HH_ */
