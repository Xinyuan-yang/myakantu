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
#include "mIIasym_contact.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
MIIASYMContact::MIIASYMContact(SolidMechanicsModel & model, const ID & id)
    : NTRFContact(model, id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MIIASYMContact::updateImpedance() {
  AKANTU_DEBUG_IN();

  NTRFContact::updateImpedance();

  for (Int i = 0; i < this->impedance.size(); ++i) {
    this->impedance(i) *= 0.5;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
/// WARNING: this is only valid for the acceleration in equilibrium
void MIIASYMContact::computeRelativeNormalField(
    const Array<Real> & field, Array<Real> & rel_normal_field) const {
  AKANTU_DEBUG_IN();

  NTRFContact::computeRelativeNormalField(field, rel_normal_field);

  for (auto it_rtfield = rel_normal_field.begin();
       it_rtfield != rel_normal_field.end(); ++it_rtfield) {

    // in the anti-symmetric case
    *it_rtfield *= 2.;
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void MIIASYMContact::computeRelativeTangentialField(
    const Array<Real> & field, Array<Real> & rel_tang_field) const {
  NTRFContact::computeRelativeTangentialField(field, rel_tang_field);

  auto dim = this->model.getSpatialDimension();

  for (auto && rt : make_view(rel_tang_field, dim)) {
    // in the anti-symmetric case, the tangential fields become twice as large
    rt *= 2.;
  }
}

/* -------------------------------------------------------------------------- */
void MIIASYMContact::computeContactPressureInEquilibrium() {
  NTRFContact::computeContactPressure();
}

/* -------------------------------------------------------------------------- */
void MIIASYMContact::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "MIIASYMContact [" << std::endl;
  NTRFContact::printself(stream, indent);
  stream << space << "]" << std::endl;
}

} // namespace akantu
