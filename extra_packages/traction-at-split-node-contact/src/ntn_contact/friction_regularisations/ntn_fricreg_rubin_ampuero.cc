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
#include "ntn_fricreg_rubin_ampuero.hh"
#include "dumper_nodal_field.hh"
#include "dumper_text.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
NTNFricRegRubinAmpuero::NTNFricRegRubinAmpuero(NTNBaseContact & contact,
                                               const ID & id)
    : NTNFricRegNoRegularisation(contact, id),
      t_star(0, 1, 0., id + ":t_star", 0., "t_star") {
  AKANTU_DEBUG_IN();

  NTNFricRegNoRegularisation::registerSynchronizedArray(this->t_star);

  this->registerParam("t_star", this->t_star, _pat_parsmod,
                      "time scale of regularization");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
const SynchronizedArray<Real> &
NTNFricRegRubinAmpuero::internalGetContactPressure() {
  AKANTU_DEBUG_IN();

  auto & model = this->contact.getModel();
  auto dim = model.getSpatialDimension();
  auto delta_t = model.getTimeStep();

  // get contact arrays
  const auto & is_in_contact = this->internalGetIsInContact();
  const auto & pressure = this->contact.getContactPressure().getArray();
  auto it = pressure.begin(dim);

  auto nb_contact_nodes = this->contact.getNbContactNodes();
  for (Int n = 0; n < nb_contact_nodes; ++n) {
    // node pair is NOT in contact
    if (not is_in_contact(n)) {
      this->frictional_contact_pressure(n) = 0.;

      // if t_star is too small compute like Coulomb friction (without
      // regularization)
    } else if (Math::are_float_equal(this->t_star(n), 0.)) {
      const auto & pres = it[n];
      this->frictional_contact_pressure(n) = pres.norm();
    }

    else {
      // compute frictional contact pressure
      // backward euler method: first order implicit numerical integration
      // method
      // \reg_pres_n+1 = (\reg_pres_n + \delta_t / \t_star * \cur_pres)
      //               / (1 + \delta_t / \t_star)
      auto alpha = delta_t / this->t_star(n);
      const auto & pres = it[n];
      this->frictional_contact_pressure(n) += alpha * pres.norm();
      this->frictional_contact_pressure(n) /= 1 + alpha;
    }
  }

  AKANTU_DEBUG_OUT();
  return this->frictional_contact_pressure;
}

/* -------------------------------------------------------------------------- */
void NTNFricRegRubinAmpuero::setToSteadyState() {
  NTNFricRegNoRegularisation::computeFrictionalContactPressure();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegRubinAmpuero::registerSynchronizedArray(
    SynchronizedArrayBase & array) {
  this->t_star.registerDependingArray(array);
}

/* -------------------------------------------------------------------------- */
void NTNFricRegRubinAmpuero::dumpRestart(const std::string & file_name) const {
  AKANTU_DEBUG_IN();

  this->t_star.dumpRestartFile(file_name);

  NTNFricRegNoRegularisation::dumpRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegRubinAmpuero::readRestart(const std::string & file_name) {
  AKANTU_DEBUG_IN();

  this->t_star.readRestartFile(file_name);

  NTNFricRegNoRegularisation::readRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegRubinAmpuero::printself(std::ostream & stream,
                                       int indent) const {
  AKANTU_DEBUG_IN();
  std::string space(indent, AKANTU_INDENT);

  stream << space << "NTNFricRegRubinAmpuero [" << std::endl;
  NTNFricRegNoRegularisation::printself(stream, ++indent);
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NTNFricRegRubinAmpuero::addDumpFieldToDumper(
    const std::string & dumper_name, const std::string & field_id) {
  AKANTU_DEBUG_IN();

  if (field_id == "t_star") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(this->t_star.getArray()));
  } else {
    NTNFricRegNoRegularisation::addDumpFieldToDumper(dumper_name, field_id);
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
