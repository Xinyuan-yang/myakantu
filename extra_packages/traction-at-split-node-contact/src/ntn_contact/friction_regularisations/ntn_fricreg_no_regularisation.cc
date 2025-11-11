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
#include "ntn_fricreg_no_regularisation.hh"
#include "dumper_nodal_field.hh"
#include "dumper_text.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
NTNFricRegNoRegularisation::NTNFricRegNoRegularisation(NTNBaseContact & contact,
                                                       const ID & id)
    : NTNBaseFriction(contact, id),
      frictional_contact_pressure(0, 1, 0., id + ":frictional_contact_pressure",
                                  0., "frictional_contact_pressure") {
  NTNBaseFriction::registerSynchronizedArray(this->frictional_contact_pressure);

  this->registerParam("frictional_contact_pressure",
                      this->frictional_contact_pressure, _pat_internal,
                      "contact pressure used for friction law");
}

/* -------------------------------------------------------------------------- */
const SynchronizedArray<Real> &
NTNFricRegNoRegularisation::internalGetContactPressure() {
  this->computeFrictionalContactPressure();
  return this->frictional_contact_pressure;
}

/* -------------------------------------------------------------------------- */
void NTNFricRegNoRegularisation::computeFrictionalContactPressure() {
  auto & model = this->contact.getModel();
  auto dim = model.getSpatialDimension();

  // get contact arrays
  const auto & is_in_contact = this->internalGetIsInContact();
  const auto & pressure = this->contact.getContactPressure().getArray();
  auto it = pressure.begin(dim);

  auto nb_contact_nodes = this->contact.getNbContactNodes();
  for (Int n = 0; n < nb_contact_nodes; ++n) {
    // node pair is NOT in contact
    if (not is_in_contact(n)) {
      this->frictional_contact_pressure(n) = 0.;

      // node pair is in contact
    } else {
      // compute frictional contact pressure
      const auto & pres = it[n];
      this->frictional_contact_pressure(n) = pres.norm();
    }
  }
}

/* -------------------------------------------------------------------------- */
void NTNFricRegNoRegularisation::registerSynchronizedArray(
    SynchronizedArrayBase & array) {
  this->frictional_contact_pressure.registerDependingArray(array);
}

/* -------------------------------------------------------------------------- */
void NTNFricRegNoRegularisation::dumpRestart(
    const std::string & file_name) const {
  this->frictional_contact_pressure.dumpRestartFile(file_name);

  NTNBaseFriction::dumpRestart(file_name);
}

/* -------------------------------------------------------------------------- */
void NTNFricRegNoRegularisation::readRestart(const std::string & file_name) {
  this->frictional_contact_pressure.readRestartFile(file_name);

  NTNBaseFriction::readRestart(file_name);
}

/* -------------------------------------------------------------------------- */
void NTNFricRegNoRegularisation::printself(std::ostream & stream,
                                           int indent) const {
  std::string space(indent, AKANTU_INDENT);
  stream << space << "NTNFricRegNoRegularisation [" << std::endl;
  NTNBaseFriction::printself(stream, ++indent);
  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
void NTNFricRegNoRegularisation::addDumpFieldToDumper(
    const std::string & dumper_name, const std::string & field_id) {
  AKANTU_DEBUG_IN();

  if (field_id == "frictional_contact_pressure") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(
            this->frictional_contact_pressure.getArray()));
  } else {
    NTNBaseFriction::addDumpFieldToDumper(dumper_name, field_id);
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
