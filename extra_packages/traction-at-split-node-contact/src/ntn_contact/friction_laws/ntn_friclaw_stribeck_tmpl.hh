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
#include "dumper_text.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Regularisation>
NTNFricLawStribeck<Regularisation>::NTNFricLawStribeck(
    NTNBaseContact & contact, const ID & id)
    : NTNFricLawCoulomb<Regularisation>(contact, id),
      mu_s(0, 1, 0., id + ":mu_s", 0., "mu_s"),
      mu_k(0, 1, 0., id + ":mu_k", 0., "mu_k"),
      v_0(0, 1, 0., id + ":v_0", 0., "v_0") {
  AKANTU_DEBUG_IN();
  //register Stribeck parameters 
  NTNFricLawCoulomb<Regularisation>::registerSynchronizedArray(this->mu_s);
  NTNFricLawCoulomb<Regularisation>::registerSynchronizedArray(this->mu_k);
  NTNFricLawCoulomb<Regularisation>::registerSynchronizedArray(this->v_0);

  this->registerParam("mu_s", this->mu_s, _pat_parsmod,
                      "static friction coefficient");
  this->registerParam("mu_k", this->mu_k, _pat_parsmod,
                      "kinetic friction coefficient");
  this->registerParam("v_0", this->v_0, _pat_parsmod, "critial velocity");

  this->setParameterAccessType("mu", _pat_readable);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawStribeck<
    Regularisation>::computeFrictionalStrength() {
  AKANTU_DEBUG_IN();

  computeFrictionCoefficient();
  NTNFricLawCoulomb<Regularisation>::computeFrictionalStrength();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawStribeck<
    Regularisation>::computeFrictionCoefficient() {
  AKANTU_DEBUG_IN();

  // get arrays
  const SynchronizedArray<bool> & stick = this->internalGetIsSticking();
  const SynchronizedArray<Real> & slip = this->internalGetSlip();
  const SynchronizedArray<Real> & slip_velocity = this->internalGetSlipVelocity();

  auto nb_contact_nodes = this->contact.getNbContactNodes();
  for (Int n = 0; n < nb_contact_nodes; ++n) {
    if (stick(n)) {
      this->mu(n) = this->mu_s(n);
    } 
    else {
      this->mu(n) = this->mu_k(n) + (this->mu_s(n) - this->mu_k(n)) *
                                  std::exp(-slip_velocity(n) / this->v_0(n));
      }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawStribeck<Regularisation>::registerSynchronizedArray(
    SynchronizedArrayBase & array) {
  AKANTU_DEBUG_IN();

  this->mu_s.registerDependingArray(array);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawStribeck<Regularisation>::dumpRestart(
    const std::string & file_name) const {
  AKANTU_DEBUG_IN();

  this->mu_s.dumpRestartFile(file_name);
  this->mu_k.dumpRestartFile(file_name);
  this->v_0.dumpRestartFile(file_name);

  NTNFricLawCoulomb<Regularisation>::dumpRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawStribeck<Regularisation>::readRestart(
    const std::string & file_name) {
  AKANTU_DEBUG_IN();

  this->mu_s.readRestartFile(file_name);
  this->mu_k.readRestartFile(file_name);
  this->v_0.readRestartFile(file_name);

  NTNFricLawCoulomb<Regularisation>::readRestart(file_name);

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawStribeck<Regularisation>::printself(
    std::ostream & stream, int indent) const {
  AKANTU_DEBUG_IN();
  std::string space;
  for (Int i = 0; i < indent; i++, space += AKANTU_INDENT) {
    ;
  }

  stream << space << "NTNFricLawStribeck [" << std::endl;
  NTNFricLawCoulomb<Regularisation>::printself(stream, ++indent);
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawStribeck<Regularisation>::addDumpFieldToDumper(
    const std::string & dumper_name, const std::string & field_id) {
  AKANTU_DEBUG_IN();

  if (field_id == "mu_s") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(this->mu_s.getArray()));
  } else if (field_id == "mu_k") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(this->mu_k.getArray()));
  } else if (field_id == "v_0") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(this->v_0.getArray()));
  } else {
    NTNFricLawCoulomb<Regularisation>::addDumpFieldToDumper(dumper_name,
                                                            field_id);
  }

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
