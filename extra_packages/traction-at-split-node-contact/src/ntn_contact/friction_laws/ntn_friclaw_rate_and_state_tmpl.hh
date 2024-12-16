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
#include "aka_common.hh"
#include "dumper_nodal_field.hh"
#include "dumper_text.hh"
#include "ntn_base_contact.hh"
#include "ntn_friclaw_coulomb.hh"
#include "ntn_friclaw_rate_and_state.hh"
#include "synchronized_array.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
#include <map>
#include <ostream>
#include <string>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class Regularisation>
NTNFricLawRateAndState<Regularisation>::NTNFricLawRateAndState(
    NTNBaseContact & contact, const ID & id)
    : parent(contact, id), mu_0(0, 1, 0., id + ":mu_0", 0., "mu_0"),
      a(0, 1, 0., id + ":a", 0., "a"), b(0, 1, 0., id + ":b", 0., "b"),
      v_star(0, 1, 0., id + ":v_star", 0., "v_star"),
      phi_star(0, 1, 0., id + ":phi_star", 0., "phi_star"),
      d_c(0, 1, 0., id + ":d_c", 0., "d_c"),
      v0(0, 1, 0., id + ":v0", 0., "v0") {
  parent::registerSynchronizedArray(this->mu_0);
  parent::registerSynchronizedArray(this->a);
  parent::registerSynchronizedArray(this->b);
  parent::registerSynchronizedArray(this->v_star);
  parent::registerSynchronizedArray(this->phi_star);
  parent::registerSynchronizedArray(this->d_c);
  parent::registerSynchronizedArray(this->v0);

  this->registerParam("mu_0", this->mu_0, _pat_parsmod,
                      "first coefficient of friction");
  this->registerParam("a", this->a, _pat_parsmod, "coefficient of second term");
  this->registerParam("b", this->b, _pat_parsmod, "coefficient of third term");
  this->registerParam("v_star", this->v_star, _pat_parsmod,
                      "velocity coefficient");
  this->registerParam("phi_star", this->phi_star, _pat_parsmod,
                      "state coefficient");
  this->registerParam("friction_state", this->friction_state, _pat_parsmod,
                      "state of friction");
  this->registerParam("d_c", this->d_c, _pat_parsmod, "characteristic length");
  this->registerParam("v0", this->v0, _pat_parsmod, "R&S parameter");

  this->setParameterAccessType("mu", _pat_readable);
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawRateAndState<Regularisation>::computeFrictionalStrength() {
  computeFrictionCoefficient();
  parent::computeFrictionalStrength();
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawRateAndState<Regularisation>::addNoiseToStateVariable(
    std::map<Real, Real> & noise) {
  std::cout << "Adding noise to the state field\n";

  const auto & cur_pos = this->contact->getModel().getCurrentPosition();
  const auto & slaves = this->internalGetSlaves();
  auto nb_contact_nodes = this->contact->getNbContactNodes();

  for (UInt n = 0; n < nb_contact_nodes; ++n) {
    UInt slave = slaves(n);

    auto it = noise.lower_bound(cur_pos(slave));

    Real perturbation = it->second;

    std::cout << "Slave :" << slave << " " << cur_pos(slave)
              << " : perturbation is " << perturbation << "\n";

    this->friction_state(n) = this->friction_state(n) * (1 + perturbation);
  }
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawRateAndState<Regularisation>::registerSynchronizedArray(
    SynchronizedArrayBase & array) {
  this->mu_0.registerDependingArray(array);
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawRateAndState<Regularisation>::dumpRestart(
    const std::string & file_name) const {

  this->mu_0.dumpRestartFile(file_name);
  this->a.dumpRestartFile(file_name);
  this->b.dumpRestartFile(file_name);
  this->phi_star.dumpRestartFile(file_name);
  this->v_star.dumpRestartFile(file_name);
  this->d_c.dumpRestartFile(file_name);
  this->v0.dumpRestartFile(file_name);
  this->friction_state.dumpRestartFile(file_name);

  parent::dumpRestart(file_name);
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawRateAndState<Regularisation>::readRestart(
    const std::string & file_name) {

  this->mu_0.readRestartFile(file_name);
  this->a.readRestartFile(file_name);
  this->b.readRestartFile(file_name);
  this->phi_star.readRestartFile(file_name);
  this->v_star.readRestartFile(file_name);
  this->d_c.readRestartFile(file_name);
  this->v0.readRestartFile(file_name);
  this->friction_state.readRestartFile(file_name);

  parent::readRestart(file_name);
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawRateAndState<Regularisation>::printself(std::ostream & stream,
                                                       int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "NTNFricLawRateAndState ["
         << "\n";
  NTNFricLawCoulomb<Regularisation>::printself(stream, ++indent);
  stream << space << "]"
         << "\n";
}

/* -------------------------------------------------------------------------- */
template <class Regularisation>
void NTNFricLawRateAndState<Regularisation>::addDumpFieldToDumper(
    const std::string & dumper_name, const std::string & field_id) {

  if (field_id == "mu_0") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(this->mu_0.getArray()));
  } else if (field_id == "a") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(this->a.getArray()));
  } else if (field_id == "b") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(this->b.getArray()));
  } else if (field_id == "phi_star") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(this->phi_star.getArray()));
  } else if (field_id == "v_star") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(this->v_star.getArray()));
  } else if (field_id == "friction_state") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(
            this->friction_state.getArray()));
  } else if (field_id == "d_c") {
    this->internalAddDumpFieldToDumper(
        dumper_name, field_id,
        std::make_unique<dumpers::NodalField<Real>>(this->d_c.getArray()));
  } else {
    NTNFricLawCoulomb<Regularisation>::addDumpFieldToDumper(dumper_name,
                                                            field_id);
  }
}

/* -------------------------------------------------------------------------- */

} // namespace akantu
