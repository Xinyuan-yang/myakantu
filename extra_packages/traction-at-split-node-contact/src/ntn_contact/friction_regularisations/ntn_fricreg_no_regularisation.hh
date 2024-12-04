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

#ifndef AST_NTN_FRICREG_NO_REGULARISATION_HH_
#define AST_NTN_FRICREG_NO_REGULARISATION_HH_

/* -------------------------------------------------------------------------- */
// simtools
#include "ntn_base_friction.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
class NTNFricRegNoRegularisation : public NTNBaseFriction {
  using parent = NTNBaseFriction;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NTNFricRegNoRegularisation(NTNBaseContact & contact,
                             const ID & id = "no_regularisation");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  /// set to steady state for no regularisation -> do nothing
  void setToSteadyState() override{};

  void registerSynchronizedArray(SynchronizedArrayBase & array) override;
  void dumpRestart(const std::string & file_name) const override;
  void readRestart(const std::string & file_name) override;

  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

protected:
  virtual void computeFrictionalContactPressure();

  /// compute frictional strength according to friction law
  void computeFrictionalStrength() override{};

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  void addDumpFieldToDumper(const std::string & dumper_name,
                            const std::string & field_id) override;

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
  auto getArguments() -> decltype(auto) {
    using namespace tuple;
    return zip_append(parent::getArguments(),
                      "frictional_contact_pressure"_n =
                          this->frictional_contact_pressure);
  }

protected:
  /// get the is_in_contact array
  virtual const SynchronizedArray<bool> & internalGetIsInContact() {
    return this->contact.getIsInContact();
  };

  /// get the contact pressure (the norm: scalar value)
  virtual const SynchronizedArray<Real> & internalGetContactPressure();

  /// get the frictional strength array
  virtual SynchronizedArray<Real> & internalGetFrictionalStrength() {
    return this->frictional_strength;
  };

  /// get the is_sticking array
  virtual SynchronizedArray<bool> & internalGetIsSticking() {
    return this->is_sticking;
  }

  /// get the slip array
  virtual SynchronizedArray<Real> & internalGetSlip() { return this->slip; }

  /// get the slip array
  virtual SynchronizedArray<Real> & internalGetCumulativeSlip() {
    return this->cumulative_slip;
  }

  /// get the slip array
  virtual SynchronizedArray<Real> & internalGetSlipVelocity() {
    return this->slip_velocity_norm;
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
  // contact pressure (absolut value) for computation of friction
  SynchronizedArray<Real> frictional_contact_pressure;
};

} // namespace akantu

#endif /* AST_NTN_FRICREG_NO_REGULARISATION_HH_ */
