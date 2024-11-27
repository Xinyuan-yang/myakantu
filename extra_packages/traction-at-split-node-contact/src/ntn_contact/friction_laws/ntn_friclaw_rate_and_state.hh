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
#ifndef AST_NTN_FRICLAW_RATE_AND_STATE_HH
#define AST_NTN_FRICLAW_RATE_AND_STATE_HH

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "boundary_condition.hh"
#include "ntn_base_contact.hh"
#include "ntn_friclaw_coulomb.hh"
#include "ntn_fricreg_no_regularisation.hh"
#include "synchronized_array.hh"
/* -------------------------------------------------------------------------- */
#include <map>
#include <ostream>
#include <string>
/* -------------------------------------------------------------------------- */

namespace akantu {

struct NoOPRaSFunctions {
  template <class Args> inline static void stick(Args && args) {
    args["mu"_n] = 0;
  }

  template <class Args> inline static void slip(Args && args) {
    args["mu"_n] = 0;
  }
};

/* -------------------------------------------------------------------------- */
template <class RaSFunction, class Regularisation = NTNFricRegNoRegularisation>
class NTNFricLawRateAndState : public NTNFricLawCoulomb<Regularisation> {
  using parent = NTNFricLawCoulomb<Regularisation>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NTNFricLawRateAndState(NTNBaseContact & contact,
                         const ID & id = "rate_and_state_dieterich");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  /// register synchronizedarrays for sync
  void registerSynchronizedArray(SynchronizedArrayBase & array) override;

  /// dump restart file
  void dumpRestart(const std::string & file_name) const override;

  /// read restart file
  void readRestart(const std::string & file_name) override;

  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /// add noise to the state variable. A map is passed to each process
  /// containing the full noise profil to avoid any miscomunication errors
  void addNoiseToStateVariable(std::map<Real, Real> & noise);

protected:
  /// compute frictional strength according to friction law
  void computeFrictionalStrength() override;
  /// computes the friction coefficient as a function of slip
  virtual void computeFrictionCoefficient();

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
    return zip_append(parent::getArguments(), "mu_0"_n = mu_0, "a"_n = a,
                      "b"_n = b, "v_star"_n = v_star, "phi_star"_n = phi_star,
                      "d_c"_n = d_c, "v0"_n = v0,
                      "time_step"_n =
                          broadcast(this->contact.getModel().getTimeStep(),
                                    this->contact.getNbContactNodes()));
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  Real time{0.};
  Real step_count{0.};

  // first coefficient of friction
  SynchronizedArray<Real> mu_0;

  // coefficient of second term
  SynchronizedArray<Real> a;

  // coefficient of third term
  SynchronizedArray<Real> b;

  // velocity coefficient
  SynchronizedArray<Real> v_star;

  // state coefficient
  SynchronizedArray<Real> phi_star;

  // Dc the characteristic length
  SynchronizedArray<Real> d_c;

  // Revised R7S parameter
  SynchronizedArray<Real> v0;
};

} // namespace akantu

#include "ntn_friclaw_rate_and_state_tmpl.hh"

#endif /* AST_NTN_FRICLAW_RATE_AND_STATE_HH */
