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

/* -------------------------------------------------------------------------- */
template <class Regularisation = NTNFricRegNoRegularisation>
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

/* -------------------------------------------------------------------------- */
template <class FrictionFunction, class StateFunction,
          class Regularisation = NTNFricRegNoRegularisation>
class NTNFricLawRateAndStateTemplated
    : public NTNFricLawRateAndState<Regularisation> {
  using parent = NTNFricLawRateAndState<Regularisation>;

public:
  using parent::parent;

  void computeFrictionCoefficient() override {
    for (auto && args : parent::getArguments()) {
      auto && mu = args["mu"_n];
      if (args["is_sticking"_n]) {
        mu = FrictionFunction::frictionCoefficient(
            std::forward<decltype(args)>(args));
      } else {
        args["friction_state"_n] +=
            StateFunction::stateIncrement(std::forward<decltype(args)>(args));
        mu = FrictionFunction::frictionCoefficient(
            std::forward<decltype(args)>(args));
      }
    }
  }
};

/* -------------------------------------------------------------------------- */
struct AgingStateFunction {
  template <class Args>
  inline static auto stateIncrement(Args && args) -> Real {
    return (1 - args["friction_state"_n] * args["slip_velocity_norm"_n] /
                    args["d_c"_n]) *
           args["time_step"_n];
  }
};

struct SlipStateFunction {
  template <class Args>
  inline static auto stateIncrement(Args && args) -> Real {
    using std::log;
    return (-args["friction_state"_n] * args["slip_velocity_norm"_n] /
            args["d_c"_n] *
            log(args["friction_state"_n] * args["slip_velocity_norm"_n] /
                args["d_c"_n])) *
           args["time_step"_n];
  }
};

struct AgingRegularizedStateFunction {
  template <class Args>
  inline static auto stateIncrement(Args && args) -> Real {
    using std::pow;
    using std::sqrt;
    return (1 - args["friction_state"_n] *
                    sqrt(pow(args["slip_velocity_norm"_n], 2) +
                         pow(args["v0"_n], 2)) /
                    args["d_c"_n]) *
           args["time_step"_n];
  }
};

} // namespace akantu

#include "ntn_friclaw_rate_and_state_tmpl.hh"

#include "ntn_initiation_function.hh"

namespace akantu {
#define REGISTER_RATE_AND_STATE(FricCoef, name)                                \
  namespace {                                                                  \
    template <class FrictionReg>                                               \
    using RASAging_##name =                                                    \
        NTNFricLawRateAndStateTemplated<FricCoef, AgingStateFunction,          \
                                        FrictionReg>;                          \
    inline const bool rate_and_state_aging_##name [[maybe_unused]] =           \
        registerFriction<RASAging_##name>("rate_and_state:" #name ":aging");   \
    template <class FrictionReg>                                               \
    using RASSlip_##name =                                                     \
        NTNFricLawRateAndStateTemplated<FricCoef, SlipStateFunction,           \
                                        FrictionReg>;                          \
    inline const bool rate_and_state_slip_##name [[maybe_unused]] =            \
        registerFriction<RASSlip_##name>("rate_and_state:" #name ":slip");     \
                                                                               \
    template <class FrictionReg>                                               \
    using RASAgingRegularized_##name = NTNFricLawRateAndStateTemplated<        \
        FricCoef, AgingRegularizedStateFunction, FrictionReg>;                 \
    inline const bool rate_and_state_aging_regularized_##name                  \
        [[maybe_unused]] = registerFriction<RASAgingRegularized_##name>(       \
            "rate_and_state:" #name ":aging_regularized");                     \
  } // namespace akantu

/* -------------------------------------------------------------------------- */
struct DieterichFrictionCoefficient {
  // Dieterich rate and state law
  // With regularized evolution law:
  //   A + 1 term is added inside the log of velocity to avoid NaN values for
  //   v = 0 sticking phase
  //
  /* ------------------------------------------------------------------------ */
  template <class Args>
  inline static auto frictionCoefficient(Args && args) -> Real {
    using std::log;
    using std::pow;
    return args["mu_0"_n] +
           args["a"_n] *
               log(1 + args["slip_velocity_norm"_n] / args["v_star"_n]) +
           args["b"_n] * log(args["friction_state"_n] / args["phi_star"_n]) *
               args["mu_0"_n];
  }
};

REGISTER_RATE_AND_STATE(DieterichFrictionCoefficient, dietrich)

/* -------------------------------------------------------------------------- */
struct VelocityWeakeningFrictionCoefficient {
  // Revised Dieterich rate and state law
  // Make sure that b and f_0 are consistent with the Dieterich law
  /* ------------------------------------------------------------------------ */
  template <class Args>
  inline static auto frictionCoefficient(Args && args) -> Real {
    using std::log;
    return args["mu_0"_n] +
           args["a"_n] *
               log(1 + args["slip_velocity_norm"_n] / args["v_star"_n]) +
           args["b"_n] *
               log(1 + args["friction_state"_n] / args["phi_star"_n]) *
               args["mu_0"_n];
  }
};

REGISTER_RATE_AND_STATE(VelocityWeakeningFrictionCoefficient,
                        velocity_weakening)

/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AST_NTN_FRICLAW_RATE_AND_STATE_HH */
