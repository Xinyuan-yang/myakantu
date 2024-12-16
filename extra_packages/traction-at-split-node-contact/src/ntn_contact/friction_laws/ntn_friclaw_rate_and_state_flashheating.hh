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
#ifndef AST_NTN_FRICLAW_RATE_AND_STATE_FLASHHEATING_HH
#define AST_NTN_FRICLAW_RATE_AND_STATE_FLASHHEATING_HH

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "aka_iterators.hh"
#include "dumper_nodal_field.hh"
#include "ntn_base_contact.hh"
#include "ntn_friclaw_coulomb.hh"
#include "ntn_friclaw_rate_and_state.hh"
#include "ntn_fricreg_no_regularisation.hh"
#include "synchronized_array.hh"
/* -------------------------------------------------------------------------- */
#include <cmath>
#include <string>
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
namespace akantu {

struct FlashHeatingRaSFunctions {
  template <class Args> inline static void stick(Args && args) {
    using namespace tuple;

    args["mu"_n] = args["mu_0"_n];
  }

  template <class Args> inline static void slip(Args && args) {
    using namespace tuple;
    using std::log;

    if (args["slip_velocity_norm"_n] < args["v_star"_n]) {
      args["friction_state"_n] = args["d_c"_n] / args["v_star"_n];
      args["mu"_n] =
          args["mu_0"_n] +
          args["a"_n] * log(args["v_star"_n] / args["v_star"_n]) +
          args["b"_n] * log(args["friction_state"_n] / args["phi_star"_n]);
    } else {
      args["friction_state"_n] +=
          (1 - args["friction_state"_n] * args["slip_velocity_norm"_n] /
                   args["d_c"_n]) *
          args["time_step"_n];

      auto mu =
          args["mu_0"_n] +
          args["a"_n] * log(args["slip_velocity_norm"_n] / args["v_star"_n]) +
          args["b"_n] * log(args["friction_state"_n] / args["phi_star"_n]);

      args["mu"_n] =
          args["mu_w"_n] +
          (mu - args["mu_w"_n]) /
              (1 + args["d_c"_n] / (args["friction_state"_n] * args["v_w"_n]));
    }
  }
};

/* -------------------------------------------------------------------------- */
template <class Regularisation = NTNFricRegNoRegularisation>
class NTNFricLawRateAndStateFlashHeating
    : public NTNFricLawRateAndState<Regularisation> {

  using parent = NTNFricLawRateAndState<Regularisation>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NTNFricLawRateAndStateFlashHeating(NTNBaseContact & contact,
                                     const ID & id = "rate_and_state")
      : parent(contact, id), mu_w(0, 1, 0., id + ":mu_w", 0., "mu_w"),
        v_w(0, 1, 0., id + ":v_w", 0., "v_w") {
    parent::registerSynchronizedArray(this->mu_w);
    parent::registerSynchronizedArray(this->v_w);
    this->registerParam("mu_w", this->mu_w, _pat_parsmod,
                        "coefficient of friction of flash heating");
    this->registerParam("v_w", this->v_w, _pat_parsmod,
                        "velocity of flash heating");
  };

  /* ------------------------------------------------------------------------ */
  inline auto getArguments() -> decltype(auto) {
    using namespace tuple;
    return zip_append(parent::getArguments(), "mu_w"_n = mu_w, "v_w"_n = v_w);
  }

  void computeFrictionCoefficient() override {
    for (auto && args : getArguments()) {
      if (args["is_sticking"_n]) {
        FlashHeatingRaSFunctions::stick(std::forward<decltype(args)>(args));
      } else {
        FlashHeatingRaSFunctions::slip(std::forward<decltype(args)>(args));
      }
    }
  }

  /* ------------------------------------------------------------------------ */
  /* Dumpable                                                                 */
  /* ------------------------------------------------------------------------ */

  /// dump restart file
  void dumpRestart(const std::string & file_name) const override {
    this->mu_w.dumpRestartFile(file_name);
    this->v_w.dumpRestartFile(file_name);
    parent::dumpRestart(file_name);
  }

  /// read restart file
  void readRestart(const std::string & file_name) override {
    this->mu_w.readRestartFile(file_name);
    this->v_w.readRestartFile(file_name);
    parent::readRestart(file_name);
  }

  void addDumpFieldToDumper(const std::string & dumper_name,
                            const std::string & field_id) override {

    if (field_id == "mu_w") {
      this->internalAddDumpFieldToDumper(
          dumper_name, field_id,
          std::make_unique<dumpers::NodalField<Real>>(this->mu_w.getArray()));
    } else if (field_id == "v_w") {
      this->internalAddDumpFieldToDumper(
          dumper_name, field_id,
          std::make_unique<dumpers::NodalField<Real>>(this->v_w.getArray()));
    } else {
      NTNFricLawCoulomb<Regularisation>::addDumpFieldToDumper(dumper_name,
                                                              field_id);
    }
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
protected:
  // coefficient of friction of flash heating
  SynchronizedArray<Real> mu_w;

  // velocity of flash heating
  SynchronizedArray<Real> v_w;
};

namespace {
  inline const bool rate_and_state_flashheating_registered =
      registerFriction<NTNFricLawRateAndStateFlashHeating>(
          "rate_and_state:flashheating");

} // namespace

} // namespace akantu

#endif /* __AST_NTN_FRICLAW_RATE_AND_STATE_FLASHHEATING_HH__ */
