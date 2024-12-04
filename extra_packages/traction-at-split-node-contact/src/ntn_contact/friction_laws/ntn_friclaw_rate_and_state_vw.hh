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
#ifndef AST_NTN_FRICLAW_RATE_AND_STATE_VW_HH
#define AST_NTN_FRICLAW_RATE_AND_STATE_VW_HH

/* -------------------------------------------------------------------------- */
#include "aka_tuple_tools.hh"
#include "ntn_friclaw_rate_and_state.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {
/* -------------------------------------------------------------------------- */
struct VWRaSFunctions {
  // Revised Dieterich rate and state law
  template <class Args> inline static void stick(Args && args) {
    using namespace tuple;
    using std::log;
    args["friction_state"_n] += args["time_step"_n];
    args["mu"_n] =
        args["mu_0"_n] +
        args["a"_n] * log(1 + args["slip_velocity_norm"_n] / args["v_star"_n]) +
        args["b"_n] * log(1 + args["friction_state"_n] / args["phi_star"_n]) *
            args["mu_0"_n];
  }

  /* ------------------------------------------------------------------------ */
  template <class Args> inline static void slip(Args && args) {
    using std::log;
    args["friction_state"_n] +=
        (1 - args["friction_state"_n] * args["slip_velocity_norm"_n] /
                 args["d_c"_n]) *
        args["time_step"_n];
    args["mu"_n] =
        args["mu_0"_n] +
        args["a"_n] * log(1 + args["slip_velocity_norm"_n] / args["v_star"_n]) +
        args["b"_n] * log(1 + args["friction_state"_n] / args["phi_star"_n]);
  }
};

template <class Regularisation = NTNFricRegNoRegularisation>
using NTNFricLawRateAndStateVW =
    NTNFricLawRateAndState<VWRaSFunctions, Regularisation>;

} // namespace akantu

#endif /* AST_NTN_FRICLAW_RATE_AND_STATE_VW_HH */
