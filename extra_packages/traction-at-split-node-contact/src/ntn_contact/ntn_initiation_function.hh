/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_common.hh"
#include "boundary_condition.hh"
#include "ntn_base_contact.hh"
#include "ntn_base_friction.hh"
#include "ntn_friction.hh"
#include "ntrf_friction.hh"

// friction regularisations
#include "ntn_fricreg_rubin_ampuero.hh"
#include "ntn_fricreg_simplified_prakash_clifton.hh"
/* -------------------------------------------------------------------------- */
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>

#ifndef AKANTU_NTN_INITIATION_FUNCTION_HH
#define AKANTU_NTN_INITIATION_FUNCTION_HH

namespace akantu {

namespace {
using NTNFactory = Factory<NTNBaseFriction, ID, ID, NTNBaseContact &>;

template <bool is_ntn, template <class> class FrictionLaw, class FrictionReg>
using Contact_ =
    std::conditional_t<is_ntn, NTNFriction<FrictionLaw, FrictionReg>,
                       NTRFFriction<FrictionLaw, FrictionReg>>;

template <template <class> class FrictionLaw>
inline auto registerFriction(const ID & friction_law) -> bool {
  NTNFactory::getInstance().registerAllocator(
      friction_law,
      [](const ID & friction_reg,
         NTNBaseContact & contact) -> std::unique_ptr<NTNBaseFriction> {
        return tuple_dispatch<std::tuple<std::true_type, std::false_type>>(
            [&](auto is_ntn_param) -> std::unique_ptr<NTNBaseFriction> {
              const bool is_ntn = is_ntn_param.value;
              if (friction_reg == "rubin_ampuero") {
                return std::make_unique<
                    Contact_<is_ntn, FrictionLaw, NTNFricRegRubinAmpuero>>(
                    contact);
              }

              if (friction_reg == "simplified_prakash_clifton") {
                return std::make_unique<Contact_<
                    is_ntn, FrictionLaw, NTNFricRegSimplifiedPrakashClifton>>(
                    contact);
              }

              return std::make_unique<
                  Contact_<is_ntn, FrictionLaw, NTNFricRegRubinAmpuero>>(
                  contact);
            },
            contact.isNTNContact());
      });
  return true;
}
} // namespace

auto initializeNTNFriction(NTNBaseContact & contact)
    -> std::unique_ptr<NTNBaseFriction>;

auto initializeNTNFriction(NTNBaseContact & contact,
                           const std::string & friction_law,
                           const std::string & friction_reg)
    -> std::unique_ptr<NTNBaseFriction>;

} // namespace akantu
/* -------------------------------------------------------------------------- */
#endif /* AKANTU_NTN_INITIATION_FUNCTION_HH */
