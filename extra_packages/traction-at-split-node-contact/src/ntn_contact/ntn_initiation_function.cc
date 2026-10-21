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
#include "ntn_initiation_function.hh"
#include "aka_common.hh"
#include "ntn_base_contact.hh"
#include "ntn_base_friction.hh"
#include "ntn_friclaw_coulomb.hh"

// friction laws
#include "ntn_friclaw_linear_cohesive.hh"
#include "ntn_friclaw_linear_slip_weakening.hh"
#include "ntn_friclaw_linear_slip_weakening_no_healing.hh"

#include "ntn_friclaw_rate_and_state.hh"              // IWYU: keep
#include "ntn_friclaw_rate_and_state_flashheating.hh" // IWYU: keep
#include "ntn_friclaw_stribeck.hh" // Stribeck Friction 

#include "aka_factory.hh"

#include <memory>
#include <ostream>
#include <string>
#include <tuple>

namespace akantu {

/* -------------------------------------------------------------------------- */
auto initializeNTNFriction(NTNBaseContact & contact)
    -> std::unique_ptr<NTNBaseFriction> {
  AKANTU_DEBUG_IN();
  auto sub_sect = getStaticParser().getSubSections(ParserType::_friction);
  auto it = sub_sect.first;
  const ParserSection & section = *it;

  std::string friction_law = section.getName();
  std::string friction_reg = section.getOption("no_regularisation");

  std::unique_ptr<NTNBaseFriction> friction =
      NTNFactory::getInstance().allocate(friction_law, friction_reg, contact);

  friction->parseSection(section);

  if (++it != sub_sect.second) {
    AKANTU_DEBUG_WARNING("There were several friction sections in input file. "
                         << "Only first one was used and all others ignored.");
  }

  AKANTU_DEBUG_OUT();
  return friction;
}

namespace {

  auto registerFrictionLaws() -> bool {
    registerFriction<NTNFricLawCoulomb>("coulomb");
    registerFriction<NTNFricLawLinearSlipWeakening>("linear_slip_weakening");
    registerFriction<NTNFricLawLinearSlipWeakeningNoHealing>(
        "linear_slip_weakening_no_healing");
    registerFriction<NTNFricLawLinearCohesive>("linear_cohesive");
    registerFriction<NTNFricLawStribeck>("stribeck");

    return true;
  }

  const bool friction_laws_are_allocated [[maybe_unused]] =
      registerFrictionLaws();

} // namespace
} // namespace akantu

namespace std {
inline auto operator<<(ostream & stream,
                       const tuple<bool, string, string> & _this) -> ostream & {
  stream << "<" << get<0>(_this) << ", " << get<1>(_this) << ", "
         << get<2>(_this) << ">";
  return stream;
}
} // namespace std
