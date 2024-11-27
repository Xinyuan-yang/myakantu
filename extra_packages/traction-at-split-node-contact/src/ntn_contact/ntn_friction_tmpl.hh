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
#include "ntn_contact.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
template <template <class> class FrictionLaw, class Regularisation>
NTNFriction<FrictionLaw, Regularisation>::NTNFriction(NTNBaseContact & contact,
                                                      const ID & id)
    : FrictionLaw<Regularisation>(contact, id) {
  AKANTU_DEBUG_IN();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <class> class FrictionLaw, class Regularisation>
void NTNFriction<FrictionLaw,
                 Regularisation>::assembleGlobalFrictionTraction() {
  AKANTU_DEBUG_IN();

  auto & ntn_contact = dynamic_cast<NTNContact &>(this->contact);
  auto & model = ntn_contact.getModel();
  auto dim = model.getSpatialDimension();

  auto & global_friction_traction = this->getGlobalFrictionTraction();
  global_friction_traction.zero();

  const auto & masters = ntn_contact.getMasters();
  const auto & slaves = ntn_contact.getSlaves();
  const auto & l_boundary_slaves = ntn_contact.getLumpedBoundarySlaves();
  const auto & l_boundary_masters = ntn_contact.getLumpedBoundaryMasters();

  Int nb_contact_nodes = ntn_contact.getNbContactNodes();
  for (Int n = 0; n < nb_contact_nodes; ++n) {
    auto master = masters(n);
    auto slave = slaves(n);

    for (Int d = 0; d < dim; ++d) {
      global_friction_traction(master, d) +=
          l_boundary_masters(n) * this->friction_traction(n, d);
      global_friction_traction(slave, d) -=
          l_boundary_slaves(n) * this->friction_traction(n, d);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <template <class> class FrictionLaw, class Regularisation>
void NTNFriction<FrictionLaw, Regularisation>::printself(std::ostream & stream,
                                                         int indent) const {
  AKANTU_DEBUG_IN();
  std::string space(indent, AKANTU_INDENT);

  stream << space << "NTNFriction [" << std::endl;
  FrictionLaw<Regularisation>::printself(stream, ++indent);
  stream << space << "]" << std::endl;

  AKANTU_DEBUG_OUT();
}

} // namespace akantu
