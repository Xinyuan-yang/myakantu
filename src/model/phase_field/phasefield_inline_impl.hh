/**
 * @file   phasefield_inline_impl.cc
 *
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation: Fri Jun 19 2020
 * @date last modification: Fri Apr 02 2021
 *
 * @brief  Phase field implementation of inline functions
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "phase_field_model.hh"
/* -------------------------------------------------------------------------- */

#ifndef __AKANTU_PHASEFIELD_INLINE_IMPL_HH__
#define __AKANTU_PHASEFIELD_INLINE_IMPL_HH__

namespace akantu { // namespace akantu
/* -------------------------------------------------------------------------- */
inline Int PhaseField::getNbData(const Array<Element> & elements,
                                 const SynchronizationTag & tag) const {
  if (tag == SynchronizationTag::_pfm_damage) {
    return (spatial_dimension + 2) * Int(sizeof(Real)) *
           this->getHandler().getNbIntegrationPoints(elements);
  }
  return 0;
}

/* -------------------------------------------------------------------------- */
inline void PhaseField::packData(CommunicationBuffer & buffer,
                                 const Array<Element> & elements,
                                 const SynchronizationTag & tag) const {

  if (tag == SynchronizationTag::_pfm_damage) {
    packInternalFieldHelper(driving_force, buffer, elements);
    packInternalFieldHelper(driving_energy, buffer, elements);
    packInternalFieldHelper(damage_on_qpoints, buffer, elements);
  }
}

/* -------------------------------------------------------------------------- */
inline void PhaseField::unpackData(CommunicationBuffer & buffer,
                                   const Array<Element> & elements,
                                   const SynchronizationTag & tag) {
  if (tag == SynchronizationTag::_pfm_damage) {
    unpackInternalFieldHelper(driving_force, buffer, elements);
    unpackInternalFieldHelper(driving_energy, buffer, elements);
    unpackInternalFieldHelper(damage_on_qpoints, buffer, elements);
  }
}
} // namespace akantu

#endif
