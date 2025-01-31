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
#include "periodic_node_synchronizer.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PERIODIC_NODE_SYNCHRONIZER_TMPL_HH
#define AKANTU_PERIODIC_NODE_SYNCHRONIZER_TMPL_HH

namespace akantu {

/* -------------------------------------------------------------------------- */
template <typename T>
void PeriodicNodeSynchronizer::synchronizeArrayLocal(Array<T> & array) const {
  auto tag = SynchronizationTag::_whatever;

  SimpleIdxDataAccessor<T> data_accessor(array, SynchronizationTag::_whatever);

  auto size = data_accessor.getNbData(masters_list, tag);
  CommunicationBuffer buffer(size);

  data_accessor.packData(buffer, masters_list, SynchronizationTag::_whatever);
  data_accessor.unpackData(buffer, slaves_list, SynchronizationTag::_whatever);
}

/* -------------------------------------------------------------------------- */
template <typename T>
void PeriodicNodeSynchronizer::synchronizeArrayPBC(Array<T> & array) const {
  this->synchronizeArray(array);

  this->synchronizeArrayLocal(array);
}

/* -------------------------------------------------------------------------- */
template <template <class> class Op, typename T>
void PeriodicNodeSynchronizer::reduceSynchronizeArrayLocal(
    Array<T> & array) const {
  auto tag = SynchronizationTag::_whatever;

  ReduceDataAccessor<Idx, Op, T> data_accessor(array, tag);

  auto size = data_accessor.getNbData(slaves_list, tag);
  CommunicationBuffer buffer(size);

  data_accessor.packData(buffer, slaves_list, tag);
  data_accessor.unpackData(buffer, masters_list, tag);
}

/* -------------------------------------------------------------------------- */
template <template <class> class Op, typename T>
void PeriodicNodeSynchronizer::reduceSynchronizeArrayPBC(
    Array<T> & array) const {
  this->reduceSynchronizeArrayLocal<Op>(array);

  this->reduceSynchronizeArray<Op>(array);

  this->synchronizeArrayLocal(array);
}

} // namespace akantu

#endif /* AKANTU_PERIODIC_NODE_SYNCHRONIZER_TMPL_HH */
