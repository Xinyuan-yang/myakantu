/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "node_synchronizer.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_PERIODIC_NODE_SYNCHRONIZER_HH_
#define AKANTU_PERIODIC_NODE_SYNCHRONIZER_HH_

namespace akantu {

class PeriodicNodeSynchronizer : public NodeSynchronizer {
public:
  PeriodicNodeSynchronizer(
      Mesh & mesh, const ID & id = "periodic_node_synchronizer",
      bool register_to_event_manager = true,
      EventHandlerPriority event_priority = _ehp_synchronizer);

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  void update();

  /// Uses the synchronizer to perform a reduction on the vector
  template <template <class> class Op, typename T>
  void reduceSynchronizeWithPBCSlaves(Array<T> & array) const;

  /// synchronize ghosts without state
  void synchronizeOnceImpl(DataAccessor<Idx> & data_accessor,
                           const SynchronizationTag & tag) const override;

  /// wait end of asynchronous synchronization of ghosts
  void waitEndSynchronizeImpl(DataAccessor<Idx> & data_accessor,
                              const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  // NodeSynchronizer master_to_slaves_synchronizer;
  Array<Idx> masters_list;
  Array<Idx> slaves_list;
};

/* -------------------------------------------------------------------------- */
template <template <class> class Op, typename T>
void PeriodicNodeSynchronizer::reduceSynchronizeWithPBCSlaves(
    Array<T> & array) const {
  ReduceDataAccessor<Idx, Op, T> data_accessor(array,
                                               SynchronizationTag::_whatever);
  auto size =
      data_accessor.getNbData(slaves_list, SynchronizationTag::_whatever);
  CommunicationBuffer buffer(size);

  data_accessor.packData(buffer, slaves_list, SynchronizationTag::_whatever);
  data_accessor.unpackData(buffer, masters_list, SynchronizationTag::_whatever);

  this->reduceSynchronizeArray<IdentityOperation>(array);
}

} // namespace akantu

#endif /* AKANTU_PERIODIC_NODE_SYNCHRONIZER_HH_ */
