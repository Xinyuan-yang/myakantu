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
#include "periodic_node_synchronizer.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
PeriodicNodeSynchronizer::PeriodicNodeSynchronizer(
    Mesh & mesh, const ID & id, const bool register_to_event_manager,
    EventHandlerPriority event_priority)
    : NodeSynchronizer(mesh, id + ":masters", register_to_event_manager,
                       event_priority) {}

/* -------------------------------------------------------------------------- */
void PeriodicNodeSynchronizer::update() {
  static int count = 0;
  const auto & masters_to_slaves = this->mesh.getPeriodicMasterSlaves();
  masters_list.resize(0);
  masters_list.reserve(masters_to_slaves.size());

  slaves_list.resize(0);
  slaves_list.reserve(masters_to_slaves.size());

  reset();

  std::set<Idx> masters_to_receive;
  for (auto && [master, slave] : masters_to_slaves) {
    masters_list.push_back(master);
    slaves_list.push_back(slave);

    if (not(mesh.isMasterNode(master) or mesh.isLocalNode(master))) {
      masters_to_receive.insert(master);
    }
  }

  if (not mesh.isDistributed() or nb_proc == 1) {
    return;
  }

  std::map<Int, Array<Idx>> buffers;
  for (auto node : masters_to_receive) {
    auto && proc = mesh.getNodePrank(node);
    auto && scheme = this->communications.createRecvScheme(proc);
    scheme.push_back(node);

    buffers[proc].push_back(mesh.getNodeGlobalId(node));
  }

  auto tag = Tag::genTag(0, count, Tag::_modify_scheme);
  std::vector<CommunicationRequest> requests;
  for (auto && [proc, buffer] : buffers) {
    requests.push_back(communicator.asyncSend(buffer, proc, tag,
                                              CommunicationMode::_synchronous));
    std::cout << "Recv from proc : " << proc << " -> "
              << this->communications.getScheme(proc, _recv).size()
              << std::endl;
  }

  communicator.receiveAnyNumber<Idx>(
      requests,
      [&](auto && proc, auto && msg) {
        auto && scheme = this->communications.createSendScheme(proc);
        for (auto node : msg) {
          scheme.push_back(mesh.getNodeLocalId(node));
        }
        std::cout << "Send to proc : " << proc << " -> " << scheme.size()
                  << " [" << tag << "]" << std::endl;
      },
      tag);
  ++count;
}

/* -------------------------------------------------------------------------- */
void PeriodicNodeSynchronizer::synchronizeOnceImpl(
    DataAccessor<Idx> & data_accessor, const SynchronizationTag & tag) const {
  NodeSynchronizer::synchronizeOnceImpl(data_accessor, tag);

  auto size = data_accessor.getNbData(masters_list, tag);
  CommunicationBuffer buffer(size);

  data_accessor.packData(buffer, masters_list, tag);
  data_accessor.unpackData(buffer, slaves_list, tag);
}

/* -------------------------------------------------------------------------- */
void PeriodicNodeSynchronizer::waitEndSynchronizeImpl(
    DataAccessor<Idx> & data_accessor, const SynchronizationTag & tag) {
  NodeSynchronizer::waitEndSynchronizeImpl(data_accessor, tag);

  auto size = data_accessor.getNbData(masters_list, tag);
  CommunicationBuffer buffer(size);

  data_accessor.packData(buffer, masters_list, tag);
  data_accessor.unpackData(buffer, slaves_list, tag);
}

} // namespace akantu
