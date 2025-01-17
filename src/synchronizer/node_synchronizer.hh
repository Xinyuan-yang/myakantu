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
#include "mesh_events.hh"
#include "synchronizer_impl.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NODE_SYNCHRONIZER_HH_
#define AKANTU_NODE_SYNCHRONIZER_HH_

namespace akantu {

class NodeSynchronizer : public MeshEventHandler, public SynchronizerImpl<Idx> {
public:
  NodeSynchronizer(Mesh & mesh, const ID & id = "element_synchronizer",
                   bool register_to_event_manager = true,
                   EventHandlerPriority event_priority = _ehp_synchronizer);

  ~NodeSynchronizer() override;

  Int sanityCheckDataSize(const Array<Idx> & nodes,
                          const SynchronizationTag & tag,
                          bool from_comm_desc) const override;
  void packSanityCheckData(CommunicationBuffer & buffer,
                           const Array<Idx> & nodes,
                           const SynchronizationTag & /*tag*/) const override;
  void unpackSanityCheckData(CommunicationBuffer & buffer,
                             const Array<Idx> & nodes,
                             const SynchronizationTag & tag, Int proc,
                             Int rank) const override;

  /// function to implement to react on  akantu::NewNodesEvent
  void onNodesAdded(const Array<Idx> &, const NewNodesEvent &) override;

  /// function to implement to react on  akantu::RemovedNodesEvent
  void onNodesRemoved(const Array<Idx> &, const Array<Idx> &,
                      const RemovedNodesEvent &) override {}
  /// function to implement to react on  akantu::NewElementsEvent
  void onElementsAdded(const Array<Element> & /*unused*/,
                       const NewElementsEvent & /*unused*/) override {}
  /// function to implement to react on  akantu::RemovedElementsEvent
  void onElementsRemoved(const Array<Element> &,
                         const ElementTypeMapArray<Idx> &,
                         const RemovedElementsEvent &) override {}
  /// function to implement to react on  akantu::ChangedElementsEvent
  void onElementsChanged(const Array<Element> &, const Array<Element> &,
                         const ElementTypeMapArray<Idx> &,
                         const ChangedElementsEvent &) override {}

  /* ------------------------------------------------------------------------ */
  NodeSynchronizer & operator=(const NodeSynchronizer & other) {
    copySchemes(other);
    return *this;
  }

  friend class NodeInfoPerProc;

protected:
  void fillEntityToSend(Array<Idx> & entities_to_send) override;

public:
  AKANTU_GET_MACRO(Mesh, mesh, Mesh &);

  inline Int canScatterSize() override;
  inline Int gatheredSize() override;

  inline Idx localToGlobalEntity(const Idx & local) override;

protected:
  Int getRank(const Idx & node) const final;

protected:
  Mesh & mesh;
};

} // namespace akantu

#include "node_synchronizer_inline_impl.hh"

#endif /* AKANTU_NODE_SYNCHRONIZER_HH_ */
