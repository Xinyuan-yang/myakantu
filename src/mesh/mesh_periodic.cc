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
#include "aka_common.hh"
#include "aka_iterators.hh"
#include "aka_safe_enum.hh"
#include "communication_tag.hh"
#include "communicator.hh"
#include "element_group.hh"
#include "mesh.hh"
#include "periodic_node_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <array>
#include <cstdlib>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
void Mesh::makePeriodic(const SpatialDirection & direction) {
  Array<Idx> list_1;
  Array<Idx> list_2;
  Real tolerance = 1e-10;

  auto lower_bound = this->getLowerBounds();
  auto upper_bound = this->getUpperBounds();
  auto length = upper_bound(direction) - lower_bound(direction);

  const auto & positions = *nodes;

  for (auto && [node, pos] :
       enumerate(make_view(positions, spatial_dimension))) {
    if (std::abs((pos(direction) - lower_bound(direction)) / length) <
        tolerance) {
      list_1.push_back(node);
    }

    if (std::abs((pos(direction) - upper_bound(direction)) / length) <
        tolerance) {
      list_2.push_back(node);
    }
  }

  this->makePeriodic(direction, list_1, list_2);
}

/* -------------------------------------------------------------------------- */
void Mesh::makePeriodic(const SpatialDirection & direction, const ID & list_1,
                        const ID & list_2) {
  const auto & list_nodes_1 =
      mesh.getElementGroup(list_1).getNodeGroup().getNodes();
  const auto & list_nodes_2 =
      mesh.getElementGroup(list_2).getNodeGroup().getNodes();

  this->makePeriodic(direction, list_nodes_1, list_nodes_2);
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
namespace {
  struct NodeInfo {
    NodeInfo() = default;
    NodeInfo(Int spatial_dimension) : position(spatial_dimension) {}
    NodeInfo(Idx node, const Vector<Real> & position,
             const SpatialDirection & direction)
        : node(node), position(position) {
      this->direction_position = position(direction);
      this->position(direction) = 0.;
    }

    NodeInfo(const NodeInfo & other) = default;
    NodeInfo(NodeInfo && other) noexcept = default;
    auto operator=(const NodeInfo & other) -> NodeInfo & = default;
    auto operator=(NodeInfo && other) -> NodeInfo & = default;

    Idx node{0};
    Vector<Real> position;
    Real direction_position{0.};
  };

  enum Side { _left = 0, _right = 1, _side_end = 2 };
  struct SideDef {
    using type = Side;
    static const type _begin_ = _left;   // NOLINT
    static const type _end_ = _side_end; // NOLINT
  };

  using side_t = safe_enum<SideDef>;
  inline constexpr side_t sides{_side_end};

  inline auto operator!(const Side & a) -> Side {
    using under = std::underlying_type_t<Side>;
    return Side((under(a) + 1) % 2);
  }
} // namespace

/* -------------------------------------------------------------------------- */
// left is for lower values on direction and right for highest values
void Mesh::makePeriodic(const SpatialDirection & direction,
                        const Array<Idx> & list_left,
                        const Array<Idx> & list_right) {
  Real tolerance = 1e-10;

  const auto & positions = *nodes;
  auto lower_bound = this->getLowerBounds();
  auto upper_bound = this->getUpperBounds();
  auto length = upper_bound(direction) - lower_bound(direction);

  lower_bound(direction) = 0;
  upper_bound(direction) = 0;

  auto prank = communicator->whoAmI();

  std::array<std::vector<NodeInfo>, 2> node_lists;
  node_lists[_left].resize(list_left.size());
  node_lists[_right].resize(list_right.size());

  auto to_position = [&](auto && bbox, auto node) {
    Vector<Real> pos(spatial_dimension);
    for (auto s : arange(spatial_dimension)) {
      pos(s) = positions(node, s);
    }
    auto && info = NodeInfo(node, pos, direction);
    bbox += info.position;
    return std::move(info);
  };

  std::array<BBox, 2> local_bboxes{spatial_dimension, spatial_dimension};
  std::transform(
      list_left.begin(), list_left.end(), node_lists[_left].begin(),
      [&](auto node) { return to_position(local_bboxes[_left], node); });

  std::transform(
      list_right.begin(), list_right.end(), node_lists[_right].begin(),
      [&](auto node) { return to_position(local_bboxes[_right], node); });

  std::vector<Idx> new_nodes;
  if (is_distributed) {
    NewNodesEvent event(AKANTU_CURRENT_FUNCTION);

    /* ---------------------------------------------------------------------- */
    // function to send nodes in bboxes intersections
    auto extract_and_send_nodes = [&](const auto & bbox, const auto & node_list,
                                      auto & buffers, auto proc, auto side) {
      buffers.push_back(std::make_unique<DynamicCommunicationBuffer>());
      auto & buffer = *buffers.back();

      auto tag =
          Tag::genTag(prank, 10 * direction + side, Tag::_periodic_nodes);

      for (auto & info : node_list) {
        if (not bbox.contains(info.position) or
            not isLocalOrMasterNode(info.node)) {
          continue;
        }

        auto pos = info.position;
        pos(direction) = info.direction_position;

        auto flag = (*nodes_flags)(info.node) & NodeFlag::_periodic_mask;
        auto gnode = getNodeGlobalId(info.node);
        buffer << gnode;
        buffer << pos;
        buffer << flag;

        // if is slave sends master info
        if (flag == NodeFlag::_periodic_slave) {
          auto master = getNodeGlobalId(periodic_slave_master[info.node]);
          buffer << master;
        }

        // if is master sends list of slaves
        if (flag == NodeFlag::_periodic_master) {
          Int nb_slaves = periodic_master_slave.count(info.node);
          buffer << nb_slaves;

          auto slaves = periodic_master_slave.equal_range(info.node);
          for (auto && [_, slave] : range(slaves.first, slaves.second)) {
            auto gslave = getNodeGlobalId(slave);
            buffer << gslave;
          }
        }
      }

      AKANTU_DEBUG_TRACE("Sending PBC node date to " << proc << " - " << tag);
      return communicator->asyncSend(buffer, proc, tag);
    };

    /* ---------------------------------------------------------------------- */
    // function to receive nodes in bboxes intersections
    auto recv_and_extract_nodes = [&](auto & node_list, const auto proc,
                                      auto side) {
      DynamicCommunicationBuffer buffer;
      auto tag = Tag::genTag(proc, 10 * direction + side, Tag::_periodic_nodes);

      AKANTU_DEBUG_TRACE("Trying to receive PBC node data from "
                         << proc << " - " << tag);

      communicator->receive(buffer, proc, tag);

      while (not buffer.empty()) {
        Vector<Real> pos(spatial_dimension);
        Idx global_node{};
        NodeFlag flag{};
        buffer >> global_node;
        buffer >> pos;
        buffer >> flag;

        auto local_node = getNodeLocalId(global_node);

        // get the master info of is slave
        if (flag == NodeFlag::_periodic_slave) {
          Idx master_node{};
          buffer >> master_node;
        }

        // get the list of slaves if is master
        if (flag == NodeFlag::_periodic_master) {
          Int nb_slaves{};
          buffer >> nb_slaves;
          for (auto ns [[gnu::unused]] : arange(nb_slaves)) {
            Idx gslave_node{};
            buffer >> gslave_node;
          }
        }

        if (local_node != -1) {
          continue;
        }

        local_node = nodes->size();

        nodes->push_back(pos);
        nodes_global_ids->push_back(global_node);
        nodes_flags->push_back(flag | NodeFlag::_pure_ghost);
        nodes_prank[local_node] = proc;
        new_nodes.push_back(local_node);

        node_list.emplace_back(local_node, pos, direction);

        event.getList().push_back(local_node);
      }
    };

    /* ---------------------------------------------------------------------- */
    std::array intersections{
        local_bboxes[_right].intersection(local_bboxes[_left], *communicator),
        local_bboxes[_left].intersection(local_bboxes[_right], *communicator)};

    std::vector<CommunicationRequest> send_requests;
    std::vector<std::unique_ptr<DynamicCommunicationBuffer>> send_buffers;

    // if local left intersects with remote right send left nodes and vice
    // versa
    for (auto side : sides) {
      for (auto && [proc, intersection_with_proc] : intersections[not side]) {
        // Send local nodes if intersects with remote
        if (intersection_with_proc) {
          send_requests.push_back(
              extract_and_send_nodes(intersection_with_proc, node_lists[side],
                                     send_buffers, proc, side));
        }
      }
    }

    // if local right intersects with remote left receive new left nodes and
    // vice versa
    for (auto side : sides) {
      for (auto && [proc, intersection_with_proc] : intersections[side]) {
        // receive remote nodes if intersects with local
        if (intersection_with_proc) {
          recv_and_extract_nodes(node_lists[side], proc, side);
        }
      }
    }

    Communicator::waitAll(send_requests);
    Communicator::freeCommunicationRequest(send_requests);

    this->sendEvent(event);
  } // end distributed work

  // sort nodes based on their distance to lower corner
  for (auto side : sides) {
    std::sort(node_lists[side].begin(), node_lists[side].end(),
              [&](auto && info1, auto && info2) -> bool {
                return info1.position < info2.position;
              });
  }

  // function to change the master of nodes
  auto updating_master = [&](auto & old_master, auto & new_master) {
    if (old_master == new_master) {
      return;
    }

    auto slaves = periodic_master_slave.equal_range(old_master);
    AKANTU_DEBUG_ASSERT(
        isPeriodicMaster(
            old_master), // slaves.first != periodic_master_slave.end(),
        "Cannot update master " << old_master << ", its not a master node!");
    decltype(periodic_master_slave) tmp_master_slave;
    for (auto it = slaves.first; it != slaves.second; ++it) {
      auto slave = it->second;
      tmp_master_slave.insert(std::make_pair(new_master, slave));
      periodic_slave_master[slave] = new_master;
    }

    periodic_master_slave.erase(old_master);
    (*nodes_flags)[old_master] &= ~NodeFlag::_periodic_master;
    addPeriodicSlave(old_master, new_master);

    for (auto && data : tmp_master_slave) {
      addPeriodicSlave(data.second, data.first);
    }
  };

  // handling 2 nodes that are periodic
  auto match_found = [&](auto & info1, auto & info2) {
    const auto & node1 = info1.node;
    const auto & node2 = info2.node;

    auto master = node1;
    bool node1_side_master = false;
    if (isPeriodicMaster(node1)) {
      node1_side_master = true;
    } else if (isPeriodicSlave(node1)) {
      node1_side_master = true;
      master = periodic_slave_master[node1];
    }

    auto node2_master = node2;
    if (isPeriodicSlave(node2)) {
      node2_master = periodic_slave_master[node2];
    }

    if (node1_side_master) {
      if (isPeriodicSlave(node2)) {
        updating_master(node2_master, master);
        return;
      }

      if (isPeriodicMaster(node2)) {
        updating_master(node2, master);
        return;
      }

      addPeriodicSlave(node2, master);
    } else {
      if (isPeriodicSlave(node2)) {
        addPeriodicSlave(node1, node2_master);
        return;
      }

      if (isPeriodicMaster(node2)) {
        addPeriodicSlave(node1, node2);
        return;
      }

      addPeriodicSlave(node2, node1);
    }
  };

  // matching the nodes from 2 lists
  auto match_pairs = [&](auto & nodes_1, auto & nodes_2) {
    // Guillaume to Nico: It seems that the list of nodes is not sorted
    // as it was: therefore the loop cannot be truncated anymore.
    // Otherwise many pairs are missing.
    // I replaced (temporarily?) for the N^2 loop so as not to miss
    // any pbc pair.
    //

    // auto it = nodes_2.begin();
    // for every nodes in 1st list
    for (auto && info1 : nodes_1) {
      auto & pos1 = info1.position;
      // auto it_cur = it;

      // try to find a match in 2nd list
      for (auto && info2 : nodes_2) {
        // auto & info2 = *it_cur;
        auto & pos2 = info2.position;

        auto dist = pos1.distance(pos2) / length;
        if (dist < tolerance) {
          // handles the found matches
          match_found(info1, info2);
          // it = it_cur;
          break;
        }
      }
    }
  };

  match_pairs(node_lists[_left], node_lists[_right]);
  // match_pairs(nodes_right, nodes_left);

  this->updatePeriodicSynchronizer();

  this->is_periodic = true;
}

/* -------------------------------------------------------------------------- */
void Mesh::wipePeriodicInfo() {
  this->is_periodic = false;

  this->periodic_slave_master.clear();
  this->periodic_master_slave.clear();

  for (auto && flags : *nodes_flags) {
    flags &= ~NodeFlag::_periodic_mask;
  }
}

/* -------------------------------------------------------------------------- */
void Mesh::updatePeriodicSynchronizer() {
  if (not this->periodic_node_synchronizer) {
    this->periodic_node_synchronizer =
        std::make_unique<PeriodicNodeSynchronizer>(
            *this, this->getID() + ":periodic_synchronizer", false);
  }

  this->periodic_node_synchronizer->update();
}

/* -------------------------------------------------------------------------- */
template <class T>
void Mesh::makeConsistentForPeriodicity(Array<T> & array,
                                        ReduceOperation operation) const {
  if (not isPeriodic()) {
    return;
  }

  auto zero_pure_ghosts = [&](auto & array) {
    auto periodic_remote =
        NodeFlag::_periodic | NodeFlag::_pure_ghost | NodeFlag::_slave;
    auto periodic_ghost = NodeFlag::_periodic | NodeFlag::_pure_ghost;
    auto periodic_slave = NodeFlag::_periodic | NodeFlag::_slave;

    for (auto && [flag, val] :
         zip(*nodes_flags, make_view(array, array.getNbComponent()))) {
      if (auto mask = flag & periodic_remote;
          mask == periodic_slave or mask == periodic_ghost) {
        val.fill(T());
      }
    }
  };

  switch (operation) {
  case ReduceOperation::_id:
    periodic_node_synchronizer->synchronizeArrayPBC(array);
    break;
  case ReduceOperation::_sum:
    zero_pure_ghosts(array);
    periodic_node_synchronizer->reduceSynchronizeArrayPBC<AddOperation>(array);
    break;
  case ReduceOperation::_or:
    zero_pure_ghosts(array);
    periodic_node_synchronizer->reduceSynchronizeArrayPBC<OrOperation>(array);
    break;
  }
}

template void
Mesh::makeConsistentForPeriodicity<Real>(Array<Real> & array,
                                         ReduceOperation operation) const;
template void
Mesh::makeConsistentForPeriodicity<bool>(Array<bool> & array,
                                         ReduceOperation operation) const;

} // namespace akantu
