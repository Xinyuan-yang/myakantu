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
#include "group_manager.hh"
#include "aka_csr.hh"
#include "data_accessor.hh"
#include "element_group.hh"
#include "element_synchronizer.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "mesh_global_data_updater.hh"
#include "mesh_iterators.hh"
#include "mesh_utils.hh"
#include "node_group.hh"
#include "node_synchronizer.hh"
#include "periodic_node_synchronizer.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <iterator>
#include <list>
#include <numeric>
#include <queue>
#include <sstream>
#include <utility>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
GroupManager::GroupManager(Mesh & mesh, const ID & id) : id(id), mesh(mesh) {

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
GroupManager::~GroupManager() = default;

/* -------------------------------------------------------------------------- */
NodeGroup & GroupManager::createNodeGroup(const std::string & group_name,
                                          bool replace_group) {
  AKANTU_DEBUG_IN();

  auto it = node_groups.find(group_name);

  if (it != node_groups.end()) {
    if (replace_group) {
      it->second.reset();
    } else {
      AKANTU_EXCEPTION(
          "Trying to create a node group that already exists:" << group_name);
    }
  }

  std::stringstream sstr;
  sstr << this->id << ":" << group_name << "_node_group";

  auto && ptr = std::make_unique<NodeGroup>(group_name, mesh, sstr.str());

  auto & node_group = *ptr;

  // \todo insert_or_assign in c++17
  if (it != node_groups.end()) {
    it->second = std::move(ptr);
  } else {
    node_groups[group_name] = std::move(ptr);
  }

  AKANTU_DEBUG_OUT();

  return node_group;
}

/* -------------------------------------------------------------------------- */
template <typename T>
NodeGroup &
GroupManager::createFilteredNodeGroup(const std::string & group_name,
                                      const NodeGroup & source_node_group,
                                      T & filter) {
  AKANTU_DEBUG_IN();

  NodeGroup & node_group = this->createNodeGroup(group_name);
  node_group.append(source_node_group);
  if (T::type == FilterFunctor::_node_filter_functor) {
    node_group.applyNodeFilter(filter);
  } else {
    AKANTU_ERROR("ElementFilter cannot be applied to NodeGroup yet."
                 << " Needs to be implemented.");
  }

  AKANTU_DEBUG_OUT();
  return node_group;
}

/* -------------------------------------------------------------------------- */
ElementGroup & GroupManager::createElementGroup(const std::string & group_name,
                                                Int dimension,
                                                bool replace_group) {
  AKANTU_DEBUG_IN();

  auto it = element_groups.find(group_name);
  if (it != element_groups.end()) {
    if (replace_group) {
      it->second.reset();
    } else {
      AKANTU_EXCEPTION("Trying to create a element group that already exists:"
                       << group_name);
    }
  }

  auto & new_node_group = createNodeGroup(group_name + "_nodes", replace_group);

  auto && ptr = std::make_unique<ElementGroup>(
      group_name, mesh, new_node_group, dimension,
      this->id + ":" + group_name + "_element_group");

  auto & element_group = *ptr;

  if (it != element_groups.end()) {
    it->second = std::move(ptr);
  } else {
    element_groups[group_name] = std::move(ptr);
  }

  AKANTU_DEBUG_OUT();

  return element_group;
}

/* -------------------------------------------------------------------------- */
void GroupManager::destroyElementGroup(const std::string & group_name,
                                       bool destroy_node_group) {
  AKANTU_DEBUG_IN();

  auto eit = element_groups.find(group_name);
  if (eit != element_groups.end()) {
    if (destroy_node_group) {
      destroyNodeGroup(eit->second->getNodeGroup().getName());
    }
    element_groups.erase(eit);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GroupManager::destroyNodeGroup(const std::string & group_name) {
  AKANTU_DEBUG_IN();

  auto nit = node_groups.find(group_name);
  if (nit != node_groups.end()) {
    node_groups.erase(nit);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ElementGroup & GroupManager::createElementGroup(const std::string & group_name,
                                                Int dimension,
                                                NodeGroup & node_group) {
  AKANTU_DEBUG_IN();

  if (element_groups.find(group_name) != element_groups.end()) {
    AKANTU_EXCEPTION(
        "Trying to create a element group that already exists:" << group_name);
  }

  auto && ptr =
      std::make_unique<ElementGroup>(group_name, mesh, node_group, dimension,
                                     id + ":" + group_name + "_element_group");

  auto & element_group = *ptr;
  element_groups[group_name] = std::move(ptr);

  AKANTU_DEBUG_OUT();

  return element_group;
}

/* -------------------------------------------------------------------------- */
template <typename T>
ElementGroup & GroupManager::createFilteredElementGroup(
    const std::string & group_name, Int dimension, const NodeGroup & node_group,
    T & filter) {
  AKANTU_DEBUG_IN();

  if (T::type == FilterFunctor::_node_filter_functor) {
    auto & filtered_node_group = this->createFilteredNodeGroup(
        group_name + "_nodes", node_group, filter);

    AKANTU_DEBUG_OUT();
    return this->createElementGroup(group_name, dimension, filtered_node_group);
  }
  if (T::type == FilterFunctor::_element_filter_functor) {
    AKANTU_ERROR(
        "Cannot handle an ElementFilter yet. Needs to be implemented.");
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
class ClusterSynchronizer : public DataAccessor<Element> {
  using DistantIDs = std::set<std::pair<Idx, Idx>>;

public:
  ClusterSynchronizer(GroupManager & group_manager, Int element_dimension,
                      std::string cluster_name_prefix,
                      ElementTypeMapArray<Idx> & element_to_fragment,
                      const ElementSynchronizer & element_synchronizer,
                      Int nb_cluster)
      : group_manager(group_manager), element_dimension(element_dimension),
        cluster_name_prefix(std::move(cluster_name_prefix)),
        element_to_fragment(element_to_fragment),
        element_synchronizer(element_synchronizer), nb_cluster(nb_cluster) {}

  auto synchronize() {
    auto & comm = Communicator::getStaticCommunicator();
    auto rank = comm.whoAmI();
    auto nb_proc = comm.getNbProc();

    /// find starting index to renumber local clusters
    starting_index = nb_cluster;
    comm.exclusiveScan(starting_index);

    auto global_nb_fragment = starting_index + nb_cluster;
    comm.broadcast(global_nb_fragment, nb_proc - 1);

    /// create the local to distant cluster pairs with neighbors
    element_synchronizer.synchronizeOnce(*this,
                                         SynchronizationTag::_gm_clusters);

    /// count total number of pairs
    Array<Int> nb_pairs(nb_proc); // This is potentially a bug for more than
    // 2**31 pairs, but due to a all gatherv after
    // it must be int to match MPI interfaces
    nb_pairs(rank) = Int(distant_ids.size());
    comm.allGather(nb_pairs);

    auto local_pair_index =
        std::accumulate(nb_pairs.begin(), nb_pairs.begin() + rank, 0);

    auto total_nb_pairs = std::accumulate(nb_pairs.begin() + rank,
                                          nb_pairs.end(), local_pair_index);

    /// generate pairs global array
    Array<Int> total_pairs(total_nb_pairs, 2);

    for (const auto & ids : distant_ids) {
      total_pairs(local_pair_index, 0) = ids.first;
      total_pairs(local_pair_index, 1) = ids.second;
      ++local_pair_index;
    }

    /// communicate pairs to all processors
    nb_pairs *= 2;
    comm.allGatherV(total_pairs, nb_pairs);

    /// renumber clusters

    /// generate fragment list
    std::vector<std::set<Int>> global_clusters;
    Int total_nb_cluster = 0;

    Array<bool> is_fragment_in_cluster(global_nb_fragment, 1, false);
    std::queue<Int> fragment_check_list;

    while (not total_pairs.empty()) {
      /// create a new cluster
      ++total_nb_cluster;
      global_clusters.resize(total_nb_cluster);
      std::set<Int> & current_cluster = global_clusters[total_nb_cluster - 1];

      auto first_fragment = total_pairs(0, 0);
      auto second_fragment = total_pairs(0, 1);
      total_pairs.erase(0);

      fragment_check_list.push(first_fragment);
      fragment_check_list.push(second_fragment);

      while (!fragment_check_list.empty()) {
        auto current_fragment = fragment_check_list.front();
        auto * total_pairs_end = total_pairs.data() + total_pairs.size() * 2;
        auto * fragment_found =
            std::find(total_pairs.data(), total_pairs_end, current_fragment);

        if (fragment_found != total_pairs_end) {
          auto position = fragment_found - total_pairs.data();
          auto pair = position / 2;
          auto other_index = (position + 1) % 2;
          fragment_check_list.push(total_pairs(pair, other_index));
          total_pairs.erase(pair);
        } else {
          fragment_check_list.pop();
          current_cluster.insert(current_fragment);
          is_fragment_in_cluster(current_fragment) = true;
        }
      }
    }

    /// add to FragmentToCluster all local fragments
    for (auto c : arange(global_nb_fragment)) {
      if (!is_fragment_in_cluster(c)) {
        ++total_nb_cluster;
        global_clusters.resize(total_nb_cluster);
        auto & current_cluster = global_clusters[total_nb_cluster - 1];

        current_cluster.insert(c);
      }
    }

    /// reorganize element groups to match global clusters
    for (auto c : arange(global_clusters.size())) {

      /// create new element group corresponding to current cluster
      auto & cluster = group_manager.createElementGroup(
          cluster_name_prefix + "_" + std::to_string(c), element_dimension,
          true);

      /// append to current element group all fragments that belong to
      /// the same cluster if they exist
      for (auto gc : global_clusters[c]) {
        Int local_index = gc - starting_index;

        if (local_index < 0 || local_index >= Int(nb_cluster)) {
          continue;
        }

        auto id =
            "tmp_" + cluster_name_prefix + "_" + std::to_string(local_index);
        AKANTU_DEBUG_ASSERT(group_manager.elementGroupExists(id),
                            "Temporary fragment \"" << id << "\" not found");

        cluster.append(group_manager.getElementGroup(id));
        group_manager.destroyElementGroup(id, true);
      }
    }

    return total_nb_cluster;
  }

private:
  /// functions for parallel communications
  inline Int getNbData(const Array<Element> & elements,
                       const SynchronizationTag & tag) const override {
    if (tag == SynchronizationTag::_gm_clusters) {
      return elements.size() * sizeof(UInt);
    }

    return 0;
  }

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override {
    if (tag != SynchronizationTag::_gm_clusters) {
      return;
    }

    for (const auto & el : elements) {
      /// for each element pack its global cluster index
      buffer << element_to_fragment(el) + starting_index;
    }
  }

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override {
    if (tag != SynchronizationTag::_gm_clusters) {
      return;
    }

    for (const auto & el : elements) {
      Idx distant_cluster;
      buffer >> distant_cluster;

      auto local_cluster = element_to_fragment(el);
      if (local_cluster == -1) {
        continue;
      }

      distant_ids.insert(
          std::make_pair(local_cluster + starting_index, distant_cluster));
    }
  } // namespace akantu

private:
  GroupManager & group_manager;
  Int element_dimension{-1};
  std::string cluster_name_prefix;
  ElementTypeMapArray<Idx> & element_to_fragment;
  const ElementSynchronizer & element_synchronizer;

  Int nb_cluster{0};
  DistantIDs distant_ids;

  Idx starting_index{0};
}; // namespace akantu

/* -------------------------------------------------------------------------- */
/// \todo this function doesn't work in 1D
Int GroupManager::createBoundaryGroupFromGeometry() {
  auto spatial_dimension = mesh.getSpatialDimension();

  return createClusters(spatial_dimension - 1, "boundary");
}

/* -------------------------------------------------------------------------- */
Int GroupManager::createClusters(
    Int element_dimension, Mesh & mesh_facets, std::string cluster_name_prefix,
    const GroupManager::ClusteringFilter & filter) {
  return createClusters(element_dimension, cluster_name_prefix, filter,
                        mesh_facets);
}

/* -------------------------------------------------------------------------- */
Int GroupManager::createClusters(
    Int element_dimension, std::string cluster_name_prefix,
    const GroupManager::ClusteringFilter & filter) {

  MeshAccessor mesh_accessor(mesh);
  auto mesh_facets = std::make_unique<Mesh>(mesh.getSpatialDimension(),
                                            mesh_accessor.getNodesSharedPtr(),
                                            "mesh_facets_for_clusters");

  mesh_facets->defineMeshParent(mesh);
  MeshUtils::buildAllFacets(mesh, *mesh_facets, element_dimension,
                            element_dimension - 1);

  return createClusters(element_dimension, cluster_name_prefix, filter,
                        *mesh_facets);
}

/* -------------------------------------------------------------------------- */
//// \todo if needed element list construction can be optimized by
//// templating the filter class
Int GroupManager::createClusters(Int element_dimension,
                                 const std::string & cluster_name_prefix,
                                 const GroupManager::ClusteringFilter & filter,
                                 Mesh & mesh_facets) {
  AKANTU_DEBUG_IN();

  auto nb_proc = mesh.getCommunicator().getNbProc();
  std::string tmp_cluster_name_prefix = cluster_name_prefix;

  std::unique_ptr<ElementTypeMapArray<Idx>> element_to_fragment;

  if (nb_proc > 1 && mesh.isDistributed()) {
    element_to_fragment =
        std::make_unique<ElementTypeMapArray<Idx>>("element_to_fragment", id);

    element_to_fragment->initialize(
        mesh, _nb_component = 1, _spatial_dimension = element_dimension,
        _element_kind = _ek_not_defined, _with_nb_element = true,
        _default_value = -1);
    tmp_cluster_name_prefix = "tmp_" + tmp_cluster_name_prefix;
  }

  ElementTypeMapArray<bool> seen_elements("seen_elements", id);
  seen_elements.initialize(mesh, _spatial_dimension = element_dimension,
                           _element_kind = _ek_not_defined,
                           _with_nb_element = true, _default_value = false);

  for_each_element(
      mesh,
      [&filter, &seen_elements](auto && el) {
        seen_elements(el) = not filter(el);
      },
      _spatial_dimension = element_dimension);

  Int nb_cluster = 0;

  std::vector<std::string> created_groups;

  for (auto ghost_type : ghost_types) {
    Element uns_el;
    uns_el.ghost_type = ghost_type;
    for (auto type :
         mesh.elementTypes(_spatial_dimension = element_dimension,
         _ghost_type = ghost_type, _element_kind = _ek_not_defined)) {
      uns_el.type = type;
      auto & seen_elements_vec = seen_elements(uns_el.type, uns_el.ghost_type);

      for (auto && [e, seen_element] : enumerate(seen_elements_vec)) {
        // skip elements that have been already seen
        if (seen_element) {
          continue;
        }

        // set current element
        uns_el.element = e;
        seen_element = true;

        auto group_name =
            tmp_cluster_name_prefix + "_" + std::to_string(nb_cluster);

        created_groups.push_back(group_name);

        /// create a new cluster
        auto & cluster =
            createElementGroup(group_name, element_dimension, true);
        ++nb_cluster;

        // point element are cluster by themself
        if (element_dimension == 0) {
          cluster.add(uns_el, true, false);
          continue;
        }

        std::queue<Element> element_to_add;
        element_to_add.push(uns_el);

        /// keep looping until current cluster is complete (no more
        /// connected elements)
        while (not element_to_add.empty()) {

          /// take first element and erase it in the queue
          auto el = element_to_add.front();
          element_to_add.pop();

          /// if parallel, store cluster index per element
          if (nb_proc > 1 && mesh.isDistributed()) {
            (*element_to_fragment)(el) = nb_cluster - 1;
          }

          /// add current element to the cluster
          cluster.add(el, true, false);

          const auto & element_to_facets =
              mesh_facets.getSubelementToElement(el);

          for (auto && facet : element_to_facets) {
            if (facet == ElementNull) {
              continue;
            }

            const auto & connected_elements =
                mesh_facets.getElementToSubelement(facet);

            for (const auto & check_el : connected_elements) {

              // check if this element has to be skipped
              if (check_el == ElementNull || check_el == el) {
                continue;
              }

              auto & seen_elements_current = seen_elements(check_el);

              if (not seen_elements_current) {
                seen_elements_current = true;
                element_to_add.push(check_el);
              }
            }
          }
        }
      }
    }
  }

  for (auto && gid : created_groups) {
    this->getElementGroup(gid).getNodeGroup().optimize();
  }

  if (nb_proc > 1 && mesh.isDistributed()) {
    ClusterSynchronizer cluster_synchronizer(
        *this, element_dimension, cluster_name_prefix, *element_to_fragment,
        this->mesh.getElementSynchronizer(), nb_cluster);
    nb_cluster = cluster_synchronizer.synchronize();
  }

  if (mesh.isDistributed()) {
    this->synchronizeGroupNames();
  }

  AKANTU_DEBUG_OUT();
  return nb_cluster;
}

/* -------------------------------------------------------------------------- */
template <typename T>
void GroupManager::createGroupsFromMeshData(const std::string & dataset_name) {
  std::set<std::string> group_names;
  const auto & datas = mesh.getData<T>(dataset_name);

  std::map<std::string, Int> group_dim;

  for (auto ghost_type : ghost_types) {
    for (auto type : datas.elementTypes(_ghost_type = ghost_type)) {
      const auto & dataset = datas(type, ghost_type);
      auto nb_element = mesh.getNbElement(type, ghost_type);
      AKANTU_DEBUG_ASSERT(dataset.size() == nb_element,
                          "Not the same number of elements ("
                              << type << ":" << ghost_type
                              << ") in the map from MeshData ("
                              << dataset.size() << ") " << dataset_name
                              << " and in the mesh (" << nb_element << ")!");
      for (Int e(0); e < nb_element; ++e) {
        std::stringstream sstr;
        sstr << dataset(e);
        std::string gname = sstr.str();
        group_names.insert(gname);

        auto it = group_dim.find(gname);
        if (it == group_dim.end()) {
          group_dim[gname] = mesh.getSpatialDimension(type);
        } else {
          it->second = std::max(it->second, mesh.getSpatialDimension(type));
        }
      }
    }
  }

  for (auto && name : group_names) {
    createElementGroup(name, group_dim[name]);
  }

  if (mesh.isDistributed()) {
    this->synchronizeGroupNames();
  }

  Element el;
  for (auto ghost_type : ghost_types) {
    el.ghost_type = ghost_type;
    for (auto type : datas.elementTypes(_ghost_type = ghost_type)) {
      el.type = type;
      const auto & dataset = datas(type, ghost_type);
      auto nb_element = mesh.getNbElement(type, ghost_type);
      AKANTU_DEBUG_ASSERT(dataset.size() == nb_element,
                          "Not the same number of elements in the map from "
                          "MeshData and in the mesh!");

      auto nb_nodes_per_element = mesh.getNbNodesPerElement(el.type);

      auto cit =
          mesh.getConnectivity(type, ghost_type).begin(nb_nodes_per_element);

      for (Int e(0); e < nb_element; ++e, ++cit) {
        el.element = e;
        std::stringstream sstr;
        sstr << dataset(e);
        auto & group = getElementGroup(sstr.str());
        group.add(el, false, false);

        const auto & connect = *cit;
        for (auto node : connect) {
          group.addNode(node, false);
        }
      }
    }
  }

  for (auto && name : group_names) {
    getElementGroup(name).optimize();
  }
}

template void GroupManager::createGroupsFromMeshData<std::string>(
    const std::string & dataset_name);
template void
GroupManager::createGroupsFromMeshData<UInt>(const std::string & dataset_name);

/* -------------------------------------------------------------------------- */
void GroupManager::createElementGroupFromNodeGroup(
    const std::string & name, const std::string & node_group_name,
    Int dimension) {
  NodeGroup & node_group = getNodeGroup(node_group_name);
  ElementGroup & group = createElementGroup(name, dimension, node_group);

  group.fillFromNodeGroup();
}

/* -------------------------------------------------------------------------- */
void GroupManager::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "GroupManager [" << std::endl;

  std::set<std::string> node_group_seen;
  for (auto & group : iterateElementGroups()) {
    group.printself(stream, indent + 1);
    node_group_seen.insert(group.getNodeGroup().getName());
  }

  for (auto & group : iterateNodeGroups()) {
    if (node_group_seen.find(group.getName()) == node_group_seen.end()) {
      group.printself(stream, indent + 1);
    }
  }

  stream << space << "]" << std::endl;
}

/* -------------------------------------------------------------------------- */
Int GroupManager::getNbElementGroups(Int dimension) const {
  if (dimension == _all_dimensions) {
    return element_groups.size();
  }

  return std::count_if(element_groups.begin(), element_groups.end(),
                       [dimension](auto && eg) {
                         return eg.second->getDimension() == dimension;
                       });
}

/* -------------------------------------------------------------------------- */
void GroupManager::checkAndAddGroups(DynamicCommunicationBuffer & buffer) {
  AKANTU_DEBUG_IN();

  Int nb_node_group;
  buffer >> nb_node_group;
  AKANTU_DEBUG_INFO("Received " << nb_node_group << " node group names");

  for (Int ng = 0; ng < nb_node_group; ++ng) {
    std::string node_group_name;
    buffer >> node_group_name;

    if (node_groups.find(node_group_name) == node_groups.end()) {
      this->createNodeGroup(node_group_name);
    }

    AKANTU_DEBUG_INFO("Received node goup name: " << node_group_name);
  }

  Int nb_element_group;
  buffer >> nb_element_group;
  AKANTU_DEBUG_INFO("Received " << nb_element_group << " element group names");
  for (Int eg = 0; eg < nb_element_group; ++eg) {
    std::string element_group_name;
    buffer >> element_group_name;
    std::string node_group_name;
    buffer >> node_group_name;
    Int dim;
    buffer >> dim;

    AKANTU_DEBUG_INFO("Received element group name: "
                      << element_group_name << " corresponding to a "
                      << Int(dim) << "D group with node group "
                      << node_group_name);

    auto & node_group = *node_groups[node_group_name];

    if (element_groups.find(element_group_name) == element_groups.end()) {
      this->createElementGroup(element_group_name, dim, node_group);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GroupManager::fillBufferWithGroupNames(
    DynamicCommunicationBuffer & comm_buffer) const {
  AKANTU_DEBUG_IN();

  // packing node group names;
  Int nb_groups = this->node_groups.size();
  comm_buffer << nb_groups;
  AKANTU_DEBUG_INFO("Sending " << nb_groups << " node group names");

  auto nnames_it = node_groups.begin();
  auto nnames_end = node_groups.end();
  for (; nnames_it != nnames_end; ++nnames_it) {
    std::string node_group_name = nnames_it->first;
    comm_buffer << node_group_name;
    AKANTU_DEBUG_INFO("Sending node goupe name: " << node_group_name);
  }

  // packing element group names with there associated node group name
  nb_groups = this->element_groups.size();
  comm_buffer << nb_groups;
  AKANTU_DEBUG_INFO("Sending " << nb_groups << " element group names");

  for (auto && pair : this->element_groups) {
    auto & element_group = *(pair.second);
    std::string element_group_name = pair.first;
    std::string node_group_name = element_group.getNodeGroup().getName();
    Int dim = element_group.getDimension();

    comm_buffer << element_group_name;
    comm_buffer << node_group_name;
    comm_buffer << dim;

    AKANTU_DEBUG_INFO("Sending element group name: "
                      << element_group_name << " corresponding to a "
                      << Int(dim) << "D group with the node group "
                      << node_group_name);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void GroupManager::synchronizeGroupNames() {
  AKANTU_DEBUG_IN();

  const auto & comm = mesh.getCommunicator();
  auto nb_proc = comm.getNbProc();
  auto my_rank = comm.whoAmI();

  if (nb_proc == 1) {
    return;
  }

  if (my_rank == 0) {
    for (Int p = 1; p < nb_proc; ++p) {
      DynamicCommunicationBuffer recv_buffer;
      auto tag = Tag::genTag(p, 0, Tag::_element_group);
      comm.receive(recv_buffer, p, tag);
      AKANTU_DEBUG_INFO("Got " << printMemorySize<char>(recv_buffer.size())
                               << " from proc " << p << " " << tag);
      this->checkAndAddGroups(recv_buffer);
    }

    DynamicCommunicationBuffer comm_buffer;
    this->fillBufferWithGroupNames(comm_buffer);

    AKANTU_DEBUG_INFO("Initiating broadcast with "
                      << printMemorySize<char>(comm_buffer.size()));
    comm.broadcast(comm_buffer);

  } else {
    DynamicCommunicationBuffer comm_buffer;
    this->fillBufferWithGroupNames(comm_buffer);

    auto tag = Tag::genTag(my_rank, 0, Tag::_element_group);
    AKANTU_DEBUG_INFO("Sending " << printMemorySize<char>(comm_buffer.size())
                                 << " to proc " << 0 << " " << tag);
    comm.send(comm_buffer, 0, tag);

    DynamicCommunicationBuffer recv_buffer;
    comm.broadcast(recv_buffer);
    AKANTU_DEBUG_INFO("Receiving broadcast with "
                      << printMemorySize<char>(recv_buffer.size()));

    this->checkAndAddGroups(recv_buffer);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
const ElementGroup &
GroupManager::getElementGroup(const std::string & name) const {
  auto it = element_groups.find(name);
  if (it == element_groups.end()) {
    AKANTU_EXCEPTION("There are no element groups named "
                     << name << " associated to the group manager: " << id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
ElementGroup & GroupManager::getElementGroup(const std::string & name) {
  auto it = element_groups.find(name);
  if (it == element_groups.end()) {
    AKANTU_EXCEPTION("There are no element groups named "
                     << name << " associated to the group manager: " << id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
const NodeGroup & GroupManager::getNodeGroup(const std::string & name) const {
  auto it = node_groups.find(name);
  if (it == node_groups.end()) {
    AKANTU_EXCEPTION("There are no node groups named "
                     << name << " associated to the group manager: " << id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
NodeGroup & GroupManager::getNodeGroup(const std::string & name) {
  auto it = node_groups.find(name);
  if (it == node_groups.end()) {
    AKANTU_EXCEPTION("There are no node groups named "
                     << name << " associated to the group manager: " << id);
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
template <typename GroupsType>
void GroupManager::renameGroup(GroupsType & groups, const std::string & name,
                               const std::string & new_name) {
  auto it = groups.find(name);
  if (it == groups.end()) {
    AKANTU_EXCEPTION("There are no group named "
                     << name << " associated to the group manager: " << id);
  }

  auto && group_ptr = std::move(it->second);

  group_ptr->name = new_name;

  groups.erase(it);
  groups[new_name] = std::move(group_ptr);
}

/* -------------------------------------------------------------------------- */
void GroupManager::renameElementGroup(const std::string & name,
                                      const std::string & new_name) {
  renameGroup(element_groups, name, new_name);
}

/* -------------------------------------------------------------------------- */
void GroupManager::renameNodeGroup(const std::string & name,
                                   const std::string & new_name) {
  renameGroup(node_groups, name, new_name);
}

/* -------------------------------------------------------------------------- */
void GroupManager::copyElementGroup(const std::string & name,
                                    const std::string & new_name) {
  const auto & grp = getElementGroup(name);
  auto & new_grp = createElementGroup(new_name, grp.getDimension());

  new_grp.getElements().copy(grp.getElements());
}

/* -------------------------------------------------------------------------- */
void GroupManager::copyNodeGroup(const std::string & name,
                                 const std::string & new_name) {
  const auto & grp = getNodeGroup(name);
  auto & new_grp = createNodeGroup(new_name);

  new_grp.getNodes().copy(grp.getNodes());
}

/* -------------------------------------------------------------------------- */
void GroupManager::onNodesAdded(const Array<Idx> & new_nodes,
                                const NewNodesEvent & event) {
  for (auto && group : make_values_adaptor(element_groups)) {
    group->onNodesAdded(new_nodes, event);
  }
}
/* -------------------------------------------------------------------------- */

} // namespace akantu
