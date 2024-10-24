/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_config.hh"
/* -------------------------------------------------------------------------- */
#include "element_class.hh"
#include "group_manager_inline_impl.hh"
#include "mesh.hh"
#include "mesh_global_data_updater.hh"
#include "mesh_io.hh"
#include "mesh_iterators.hh"
#include "mesh_utils.hh"
/* -------------------------------------------------------------------------- */
#include "communicator.hh"
#include "element_synchronizer.hh"
#include "facet_synchronizer.hh"
#include "mesh_utils_distribution.hh"
#include "node_synchronizer.hh"
#include "periodic_node_synchronizer.hh"
#if defined(AKANTU_COHESIVE_ELEMENT)
#include "cohesive_element_inserter.hh"
#endif
/* -------------------------------------------------------------------------- */
#include <algorithm>
/* -------------------------------------------------------------------------- */
#include "dumper_field.hh"
#include "dumper_internal_material_field.hh"
/* -------------------------------------------------------------------------- */
#include <limits>
#include <sstream>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
Mesh::Mesh(Int spatial_dimension, const ID & id, Communicator & communicator)
    : GroupManager(*this, id + ":group_manager"), MeshData("mesh_data", id),
      id(id), connectivities("connectivities", id),
      ghosts_counters("ghosts_counters", id),
      spatial_dimension(spatial_dimension),
      size(Vector<double>::Zero(spatial_dimension)), bbox(spatial_dimension),
      bbox_local(spatial_dimension), communicator(&communicator) {
  AKANTU_DEBUG_IN();
  size.fill(0.);
  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Mesh::Mesh(Int spatial_dimension, Communicator & communicator, const ID & id)
    : Mesh(spatial_dimension, id, communicator) {
  AKANTU_DEBUG_IN();

  this->nodes =
      std::make_shared<Array<Real>>(0, spatial_dimension, id + ":coordinates");
  this->nodes_flags = std::make_shared<Array<NodeFlag>>(0, 1, NodeFlag::_normal,
                                                        id + ":nodes_flags");

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
Mesh::Mesh(Int spatial_dimension, const ID & id)
    : Mesh(spatial_dimension, Communicator::getStaticCommunicator(), id) {}

/* -------------------------------------------------------------------------- */
Mesh::Mesh(Int spatial_dimension, const std::shared_ptr<Array<Real>> & nodes,
           const ID & id)
    : Mesh(spatial_dimension, id, Communicator::getStaticCommunicator()) {
  // NOLINTBEGIN(cppcoreguidelines-prefer-member-initializer)
  this->nodes = nodes;
  this->nb_global_nodes = this->nodes->size();
  // NOLINTEND(cppcoreguidelines-prefer-member-initializer)
  //
  this->nodes_to_elements.resize(nodes->size());
  for (auto & node_set : nodes_to_elements) {
    node_set = std::make_unique<std::set<Element>>();
  }

  this->computeBoundingBox();
}

/* -------------------------------------------------------------------------- */
Mesh::~Mesh() = default;

/* -------------------------------------------------------------------------- */
template <>
void Mesh::sendEvent<MeshIsDistributedEvent>(MeshIsDistributedEvent & event) {
  //    if(event.getList().size() != 0)
  EventHandlerManager<MeshEventHandler>::sendEvent(event);
}

/* -------------------------------------------------------------------------- */
template <> void Mesh::sendEvent<NewNodesEvent>(NewNodesEvent & event) {
  this->computeBoundingBox();
  this->nodes_flags->resize(this->nodes->size(), NodeFlag::_normal);

#if defined(AKANTU_COHESIVE_ELEMENT)
  if (aka::is_of_type<CohesiveNewNodesEvent>(event)) {
    // nodes might have changed in the connectivity
    const auto & mesh_to_mesh_facet =
        this->getData<Element>("mesh_to_mesh_facet");

    for (auto ghost_type : ghost_types) {
      for (auto type : connectivities.elementTypes(_spatial_dimension =
                                                       spatial_dimension - 1,
                       _ghost_type = ghost_type)) {
        auto nb_nodes_per_element = Mesh::getNbNodesPerElement(type);
        if (not mesh_to_mesh_facet.exists(type, ghost_type)) {
          continue;
        }

        const auto & mesh_to_mesh_facet_type =
            mesh_to_mesh_facet(type, ghost_type);
        auto && mesh_facet_conn_it =
            make_view(mesh_facets->connectivities(type, ghost_type),
                      nb_nodes_per_element)
                .begin();

        for (auto && [el, conn] : enumerate(make_view(
                 connectivities(type, ghost_type), nb_nodes_per_element))) {
          conn = mesh_facet_conn_it[mesh_to_mesh_facet_type(el).element];
        }
      }
    }
  }
#endif

  GroupManager::onNodesAdded(event.getList(), event);
  EventHandlerManager<MeshEventHandler>::sendEvent(event);
}

/* -------------------------------------------------------------------------- */
void Mesh::getBarycenters(Array<Real> & barycenter, ElementType type,
                          GhostType ghost_type) const {
  barycenter.resize(getNbElement(type, ghost_type));
  for (auto && data : enumerate(make_view(barycenter, spatial_dimension))) {
    getBarycenter(Element{type, Idx(std::get<0>(data)), ghost_type},
                  std::get<1>(data));
  }
}

class FacetGlobalConnectivityAccessor : public DataAccessor<Element> {
public:
  FacetGlobalConnectivityAccessor(Mesh & mesh)
      : global_connectivity("global_connectivity",
                            "facet_connectivity_synchronizer") {
    global_connectivity.initialize(
        mesh, _spatial_dimension = _all_dimensions, _with_nb_element = true,
        _with_nb_nodes_per_element = true, _element_kind = _ek_regular);
    mesh.getGlobalConnectivity(global_connectivity);
  }

  [[nodiscard]] Int getNbData(const Array<Element> & elements,
                              const SynchronizationTag & tag) const override {
    Int size = 0;
    if (tag == SynchronizationTag::_smmc_facets_conn) {
      Int nb_nodes = Mesh::getNbNodesPerElementList(elements);
      size += nb_nodes * Int(sizeof(Idx));
    }
    return size;
  }

  void packData(CommunicationBuffer & buffer, const Array<Element> & elements,
                const SynchronizationTag & tag) const override {
    if (tag == SynchronizationTag::_smmc_facets_conn) {
      for (const auto & element : elements) {
        const auto & conns =
            global_connectivity(element.type, element.ghost_type);
        for (auto n : arange(conns.getNbComponent())) {
          buffer << conns(element.element, n);
        }
      }
    }
  }

  void unpackData(CommunicationBuffer & buffer, const Array<Element> & elements,
                  const SynchronizationTag & tag) override {
    if (tag == SynchronizationTag::_smmc_facets_conn) {
      for (const auto & element : elements) {
        auto & conns = global_connectivity(element.type, element.ghost_type);
        for (auto n : arange(conns.getNbComponent())) {
          buffer >> conns(element.element, n);
        }
      }
    }
  }

  AKANTU_GET_MACRO(GlobalConnectivity, (global_connectivity), decltype(auto));

private:
  ElementTypeMapArray<Idx> global_connectivity;
};

/* -------------------------------------------------------------------------- */
// const Array<Real> & Mesh::getNormals(ElementType element_type,
//                                      GhostType ghost_type) {
//   if (this->hasData<Real>("normals", element_type, ghost_type)) {
//     return this->getData<Real>("normals", element_type, ghost_type);
//   }

//   auto & normals = getDataPointer<Real>("normals", element_type, ghost_type,
//                                         spatial_dimension, true);
//   for (auto && data [[gnu::unused]] :
//        enumerate(make_view(normals, spatial_dimension))) {
//     AKANTU_TO_IMPLEMENT();
//   }

//   AKANTU_TO_IMPLEMENT();
// }

/* -------------------------------------------------------------------------- */
Mesh & Mesh::initMeshFacets(const ID & id) {
  AKANTU_DEBUG_IN();

  if (mesh_facets) {
    AKANTU_DEBUG_OUT();
    return *mesh_facets;
  }

  mesh_facets = std::make_unique<Mesh>(spatial_dimension, this->nodes,
                                       getID() + ":" + id);
  mesh_facets->mesh_parent = this;
  mesh_facets->is_mesh_facets = true;
  mesh_facets->nodes_flags = this->nodes_flags;
  mesh_facets->nodes_global_ids = this->nodes_global_ids;

  MeshUtils::buildAllFacets(*this, *mesh_facets, 0);

  if (mesh.isDistributed()) {
    mesh_facets->is_distributed = true;
    mesh_facets->element_synchronizer = std::make_unique<FacetSynchronizer>(
        *mesh_facets, mesh.getElementSynchronizer());

    FacetGlobalConnectivityAccessor data_accessor(*mesh_facets);
    /// communicate
    mesh_facets->element_synchronizer->synchronizeOnce(
        data_accessor, SynchronizationTag::_smmc_facets_conn);

    /// flip facets
    MeshUtils::flipFacets(*mesh_facets, data_accessor.getGlobalConnectivity(),
                          _ghost);
  }

  /// transfers the the mesh physical names to the mesh facets
  if (not this->hasData("physical_names")) {
    AKANTU_DEBUG_OUT();
    return *mesh_facets;
  }

  auto & mesh_phys_data = this->getData<std::string>("physical_names");
  auto & mesh_to_mesh_facet = this->getData<Element>("mesh_to_mesh_facet");
  mesh_to_mesh_facet.initialize(*this,
                                _spatial_dimension = spatial_dimension - 1,
                                _with_nb_element = true);

  auto & phys_data = mesh_facets->getData<std::string>("physical_names");
  phys_data.initialize(*mesh_facets, _spatial_dimension = spatial_dimension - 1,
                       _with_nb_element = true);

  ElementTypeMapArray<Real> barycenters(getID(), "temporary_barycenters");
  barycenters.initialize(*mesh_facets, _nb_component = spatial_dimension,
                         _spatial_dimension = spatial_dimension - 1,
                         _with_nb_element = true);

  for (auto && ghost_type : ghost_types) {
    for (auto && type :
         barycenters.elementTypes(spatial_dimension - 1, ghost_type)) {
      mesh_facets->getBarycenters(barycenters(type, ghost_type), type,
                                  ghost_type);
    }
  }

  for_each_element(
      mesh,
      [&](auto && element) {
        Vector<Real> barycenter(spatial_dimension);
        mesh.getBarycenter(element, barycenter);
        auto norm_barycenter = barycenter.norm();
        auto tolerance = Math::getTolerance();
        if (norm_barycenter > tolerance) {
          tolerance *= norm_barycenter;
        }

        Vector<Real> barycenter_facet(spatial_dimension);

        auto range = enumerate(make_view(
            barycenters(element.type, element.ghost_type), spatial_dimension));
#ifndef AKANTU_NDEBUG
        auto min_dist = std::numeric_limits<Real>::max();
#endif
        // this is a spacial search coded the most inefficient way.
        auto facet =
            std::find_if(range.begin(), range.end(), [&](auto && data) {
              auto norm_distance = barycenter.distance(std::get<1>(data));
#ifndef AKANTU_NDEBUG
              min_dist = std::min(min_dist, norm_distance);
#endif
              return (norm_distance < tolerance);
            });

        if (facet == range.end()) {
          AKANTU_DEBUG_INFO("The element "
                            << element
                            << " did not find its associated facet in the "
                               "mesh_facets! Try to decrease math tolerance. "
                               "The closest element was at a distance of "
                            << min_dist);
          return;
        }

        // set physical name
        auto && facet_element =
            Element{element.type, std::get<0>(*facet), element.ghost_type};
        phys_data(facet_element) = mesh_phys_data(element);

        mesh_to_mesh_facet(element) = facet_element;
      },
      _spatial_dimension = spatial_dimension - 1);

  mesh_facets->createGroupsFromMeshData<std::string>("physical_names");

  AKANTU_DEBUG_OUT();
  return *mesh_facets;
}

/* -------------------------------------------------------------------------- */
void Mesh::defineMeshParent(const Mesh & mesh) {
  AKANTU_DEBUG_IN();

  this->mesh_parent = &mesh;
  this->is_mesh_facets = true;

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Mesh::read(const std::string & filename, const MeshIOType & mesh_io_type) {

  AKANTU_DEBUG_ASSERT(not is_distributed,
                      "You cannot read a mesh that is already distributed");

  MeshIO::read(filename, *this, mesh_io_type);

  auto types =
      this->elementTypes(spatial_dimension, _not_ghost, _ek_not_defined);
  auto it = types.begin();
  auto last = types.end();
  if (it == last) {
    AKANTU_DEBUG_WARNING(
        "The mesh contained in the file "
        << filename << " does not seem to be of the good dimension."
        << " No element of dimension " << spatial_dimension << " were read.");
  }

  this->makeReady();
}

/* -------------------------------------------------------------------------- */
void Mesh::write(const std::string & filename,
                 const MeshIOType & mesh_io_type) {
  MeshIO::write(filename, *this, mesh_io_type);
}

/* -------------------------------------------------------------------------- */
void Mesh::makeReady() {
  this->nb_global_nodes = this->nodes->size();
  this->computeBoundingBox();
  this->nodes_flags->resize(nodes->size(), NodeFlag::_normal);
  this->nodes_to_elements.resize(nodes->size());
  for (auto & node_set : nodes_to_elements) {
    node_set = std::make_unique<std::set<Element>>();
  }
}

/* -------------------------------------------------------------------------- */
void Mesh::printself(std::ostream & stream, int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "Mesh [\n";
  stream << space << " + id                : " << getID() << "\n";
  stream << space << " + spatial dimension : " << this->spatial_dimension
         << "\n";
  stream << space << " + nodes [\n";
  nodes->printself(stream, indent + 2);
  stream << space << " + connectivities [\n";
  connectivities.printself(stream, indent + 2);
  stream << space << " ]\n";

  GroupManager::printself(stream, indent + 1);
  stream << space << "]\n";
}

/* -------------------------------------------------------------------------- */
void Mesh::computeBoundingBox() {
  AKANTU_DEBUG_IN();

  bbox_local.reset();

  for (auto & pos : make_view(*nodes, spatial_dimension)) {
    //    if(!isPureGhostNode(i))
    bbox_local += pos;
  }

  if (this->is_distributed) {
    bbox = bbox_local.allSum(*communicator);
  } else {
    bbox = bbox_local;
  }

  size = bbox.size();

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void Mesh::getGlobalConnectivity(
    ElementTypeMapArray<Idx> & global_connectivity) {
  AKANTU_DEBUG_IN();

  for (auto && ghost_type : ghost_types) {
    for (auto type :
         global_connectivity.elementTypes(_spatial_dimension = _all_dimensions,
         _element_kind = _ek_not_defined, _ghost_type = ghost_type)) {
      if (not connectivities.exists(type, ghost_type)) {
        continue;
      }

      auto local_conn_view = make_view(connectivities(type, ghost_type));
      auto global_conn_view = make_view(global_connectivity(type, ghost_type));

      std::transform(local_conn_view.begin(), local_conn_view.end(),
                     global_conn_view.begin(),
                     [&](Idx l) -> Idx { return this->getNodeGlobalId(l); });
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
DumperIOHelper & Mesh::getGroupDumper(const std::string & dumper_name,
                                      const std::string & group_name) {
  if (group_name == "all") {
    return this->getDumper(dumper_name);
  }
  return element_groups[group_name]->getDumper(dumper_name);
}

/* -------------------------------------------------------------------------- */
template <typename T>
ElementTypeMap<Int> Mesh::getNbDataPerElem(ElementTypeMapArray<T> & arrays) {
  ElementTypeMap<Int> nb_data_per_elem;

  for (auto type : arrays.elementTypes(_element_kind = _ek_not_defined)) {
    auto nb_elements = this->getNbElement(type);
    auto & array = arrays(type);

    nb_data_per_elem(type) = array.getNbComponent() * array.size();
    nb_data_per_elem(type) /= nb_elements;
  }

  return nb_data_per_elem;
}

/* -------------------------------------------------------------------------- */
template ElementTypeMap<Int>
Mesh::getNbDataPerElem(ElementTypeMapArray<Real> & array);

template ElementTypeMap<Int>
Mesh::getNbDataPerElem(ElementTypeMapArray<Int> & array);

/* -------------------------------------------------------------------------- */
template <typename T>
std::shared_ptr<dumpers::Field>
Mesh::createFieldFromAttachedData(const std::string & field_id,
                                  const std::string & group_name,
                                  ElementKind element_kind) {

  std::shared_ptr<dumpers::Field> field;
  ElementTypeMapArray<T> * internal = nullptr;
  try {
    internal = &(this->getData<T>(field_id));
  } catch (...) {
    return nullptr;
  }

  auto && nb_data_per_elem = this->getNbDataPerElem(*internal);

  field = this->createElementalField<T, dumpers::InternalMaterialField>(
      *internal, group_name, this->spatial_dimension, element_kind,
      nb_data_per_elem);

  return field;
}

template std::shared_ptr<dumpers::Field>
Mesh::createFieldFromAttachedData<Real>(const std::string & field_id,
                                        const std::string & group_name,
                                        ElementKind element_kind);

template std::shared_ptr<dumpers::Field>
Mesh::createFieldFromAttachedData<Int>(const std::string & field_id,
                                       const std::string & group_name,
                                       ElementKind element_kind);

/* -------------------------------------------------------------------------- */
void Mesh::distribute(std::shared_ptr<MeshPartition> partition,
                      Communicator & communicator) {
  AKANTU_DEBUG_ASSERT(is_distributed == false,
                      "This mesh is already distribute");
  this->communicator = &communicator;

  this->element_synchronizer = std::make_unique<ElementSynchronizer>(
      *this, this->getID() + ":element_synchronizer", true);

  this->node_synchronizer = std::make_unique<NodeSynchronizer>(
      *this, this->getID() + ":node_synchronizer", true);

  auto psize = this->communicator->getNbProc();

  if (psize > 1) {
    auto prank = this->communicator->whoAmI();
    if (prank == 0) {
      MeshUtilsDistribution::distributeMeshCentralized(*this, 0, *partition);
    } else {
      MeshUtilsDistribution::distributeMeshCentralized(*this, 0);
    }
  }

  // if (psize > 1)
  this->is_distributed = true;

  this->computeBoundingBox();

  MeshIsDistributedEvent event(AKANTU_CURRENT_FUNCTION);
  this->sendEvent(event);
}

/* -------------------------------------------------------------------------- */
void Mesh::distributeImpl(
    Communicator & communicator,
    const std::function<Int(const Element &, const Element &)> &
        edge_weight_function [[gnu::unused]],
    const std::function<Int(const Element &)> & vertex_weight_function
    [[gnu::unused]]) {
  AKANTU_DEBUG_ASSERT(is_distributed == false,
                      "This mesh is already distribute");
  auto psize = communicator.getNbProc();

  if (psize > 1) {
#ifdef AKANTU_USE_SCOTCH
    std::shared_ptr<MeshPartition> partition;
    auto prank = communicator.whoAmI();
    if (prank == 0) {
      partition =
          std::make_shared<MeshPartitionScotch>(*this, spatial_dimension);
      partition->partitionate(psize, edge_weight_function,
                              vertex_weight_function);
    }
    this->distribute(partition, communicator);
#else
    AKANTU_ERROR("Cannot distribute a mesh without a partitioning tool");
#endif
  }
}

/* -------------------------------------------------------------------------- */
void Mesh::getAssociatedElements(const Array<Idx> & node_list,
                                 Array<Element> & elements) {
  for (const auto & node : node_list) {
    for (const auto & element : *nodes_to_elements[node]) {
      elements.push_back(element);
    }
  }
}

/* -------------------------------------------------------------------------- */
void Mesh::getAssociatedElements(const Idx & node,
                                 Array<Element> & elements) const {
  for (const auto & element : *nodes_to_elements[node]) {
    elements.push_back(element);
  }
}

/* -------------------------------------------------------------------------- */
void Mesh::fillNodesToElements(Int dimension) {
  Element element{ElementNull};

  auto nb_nodes = nodes->size();
  this->nodes_to_elements.resize(nb_nodes);

  for (Int n = 0; n < nb_nodes; ++n) {
    if (this->nodes_to_elements[n]) {
      this->nodes_to_elements[n]->clear();
    } else {
      this->nodes_to_elements[n] = std::make_unique<std::set<Element>>();
    }
  }

  for (auto ghost_type : ghost_types) {
    element.ghost_type = ghost_type;
    for (const auto & type :
         elementTypes(dimension, ghost_type, _ek_not_defined)) {
      element.type = type;

      auto nb_element = this->getNbElement(type, ghost_type);
      auto connectivity = connectivities(type, ghost_type);
      auto conn_it = connectivity.begin(connectivity.getNbComponent());

      for (Int el = 0; el < nb_element; ++el, ++conn_it) {
        element.element = el;
        const auto & conn = *conn_it;
        for (auto node : conn) {
          nodes_to_elements[node]->insert(element);
        }
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
std::tuple<Idx, Idx> Mesh::updateGlobalData(NewNodesEvent & nodes_event,
                                            NewElementsEvent & elements_event) {
  if (global_data_updater) {
    return this->global_data_updater->updateData(nodes_event, elements_event);
  }
  return std::make_tuple(nodes_event.getList().size(),
                         elements_event.getList().size());
}

/* -------------------------------------------------------------------------- */
void Mesh::registerGlobalDataUpdater(
    std::unique_ptr<MeshGlobalDataUpdater> && global_data_updater) {
  this->global_data_updater = std::move(global_data_updater);
}

/* -------------------------------------------------------------------------- */
void Mesh::eraseElements(const Array<Element> & elements) {
  ElementTypeMap<Idx> last_element;

  RemovedElementsEvent event(*this, "new_numbering", AKANTU_CURRENT_FUNCTION);
  auto & remove_list = event.getList();
  auto & new_numbering = event.getNewNumbering();

  for (auto && el : elements) {
    if (el.ghost_type != _not_ghost) {
      auto & count = ghosts_counters(el);
      --count;
      AKANTU_DEBUG_ASSERT(count >= 0,
                          "Something went wrong in the ghost element counter");

      if (count > 0) {
        continue;
      }
    }
    remove_list.push_back(el);
    if (not new_numbering.exists(el.type, el.ghost_type)) {
      auto nb_element = mesh.getNbElement(el.type, el.ghost_type);
      auto & numbering =
          new_numbering.alloc(nb_element, 1, el.type, el.ghost_type);
      for (auto && [i, numb] : enumerate(numbering)) {
        numb = i;
      }
    }

    new_numbering(el) = -1;
  }

  auto find_last_not_deleted = [](auto && array, Int start) -> Int {
    do {
      --start;
    } while (start >= 0 and array[start] == -1);

    return start;
  };

  auto find_first_deleted = [](auto && array, Int start) -> Int {
    auto begin = array.begin();
    auto it = std::find_if(begin + start, array.end(),
                           [](auto & el) { return el == -1; });
    return Int(it - begin);
  };

  for (auto ghost_type : ghost_types) {
    for (auto type : new_numbering.elementTypes(_ghost_type = ghost_type)) {
      auto & numbering = new_numbering(type, ghost_type);
      auto last_not_delete = find_last_not_deleted(numbering, numbering.size());

      if (last_not_delete < 0) {
        continue;
      }

      auto pos = find_first_deleted(numbering, 0);

      while (pos < last_not_delete) {
        std::swap(numbering[pos], numbering[last_not_delete]);
        last_not_delete = find_last_not_deleted(numbering, last_not_delete);
        pos = find_first_deleted(numbering, pos + 1);
      }
    }
  }

  this->ghosts_counters.onElementsRemoved(new_numbering);
  this->sendEvent(event);
}

} // namespace akantu
