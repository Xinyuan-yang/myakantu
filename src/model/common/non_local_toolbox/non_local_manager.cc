/**
 * Copyright (©) 2012-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "non_local_manager.hh"
#include "aka_array.hh"
#include "aka_common.hh"
#include "base_weight_function.hh"
#include "grid_synchronizer.hh"
#include "integrator.hh"
#include "mesh.hh"
#include "model.hh"
#include "non_local_neighborhood.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <utility>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
NonLocalManager::NonLocalManager(Model & model,
                                 NonLocalManagerCallback & callback,
                                 const ID & id)
    : Parsable(ParserType::_neighborhoods, id),
      spatial_dimension(model.getMesh().getSpatialDimension()), id(id),
      model(model),
      integration_points_positions("integration_points_positions", id),
      volumes("volumes", id), compute_stress_calls(0), dummy_registry(nullptr),
      dummy_grid(nullptr) {
  /// parse the neighborhood information from the input file
  const Parser & parser = getStaticParser();

  /// iterate over all the non-local sections and store them in a map
  std::pair<Parser::const_section_iterator, Parser::const_section_iterator>
      weight_sect = parser.getSubSections(ParserType::_non_local);
  Parser::const_section_iterator it = weight_sect.first;
  for (; it != weight_sect.second; ++it) {
    const ParserSection & section = *it;
    ID name = section.getName();
    this->weight_function_types[name] = section;
  }

  this->callback = &callback;
}

/* -------------------------------------------------------------------------- */
NonLocalManager::~NonLocalManager() = default;

/* -------------------------------------------------------------------------- */
void NonLocalManager::initialize() {
  volumes.initialize(this->model.getFEEngine(),
                     _spatial_dimension = spatial_dimension);

  AKANTU_DEBUG_ASSERT(this->callback,
                      "A callback should be registered prior to this call");
  this->callback->insertIntegrationPointsInNeighborhoods(_not_ghost);

  auto & mesh = this->model.getMesh();
  mesh.registerEventHandler(*this, _ehp_non_local_manager);

  /// exchange the missing ghosts for the non-local neighborhoods
  this->createNeighborhoodSynchronizers();

  /// insert the ghost quadrature points of the non-local materials into the
  /// non-local neighborhoods
  this->callback->insertIntegrationPointsInNeighborhoods(_ghost);

  FEEngine & fee = this->model.getFEEngine();
  this->updatePairLists();

  /// cleanup the unneccessary ghost elements
  this->cleanupExtraGhostElements(); // nb_ghost_protected);

  this->callback->initializeNonLocal();

  this->setJacobians(fee, _ek_regular);

  this->initNonLocalVariables();
  this->computeWeights();
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::setJacobians(const FEEngine & fe_engine,
                                   ElementKind /*kind*/) {
  jacobians = &(fe_engine.getIntegratorInterface().getJacobians());
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::createNeighborhood(const ID & weight_func,
                                         const ID & neighborhood_id) {
  auto weight_func_it = this->weight_function_types.find(weight_func);
  AKANTU_DEBUG_ASSERT(weight_func_it != weight_function_types.end(),
                      "No info found in the input file for the weight_function "
                          << weight_func << " in the neighborhood "
                          << neighborhood_id);

  const ParserSection & section = weight_func_it->second;
  const ID weight_func_type = section.getOption();

  /// create new neighborhood for given ID
  neighborhoods[neighborhood_id] =
      NonLocalNeighborhoodFactory::getInstance().allocate(
          weight_func_type, neighborhood_id, *this,
          this->integration_points_positions,
          id + ":neighborhood:" + neighborhood_id);

  neighborhoods[neighborhood_id]->parseSection(section);
  neighborhoods[neighborhood_id]->initNeighborhood();
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::createNeighborhoodSynchronizers() {
  /// exchange all the neighborhood IDs, so that every proc knows how many
  /// neighborhoods exist globally
  /// First: Compute locally the maximum ID size
  Int max_id_size = 0;
  for (const auto & id : make_keys_adaptor(neighborhoods)) {
    max_id_size = std::max(max_id_size, Int(id.size()));
  }

  /// get the global maximum ID size on each proc
  const auto & communicator = model.getMesh().getCommunicator();
  communicator.allReduce(max_id_size, SynchronizerOperation::_max);

  /// get the rank for this proc and the total nb proc
  auto prank = communicator.whoAmI();
  auto psize = communicator.getNbProc();

  /// exchange the number of neighborhoods on each proc
  Array<Int> nb_neighborhoods_per_proc(psize);
  nb_neighborhoods_per_proc(prank) = Int(neighborhoods.size());
  communicator.allGather(nb_neighborhoods_per_proc);

  /// compute the total number of neighborhoods
  Int nb_neighborhoods_global = std::accumulate(
      nb_neighborhoods_per_proc.begin(), nb_neighborhoods_per_proc.end(), 0);

  /// allocate an array of chars to store the names of all neighborhoods
  Array<char> buffer(nb_neighborhoods_global, max_id_size);

  /// starting index on this proc
  Int starting_index =
      std::accumulate(nb_neighborhoods_per_proc.begin(),
                      nb_neighborhoods_per_proc.begin() + prank, 0);

  /// store the names of local neighborhoods in the buffer
  for (auto && data : enumerate(make_keys_adaptor(neighborhoods))) {
    Int c = 0;
    Idx i = Idx(std::get<0>(data));
    const auto & id = std::get<1>(data);
    for (; c < Int(id.size()); ++c) {
      buffer(i + starting_index, c) = id[c];
    }

    for (; c < max_id_size; ++c) {
      buffer(i + starting_index, c) = char(0);
    }
  }

  /// store the nb of data to send in the all gather
  Array<Int> buffer_size(nb_neighborhoods_per_proc);
  buffer_size *= max_id_size;
  /// exchange the names of all the neighborhoods with all procs
  communicator.allGatherV(buffer, buffer_size);

  for (Int i = 0; i < nb_neighborhoods_global; ++i) {
    std::stringstream neighborhood_id;
    for (Int c = 0; c < max_id_size; ++c) {
      if (buffer(i, c) == char(0)) {
        break;
      }
      neighborhood_id << buffer(i, c);
    }
    global_neighborhoods.insert(neighborhood_id.str());
  }

  /// this proc does not know all the neighborhoods -> create dummy
  /// grid so that this proc can participate in the all gather for
  /// detecting the overlap of neighborhoods this proc doesn't know
  Vector<Real> grid_center(this->spatial_dimension);
  grid_center.fill(std::numeric_limits<Real>::max());
  Vector<Real> spacing(this->spatial_dimension);
  spacing.fill(0.);

  dummy_grid = std::make_unique<SpatialGrid<IntegrationPoint>>(
      this->spatial_dimension, spacing, grid_center);

  for (const auto & neighborhood_id : global_neighborhoods) {
    auto it = neighborhoods.find(neighborhood_id);
    if (it != neighborhoods.end()) {
      it->second->createGridSynchronizer();
    } else {
      dummy_synchronizers[neighborhood_id] = std::make_unique<GridSynchronizer>(
          this->model.getMesh(), *dummy_grid,
          std::string(this->id + ":" + neighborhood_id + ":grid_synchronizer"),
          false);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::synchronize(DataAccessor<Element> & data_accessor,
                                  const SynchronizationTag & tag) {
  for (const auto & neighborhood_id : global_neighborhoods) {
    auto it = neighborhoods.find(neighborhood_id);
    if (it != neighborhoods.end()) {
      it->second->synchronize(data_accessor, tag);
    } else {
      auto synchronizer_it = dummy_synchronizers.find(neighborhood_id);
      if (synchronizer_it == dummy_synchronizers.end()) {
        continue;
      }

      synchronizer_it->second->synchronizeOnce(data_accessor, tag);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::averageInternals(GhostType ghost_type) {
  /// update the weights of the weight function
  if (ghost_type == _not_ghost) {
    this->computeWeights();
  }

  /// loop over all neighborhoods and compute the non-local variables
  for (auto & neighborhood : neighborhoods) {
    /// loop over all the non-local variables of the given neighborhood
    for (auto & non_local_variable : non_local_variables) {
      NonLocalVariable & non_local_var = *non_local_variable.second;
      neighborhood.second->weightedAverageOnNeighbours(
          non_local_var.local, non_local_var.non_local,
          non_local_var.nb_component, ghost_type);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::computeWeights() {
  AKANTU_DEBUG_IN();

  this->updateWeightFunctionInternals();
  this->volumes.zero();

  for (const auto & global_neighborhood : global_neighborhoods) {
    auto it = neighborhoods.find(global_neighborhood);

    if (it != neighborhoods.end()) {
      it->second->updateWeights();
    } else {
      dummy_synchronizers[global_neighborhood]->synchronize(
          dummy_accessor, SynchronizationTag::_mnl_weight);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::updatePairLists() {
  AKANTU_DEBUG_IN();

  integration_points_positions.initialize(
      this->model.getFEEngine(), _nb_component = spatial_dimension,
      _spatial_dimension = spatial_dimension);

  /// compute the position of the quadrature points
  this->model.getFEEngine().computeIntegrationPointsCoordinates(
      integration_points_positions);

  for (auto & neighborhood : make_values_adaptor(neighborhoods)) {
    neighborhood->updatePairList();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::registerNonLocalVariable(const ID & variable_name,
                                               const ID & nl_variable_name,
                                               UInt nb_component) {

  AKANTU_DEBUG_IN();

  auto non_local_variable_it = non_local_variables.find(variable_name);

  if (non_local_variable_it == non_local_variables.end()) {
    non_local_variables[nl_variable_name] = std::make_unique<NonLocalVariable>(
        variable_name, nl_variable_name, this->id, nb_component);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
ElementTypeMapReal &
NonLocalManager::registerWeightFunctionInternal(const ID & field_name) {
  auto it = this->weight_function_internals.find(field_name);
  if (it == weight_function_internals.end()) {
    weight_function_internals[field_name] =
        std::make_unique<ElementTypeMapReal>(field_name, this->id);
  }

  return *(weight_function_internals[field_name]);
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::updateWeightFunctionInternals() {
  for (auto & pair : this->weight_function_internals) {
    auto & internals = *pair.second;
    internals.zero();
    for (auto ghost_type : ghost_types) {
      this->callback->updateLocalInternal(internals, ghost_type, _ek_regular);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::initNonLocalVariables() {
  /// loop over all the non-local variables
  for (auto & pair : non_local_variables) {
    auto & variable = *pair.second;
    variable.non_local.initialize(this->model.getFEEngine(),
                                  _nb_component = variable.nb_component,
                                  _spatial_dimension = spatial_dimension);
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::computeAllNonLocalContribution() {

  /// update the flattened version of the internals
  for (auto & pair : non_local_variables) {
    auto & variable = *pair.second;
    variable.local.zero();
    variable.non_local.zero();
    for (auto ghost_type : ghost_types) {
      this->callback->updateLocalInternal(variable.local, ghost_type,
                                          _ek_regular);
    }
  }

  this->volumes.zero();

  for (auto & pair : neighborhoods) {
    auto & neighborhood = *pair.second;
    neighborhood.asynchronousSynchronize(SynchronizationTag::_mnl_for_average);
  }

  this->averageInternals(_not_ghost);

  AKANTU_DEBUG_INFO("Wait distant non local stresses");

  for (auto & pair : neighborhoods) {
    auto & neighborhood = *pair.second;
    neighborhood.waitEndSynchronize(SynchronizationTag::_mnl_for_average);
  }

  this->averageInternals(_ghost);

  /// copy the results in the materials
  for (auto & pair : non_local_variables) {
    auto & variable = *pair.second;
    for (auto ghost_type : ghost_types) {
      this->callback->updateNonLocalInternal(variable.non_local, ghost_type,
                                             _ek_regular);
    }
  }

  this->callback->computeNonLocalContribution(_not_ghost);

  ++this->compute_stress_calls;
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::cleanupExtraGhostElements() {
  using ElementSet = std::set<Element>;
  ElementSet relevant_ghost_elements;

  /// loop over all the neighborhoods and get their protected ghosts
  for (auto & pair : neighborhoods) {
    auto & neighborhood = *pair.second;
    ElementSet to_keep_per_neighborhood;

    neighborhood.getRelevantGhostElements(to_keep_per_neighborhood);
    relevant_ghost_elements.insert(to_keep_per_neighborhood.begin(),
                                   to_keep_per_neighborhood.end());
  }

  for (auto & pair : neighborhoods) {
    auto & neighborhood = *pair.second;
    neighborhood.cleanupExtraGhostElements(relevant_ghost_elements);
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::onElementsRemoved(
    const Array<Element> & element_list,
    const ElementTypeMapArray<Idx> & new_numbering,
    const RemovedElementsEvent & event) {

  FEEngine & fee = this->model.getFEEngine();
  NonLocalManager::removeIntegrationPointsFromMap(
      event.getNewNumbering(), spatial_dimension, integration_points_positions,
      fee, _ek_regular);
  NonLocalManager::removeIntegrationPointsFromMap(event.getNewNumbering(), 1,
                                                  volumes, fee, _ek_regular);

  /// loop over all the neighborhoods and call onElementsRemoved
  auto global_neighborhood_it = global_neighborhoods.begin();
  NeighborhoodMap::iterator it;
  for (; global_neighborhood_it != global_neighborhoods.end();
       ++global_neighborhood_it) {
    it = neighborhoods.find(*global_neighborhood_it);
    if (it != neighborhoods.end()) {
      it->second->onElementsRemoved(element_list, new_numbering, event);
    } else {
      dummy_synchronizers[*global_neighborhood_it]->onElementsRemoved(
          element_list, new_numbering, event);
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::onElementsAdded(const Array<Element> & /*unused*/,
                                      const NewElementsEvent & /*unused*/) {
  this->resizeElementTypeMap(1, volumes, model.getFEEngine());
  this->resizeElementTypeMap(spatial_dimension, integration_points_positions,
                             model.getFEEngine());
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::resizeElementTypeMap(Int nb_component,
                                           ElementTypeMapReal & element_map,
                                           const FEEngine & fee,
                                           const ElementKind el_kind) {
  auto & mesh = this->model.getMesh();

  for (auto gt : ghost_types) {
    for (auto type : mesh.elementTypes(spatial_dimension, gt, el_kind)) {
      auto nb_element = mesh.getNbElement(type, gt);
      auto nb_quads = fee.getNbIntegrationPoints(type, gt);
      if (!element_map.exists(type, gt)) {
        element_map.alloc(nb_element * nb_quads, nb_component, type, gt);
      } else {
        element_map(type, gt).resize(nb_element * nb_quads);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::removeIntegrationPointsFromMap(
    const ElementTypeMapArray<Idx> & new_numbering, Int nb_component,
    ElementTypeMapReal & element_map, const FEEngine & fee,
    const ElementKind el_kind) {

  for (auto gt : ghost_types) {
    for (auto type : new_numbering.elementTypes(_all_dimensions, gt, el_kind)) {
      if (element_map.exists(type, gt)) {
        const auto & renumbering = new_numbering(type, gt);

        auto & vect = element_map(type, gt);
        auto nb_quad_per_elem = fee.getNbIntegrationPoints(type, gt);
        Array<Real> tmp(renumbering.size() * nb_quad_per_elem, nb_component);

        AKANTU_DEBUG_ASSERT(
            tmp.size() == vect.size(),
            "Something strange append some mater was created or disappeared in "
                << vect.getID() << "(" << vect.size() << "!=" << tmp.size()
                << ") "
                   "!!");

        auto it = make_view(vect, nb_component * nb_quad_per_elem).begin();
        auto tmp_it = make_view(tmp, nb_component * nb_quad_per_elem).begin();
        Int new_size = 0;
        for (Int i = 0; i < renumbering.size(); ++i) {
          auto new_i = renumbering(i);
          if (new_i != Int(-1)) {
            tmp_it[new_i] = it[i];
            ++new_size;
          }
        }
        tmp.resize(new_size * nb_quad_per_elem);
        vect.copy(tmp);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
Int NonLocalManager::getNbData(const Array<Element> & elements,
                               const ID & id) const {
  Int size = 0;
  auto nb_quadrature_points = this->model.getNbIntegrationPoints(elements);
  auto it = non_local_variables.find(id);

  AKANTU_DEBUG_ASSERT(it != non_local_variables.end(),
                      "The non-local variable " << id << " is not registered");

  size += it->second->nb_component * Int(sizeof(Real)) * nb_quadrature_points;
  return size;
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::packData(CommunicationBuffer & buffer,
                               const Array<Element> & elements,
                               const ID & id) const {
  auto it = non_local_variables.find(id);

  AKANTU_DEBUG_ASSERT(it != non_local_variables.end(),
                      "The non-local variable " << id << " is not registered");

  DataAccessor<Element>::packElementalDataHelper<Real>(
      it->second->local, buffer, elements, this->model.getFEEngine());
}

/* -------------------------------------------------------------------------- */
void NonLocalManager::unpackData(CommunicationBuffer & buffer,
                                 const Array<Element> & elements,
                                 const ID & id) const {
  auto it = non_local_variables.find(id);

  AKANTU_DEBUG_ASSERT(it != non_local_variables.end(),
                      "The non-local variable " << id << " is not registered");

  DataAccessor<Element>::unpackElementalDataHelper<Real>(
      it->second->local, buffer, elements, this->model.getFEEngine());
}

} // namespace akantu
