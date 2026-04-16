/**
 * Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "aka_array.hh"
#include "aka_common.hh"
#include "communicator.hh"
#include "model.hh"
#include "non_local_manager.hh"
#include "non_local_neighborhood.hh"
/* -------------------------------------------------------------------------- */
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NON_LOCAL_NEIGHBORHOOD_TMPL_HH
#define AKANTU_NON_LOCAL_NEIGHBORHOOD_TMPL_HH

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
template <class Func>
inline void
NonLocalNeighborhood<WeightFunction>::foreach_weight(GhostType ghost_type,
                                                     Func && func) {
  auto && weight_vect = *pair_weight[ghost_type];

  for (auto && [pair, weight] :
       zip(pair_list[ghost_type],
           make_view(weight_vect, weight_vect.getNbComponent()))) {
    auto && [q1, q2] = pair;
    std::forward<decltype(func)>(func)(q1, q2, weight);
  }
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
template <class Func>
inline void
NonLocalNeighborhood<WeightFunction>::foreach_weight(GhostType ghost_type,
                                                     Func && func) const {
  auto && weight_vect = *pair_weight.at(ghost_type);

  for (auto && [pair, weight] :
       zip(pair_list.at(ghost_type),
           make_view(weight_vect, weight_vect.getNbComponent()))) {
    auto && [q1, q2] = pair;
    std::forward<decltype(func)>(func)(q1, q2, weight);
  }
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
NonLocalNeighborhood<WeightFunction>::NonLocalNeighborhood(
    NonLocalManager & manager, const ElementTypeMapReal & quad_coordinates,
    const ID & id)
    : NonLocalNeighborhoodBase(manager.getModel(), quad_coordinates, id),
      non_local_manager(manager) {
  AKANTU_DEBUG_IN();

  this->weight_function = std::make_unique<WeightFunction>(manager);

  this->registerSubSection(ParserType::_weight_function, "weight_parameter",
                           *weight_function);

  for (auto ghost_type : ghost_types) {
    if (pair_weight.find(ghost_type) == pair_weight.end()) {
      pair_weight[ghost_type] = std::make_unique<Array<Real>>(0, 2);
      pair_list[ghost_type].resize(0);
    }
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
void NonLocalNeighborhood<WeightFunction>::computeWeights() {
  AKANTU_DEBUG_IN();

  this->weight_function->setRadius(this->neighborhood_radius);

  Int nb_weights_per_pair = 2; /// w1: q1->q2, w2: q2->q1

  /// get the elementtypemap for the neighborhood volume for each quadrature
  /// point
  ElementTypeMapReal & quadrature_points_volumes =
      this->non_local_manager.getVolumes();

  /// update the internals of the weight function if applicable (not
  /// all the weight functions have internals and do noting in that
  /// case)
  weight_function->updateInternals();

  for (auto ghost_type : ghost_types) {
    /// allocate the array to store the weight, if it doesn't exist already
    if (not pair_weight[ghost_type]) {
      pair_weight[ghost_type] =
          std::make_unique<Array<Real>>(0, nb_weights_per_pair);
    }

    auto && pair_lists = pair_list[ghost_type];
    auto && pair_weights = *pair_weight[ghost_type];
    pair_weights.resize(pair_lists.size());
    pair_weights.zero();

    auto && jacobians = this->non_local_manager.getJacobians();
    /// loop over all pairs in the current pair list array and their
    /// corresponding weights
    // Compute the weights
    for (auto && [weight, pairs] :
         zip(make_view(pair_weights, nb_weights_per_pair), pair_lists)) {
      const auto & [q1, q2] = pairs;

      /// get the coordinates for the given pair of quads
      auto && q1_coord = this->quad_coordinates.get(q1);
      auto && q2_coord = this->quad_coordinates.get(q2);

      auto && quad_volumes_1 = quadrature_points_volumes(q1);
      const auto & q2_wJ = jacobians(q2);

      /// compute distance between the two quadrature points
      auto r = q1_coord.distance(q2_coord);

      /// compute the weight for averaging on q1 based on the distance
      weight(0) = q2_wJ * (*this->weight_function)(r, q1, q2);

      quad_volumes_1 += weight(0);

      if (q2.ghost_type != _ghost and q1.global_num != q2.global_num) {
        const auto & q1_wJ = jacobians(q1);
        auto && quad_volumes_2 = quadrature_points_volumes(q2);

        /// compute the weight for averaging on q2
        weight(1) = q1_wJ * (*this->weight_function)(r, q2, q1);
        quad_volumes_2 += weight(1);
      } else {
        weight(1) = 0.;
      }
    }
  }

  ///  normalize the weights
  for (auto ghost_type : ghost_types) {
    foreach_weight(ghost_type,
                   [&](const auto & q1, const auto & q2, auto & weight) {
                     auto q1_volume = quadrature_points_volumes(q1);
                     auto ghost_type2 = q2.ghost_type;
                     weight(0) *= 1. / q1_volume;
                     if (ghost_type2 != _ghost) {
                       auto q2_volume = quadrature_points_volumes(q2);
                       weight(1) *= 1. / q2_volume;
                     }
                   });
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
void NonLocalNeighborhood<WeightFunction>::saveWeights(
    const std::string & filename) const {
  std::ofstream pout;

  std::stringstream sstr;

  const Communicator & comm = model.getMesh().getCommunicator();

  auto prank = comm.whoAmI();
  sstr << filename << "." << prank;

  pout.open(sstr.str().c_str());

  for (auto ghost_type : ghost_types) {
    AKANTU_DEBUG_ASSERT((pair_weight.at(ghost_type)),
                        "the weights have not been computed yet");

    const auto & weights = *pair_weight.at(ghost_type);
    for (auto && weight : make_view(weights, 2)) {
      pout << "w1: " << weight(0) << " w2: " << weight(1) << "\n";
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
void NonLocalNeighborhood<WeightFunction>::weightedAverageOnNeighbours(
    const ElementTypeMapReal & to_accumulate, ElementTypeMapReal & accumulated,
    Int /*nb_degree_of_freedom*/, GhostType ghost_type2) const {

  auto it = non_local_variables.find(accumulated.getName());
  // do averaging only for variables registered in the neighborhood
  if (it == non_local_variables.end()) {
    return;
  }

  foreach_weight(ghost_type2,
                 [ghost_type2, &to_accumulate, &accumulated](
                     const auto & q1, const auto & q2, auto & weight) {
                   auto && to_acc_2 = to_accumulate.get(q2);
                   auto && acc_1 = accumulated.get(q1);

                   acc_1 += weight(0) * to_acc_2;

                   if (ghost_type2 != _ghost) {
                     auto && to_acc_1 = to_accumulate.get(q1);
                     auto && acc_2 = accumulated.get(q2);

                     acc_2 += weight(1) * to_acc_1;
                   }
                 });
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
void NonLocalNeighborhood<WeightFunction>::updateWeights() {
  // Update the weights for the non local variable averaging
  if (this->weight_function->getUpdateRate() and
      (this->non_local_manager.getNbStressCalls() %
           this->weight_function->getUpdateRate() ==
       0)) {
    SynchronizerRegistry::synchronize(SynchronizationTag::_mnl_weight);
    this->computeWeights();
  }
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
inline auto NonLocalNeighborhood<WeightFunction>::getNbData(
    const Array<Element> & elements, const SynchronizationTag & tag) const
    -> Int {
  Int size = 0;

  if (tag == SynchronizationTag::_mnl_for_average) {
    for (auto && variable_id : non_local_variables) {
      size += this->non_local_manager.getNbData(elements, variable_id);
    }
  }

  size += this->weight_function->getNbData(elements, tag);

  return size;
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
inline void NonLocalNeighborhood<WeightFunction>::packData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & tag) const {
  if (tag == SynchronizationTag::_mnl_for_average) {
    for (auto && variable_id : non_local_variables) {
      this->non_local_manager.packData(buffer, elements, variable_id);
    }
  }

  this->weight_function->packData(buffer, elements, tag);
}

/* -------------------------------------------------------------------------- */
template <class WeightFunction>
inline void NonLocalNeighborhood<WeightFunction>::unpackData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & tag) {
  if (tag == SynchronizationTag::_mnl_for_average) {
    for (auto & variable_id : non_local_variables) {
      this->non_local_manager.unpackData(buffer, elements, variable_id);
    }
  }

  this->weight_function->unpackData(buffer, elements, tag);
}

} // namespace akantu

#endif /* __AKANTU_NON_LOCAL_NEIGHBORHOOD_TMPL__ */
