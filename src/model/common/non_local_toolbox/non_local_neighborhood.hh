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
#include "aka_common.hh"
#include "neighborhood_base.hh"
#include "non_local_neighborhood_base.hh"
#include "parsable.hh"
/* -------------------------------------------------------------------------- */
#include <map>
#include <memory>
#include <string>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_NON_LOCAL_NEIGHBORHOOD_HH
#define AKANTU_NON_LOCAL_NEIGHBORHOOD_HH

namespace akantu {
class BaseWeightFunction;
class NonLocalManager;
} // namespace akantu

namespace akantu {

template <class WeightFunction = BaseWeightFunction>
class NonLocalNeighborhood : public NonLocalNeighborhoodBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  NonLocalNeighborhood(NonLocalManager & manager,
                       const ElementTypeMapReal & quad_coordinates,
                       const ID & id = "neighborhood");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  /// compute the weights for non-local averaging
  void computeWeights() override;

  /// save the pair of weights in a file
  void saveWeights(const std::string & filename) const override;

  /// compute the non-local counter part for a given element type map
  // compute the non-local counter part for a given element type map
  void weightedAverageOnNeighbours(const ElementTypeMapReal & to_accumulate,
                                   ElementTypeMapReal & accumulated,
                                   Int nb_degree_of_freedom,
                                   GhostType ghost_type2) const override;

  /// update the weights based on the weight function
  void updateWeights() override;

  /// register a new non-local variable in the neighborhood
  // void registerNonLocalVariable(const ID & id);
protected:
  template <class Func>
  inline void foreach_weight(GhostType ghost_type, Func && func);

  template <class Func>
  inline void foreach_weight(GhostType ghost_type, Func && func) const;

  [[nodiscard]] inline auto getNbData(const Array<Element> & elements,
                                      const SynchronizationTag & tag) const
      -> Int override;

  inline void packData(CommunicationBuffer & buffer,
                       const Array<Element> & elements,
                       const SynchronizationTag & tag) const override;

  inline void unpackData(CommunicationBuffer & buffer,
                         const Array<Element> & elements,
                         const SynchronizationTag & tag) override;

  /* ------------------------------------------------------------------------ */
  /* Accessor                                                                 */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_GET_MACRO_AUTO(NonLocalManager, non_local_manager);
  AKANTU_GET_MACRO_AUTO_NOT_CONST(NonLocalManager, non_local_manager);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// Pointer to non-local manager class
  NonLocalManager & non_local_manager;

  /// the weights associated to the pairs
  std::map<GhostType, std::unique_ptr<Array<Real>>> pair_weight;

  /// weight function
  std::shared_ptr<WeightFunction> weight_function;
};

} // namespace akantu

#include "non_local_neighborhood_tmpl.hh" // IWYU pragma: keep

#endif /* AKANTU_NON_LOCAL_NEIGHBORHOOD_HH */
