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
#include "aka_math.hh"
#include "aka_array.hh"
#include "aka_iterators.hh"
/* -------------------------------------------------------------------------- */

namespace akantu::Math {
/* ---------------------------------------------------------------------- */
void compute_tangents(const Array<Real> & normals, Array<Real> & tangents) {
  if (normals.empty()) {
    return;
  }

  auto spatial_dimension = normals.getNbComponent();
  auto tangent_components = spatial_dimension * (spatial_dimension - 1);

  if (tangent_components == 0) {
    return;
  }

  AKANTU_DEBUG_ASSERT(
      tangent_components == tangents.getNbComponent(),
      "Cannot compute the tangents, the storage array for tangents"
          << " does not have the good amount of components.");

  auto nb_normals = normals.size();
  tangents.resize(nb_normals);
  tangents.zero();

  /// compute first tangent
  for (auto && [normal_v, tangent] :
       zip(make_view(normals, spatial_dimension),
           make_view(tangents, spatial_dimension, spatial_dimension - 1))) {

    if (are_float_equal(normal_v.norm(), 0.)) {
      tangent(0, 0) = 1.;
    } else {
      tangent(0) = normal(normal_v);
    }
  }

  /// compute second tangent (3D case)
  if (spatial_dimension == 3) {
    for (auto && [normal_v, tangent] :
         zip(make_view(normals, spatial_dimension),
             make_view(tangents, spatial_dimension, spatial_dimension - 1))) {

      tangent(1) = normal(normal_v, tangent(0));
    }
  }
}

/* ---------------------------------------------------------------------- */
Real reduce(Array<Real> & array) {
  UInt nb_values = array.size();
  if (nb_values == 0) {
    return 0.;
  }

  Int nb_values_to_sum = nb_values >> 1;

  std::sort(array.begin(), array.end());

  // as long as the half is not empty
  while (nb_values_to_sum != 0) {
    Int remaining = (nb_values - 2 * nb_values_to_sum);
    if (remaining != 0U) {
      array(nb_values - 2) += array(nb_values - 1);
    }

    // sum to consecutive values and store the sum in the first half
    for (Int i = 0; i < nb_values_to_sum; ++i) {
      array(i) = array(2 * i) + array(2 * i + 1);
    }

    nb_values = nb_values_to_sum;
    nb_values_to_sum >>= 1;
  }

  return array(0);
}
} // namespace akantu::Math
