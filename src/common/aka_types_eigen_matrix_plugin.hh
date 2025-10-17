/**
 * Copyright (©) 2014-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#define AKANTU_EIGEN_VERSION                                                   \
  (EIGEN_WORLD_VERSION * 10000 + EIGEN_MAJOR_VERSION * 1000 +                  \
   EIGEN_MAJOR_VERSION)

using size_type = Index;

template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<_is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE Matrix(std::initializer_list<Scalar> list)
    : Base(list) {}

#if AKANTU_EIGEN_VERSION <= 34000
template <bool _is_vector = IsVectorAtCompileTime,
          std::enable_if_t<not _is_vector> * = nullptr>
EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
Matrix(std::initializer_list<std::initializer_list<Scalar>> list)
    : Base(list) {}
#endif
