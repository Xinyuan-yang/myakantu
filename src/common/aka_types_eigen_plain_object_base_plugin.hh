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
// NOLINTBEGIN

#define AKANTU_EIGEN_VERSION                                                   \
  (EIGEN_WORLD_VERSION * 10000 + EIGEN_MAJOR_VERSION * 1000 +                  \
   EIGEN_MAJOR_VERSION)
template <bool _is_vector = IsVectorAtCompileTime, typename _Derived = Derived,
          std::enable_if_t<_is_vector> * = nullptr>
EIGEN_DEVICE_FUNC constexpr EIGEN_STRONG_INLINE
PlainObjectBase(std::initializer_list<Scalar> list) {
  static_assert(std::is_trivially_copyable<Scalar>{},
                "Cannot create a tensor on non trivial types");
  //  _check_template_params();

  if constexpr (Base::SizeAtCompileTime != 1) {
    this->template _init1<Index>(list.size());
  }

  Index i = 0;
  for (auto val : list) {
    this->operator()(i++) = val;
  }
}

#if AKANTU_EIGEN_VERSION < 34000
template <bool _is_vector = IsVectorAtCompileTime, typename _Derived = Derived,
          std::enable_if_t<not _is_vector> * = nullptr>
EIGEN_DEVICE_FUNC constexpr EIGEN_STRONG_INLINE
PlainObjectBase(std::initializer_list<std::initializer_list<Scalar>> list) {
  static_assert(std::is_trivially_copyable<Scalar>{},
                "Cannot create a tensor on non trivial types");
  Index m = list.size();
  Index n = 0;
  for (auto row : list) {
    n = std::max(n, Index(row.size()));
  }

  if (RowsAtCompileTime != -1 and RowsAtCompileTime != m) {
    throw std::range_error(
        "The size of the matrix does not correspond to the initializer_list");
  }

  if (ColsAtCompileTime != -1 and ColsAtCompileTime != n) {
    throw std::range_error(
        "The size of the matrix does not correspond to the initializer_list");
  }

  _check_template_params();
  this->template _init2<Index, Index>(m, n);
  this->fill(Scalar{});

  Index i = 0, j = 0;
  for (auto & row : list) {
    for (auto & val : row) {
      this->operator()(i, j++) = val;
    }
    ++i;
    j = 0;
  }
}
#endif

// NOLINTEND
