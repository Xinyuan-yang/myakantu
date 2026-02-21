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
#include "aka_types.hh"
/* -------------------------------------------------------------------------- */
#include <utility>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_MATH_H_
#define AKANTU_AKA_MATH_H_

namespace akantu {
/* -------------------------------------------------------------------------- */

namespace Math {
/// tolerance for functions that need one
extern Real tolerance; // NOLINT

/* ------------------------------------------------------------------------ */
/* Geometry                                                                 */
/* ------------------------------------------------------------------------ */
/// compute normal a normal to a vector
template <class D1, std::enable_if_t<aka::is_vector_v<D1>> * = nullptr>
inline Vector<Real> normal(const Eigen::MatrixBase<D1> & vec);

template <class D1, aka::enable_if_t<not aka::is_vector_v<D1>> * = nullptr>
inline Vector<Real> normal(const Eigen::MatrixBase<D1> & /*vec*/) {
  AKANTU_TO_IMPLEMENT();
}

/// compute normal a normal to a vector
template <class D1, class D2,
          std::enable_if_t<aka::are_vectors<D1, D2>::value> * = nullptr>
inline Vector<Real, 3> normal(const Eigen::MatrixBase<D1> & vec1,
                              const Eigen::MatrixBase<D2> & vec2);

/// compute the tangents to an array of normal vectors
void compute_tangents(const Array<Real> & normals, Array<Real> & tangents);

/// radius of the in-circle of a triangle in 2d space
template <class D1, class D2, class D3>
static inline Real triangle_inradius(const Eigen::MatrixBase<D1> & coord1,
                                     const Eigen::MatrixBase<D2> & coord2,
                                     const Eigen::MatrixBase<D3> & coord3);

/// radius of the in-circle of a tetrahedron
template <class D1, class D2, class D3, class D4>
static inline Real tetrahedron_inradius(const Eigen::MatrixBase<D1> & coord1,
                                        const Eigen::MatrixBase<D2> & coord2,
                                        const Eigen::MatrixBase<D3> & coord3,
                                        const Eigen::MatrixBase<D4> & coord4);
/// volume of a tetrahedron
template <class D1, class D2, class D3, class D4>
static inline Real tetrahedron_volume(const Eigen::MatrixBase<D1> & coord1,
                                      const Eigen::MatrixBase<D2> & coord2,
                                      const Eigen::MatrixBase<D3> & coord3,
                                      const Eigen::MatrixBase<D4> & coord4);

/// compute the barycenter of n points
template <class D1, class D2>
inline void barycenter(const Eigen::MatrixBase<D1> & coord,
                       Eigen::MatrixBase<D2> & barycenter);

/// test if two scalar are equal within a given tolerance
inline bool are_float_equal(Real x, Real y);

/// test if two vectors are equal within a given tolerance
inline bool are_vector_equal(Int n, Real * x, Real * y);

#ifdef isnan
#error                                                                         \
    "You probably  included <math.h> which  is incompatible with aka_math  please use\
<cmath> or add a \"#undef isnan\" before akantu includes"
#endif
/// test if a real is a NaN
inline bool isnan(Real x);

/// test if the line x and y intersects each other
inline bool intersects(Real x_min, Real x_max, Real y_min, Real y_max);

/// test if a is in the range [x_min, x_max]
inline bool is_in_range(Real a, Real x_min, Real x_max);

inline Real getTolerance() { return Math::tolerance; }
inline void setTolerance(Real tol) { Math::tolerance = tol; }

template <Int p, typename T> inline T pow(T x);

template <class T1, class T2,
          std::enable_if_t<std::is_integral<T1>::value and
                           std::is_integral<T2>::value> * = nullptr>
inline Real kronecker(T1 i, T2 j) {
  return static_cast<Real>(i == j);
}
/* --------------------------------------------------------------------------
 */
template <typename T> static inline constexpr T pow(T x, int p) {
  return p == 0 ? T(1) : (pow(x, p - 1) * x);
}

/// reduce all the values of an array, the summation is done in place and the
/// array is modified
Real reduce(Array<Real> & array);

template <class T> class NewtonRaphson {
public:
  NewtonRaphson(Real tolerance, Int max_iteration)
      : tolerance(tolerance), max_iteration(max_iteration) {}

  template <class Functor> T solve(const Functor & funct, const T & x_0);

private:
  Real tolerance;
  Int max_iteration;
};

template <class T> struct NewtonRaphsonFunctor {
  explicit NewtonRaphsonFunctor(const std::string & name) : name(name) {}

  virtual T f(const T & x) const = 0;
  virtual T f_prime(const T & x) const = 0;

  std::string name;
};
} // namespace Math
} // namespace akantu
/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "aka_math_tmpl.hh"

#endif /* AKANTU_AKA_MATH_H_ */
