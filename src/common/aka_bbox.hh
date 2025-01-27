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
#include "aka_math.hh"
#include "aka_types.hh"
#include "communicator.hh"
/* -------------------------------------------------------------------------- */
#include <algorithm>
#include <limits>
#include <map>
#include <ostream>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_BBOX_HH_
#define AKANTU_AKA_BBOX_HH_

namespace akantu {

// NOLINTNEXTLINE(cppcoreguidelines-special-member-functions)
class BBox {
public:
  BBox() = default;

  BBox(Int spatial_dimension)
      : dim(spatial_dimension), lower_bounds(spatial_dimension),
        upper_bounds(spatial_dimension) {
    lower_bounds.fill(std::numeric_limits<Real>::max());
    upper_bounds.fill(-std::numeric_limits<Real>::max());
  }

  BBox(const BBox & other) = default;
  auto operator=(const BBox & other) -> BBox & = default;

  inline auto operator+=(const Vector<Real> & position) -> BBox & {
    AKANTU_DEBUG_ASSERT(
        this->dim == position.size(),
        "You are adding a point of a wrong dimension to the bounding box");

    this->empty_ = false;

    for (auto s : arange(dim)) {
      lower_bounds(s) = std::min(lower_bounds(s), position(s));
      upper_bounds(s) = std::max(upper_bounds(s), position(s));
    }
    return *this;
  }

  /* ------------------------------------------------------------------------ */
  [[nodiscard]] inline auto intersects(const BBox & other,
                                       const SpatialDirection & direction) const
      -> bool {
    AKANTU_DEBUG_ASSERT(
        this->dim == other.dim,
        "You are intersecting bounding boxes of different dimensions");
    return Math::intersects(lower_bounds(direction), upper_bounds(direction),
                            other.lower_bounds(direction),
                            other.upper_bounds(direction));
  }

  [[nodiscard]] inline auto intersects(const BBox & other) const -> bool {
    if (this->empty() or other.empty()) {
      return false;
    }

    bool intersect = true;
    for (auto s : arange(this->dim)) {
      intersect &= this->intersects(other, SpatialDirection(s));
    }
    return intersect;
  }

  /* ------------------------------------------------------------------------ */
  [[nodiscard]] inline auto intersection(const BBox & other) const -> BBox {
    AKANTU_DEBUG_ASSERT(
        this->dim == other.dim,
        "You are intersecting bounding boxes of different dimensions");

    BBox intersect(this->dim);
    intersect.empty_ = not this->intersects(other);

    if (intersect.empty()) {
      return intersect;
    }

    for (auto s : arange(this->dim)) {
      // is lower point in range ?
      // bool point1 = Math::is_in_range(other.lower_bounds(s), lower_bounds(s),
      //                                 upper_bounds(s));

      // is upper point in range ?
      // bool point2 = Math::is_in_range(other.upper_bounds(s), lower_bounds(s),
      //                                 upper_bounds(s));
      intersect.lower_bounds(s) =
          std::clamp(other.lower_bounds(s), lower_bounds(s), upper_bounds(s));
      intersect.upper_bounds(s) = -std::clamp(
          -other.upper_bounds(s), -upper_bounds(s), -lower_bounds(s));

      // if (point1 and not point2) {
      //   // |-----------|         this (i)
      //   //       |-----------|   other(i)
      //   //       1           2
      //   intersect.lower_bounds(s) = other.lower_bounds(s);
      //   intersect.upper_bounds(s) = upper_bounds(s);
      // } else if (point1 && point2) {
      //   // |-----------------|   this (i)
      //   //   |-----------|       other(i)
      //   //   1           2
      //   intersect.lower_bounds(s) = other.lower_bounds(s);
      //   intersect.upper_bounds(s) = other.upper_bounds(s);
      // } else if (!point1 && point2) {
      //   //       |-----------|   this (i)
      //   // |-----------|         other(i)
      //   // 1           2
      //   intersect.lower_bounds(s) = this->lower_bounds(s);
      //   intersect.upper_bounds(s) = other.upper_bounds(s);
      // } else {
      //   //   |-----------|       this (i)
      //   // |-----------------|   other(i)
      //   // 1                 2
      //   intersect.lower_bounds(s) = this->lower_bounds(s);
      //   intersect.upper_bounds(s) = this->upper_bounds(s);
      // }
    }

    return intersect;
  }

  /* ------------------------------------------------------------------------ */
  template <class Derived>
  inline auto contains(const Eigen::MatrixBase<Derived> & point) const -> bool {
    return (point.array() >= lower_bounds.array()).all() and
           (point.array() <= upper_bounds.array()).all();
  }

  /* ------------------------------------------------------------------------ */
  inline void reset() {
    lower_bounds.set(std::numeric_limits<Real>::max());
    upper_bounds.set(std::numeric_limits<Real>::lowest());
    empty_ = true;
  }

  /* ------------------------------------------------------------------------ */
  template <class Derived>
  inline void getCenter(Eigen::MatrixBase<Derived> & center) {
    center = (upper_bounds + lower_bounds) / 2.;
  }

  /* ------------------------------------------------------------------------ */
  [[nodiscard]] auto getLowerBounds() const -> const Vector<Real> & {
    return lower_bounds;
  }
  [[nodiscard]] auto getUpperBounds() const -> const Vector<Real> & {
    return upper_bounds;
  }

  template <typename D>
  void setLowerBounds(const Eigen::MatrixBase<D> & lower_bounds) {
    this->lower_bounds = lower_bounds;
    this->empty_ = false;
  }
  template <typename D>
  void setUpperBounds(const Eigen::MatrixBase<D> & upper_bounds) {
    this->upper_bounds = upper_bounds;
    this->empty_ = false;
  }

  /* ------------------------------------------------------------------------ */
  [[nodiscard]] inline auto size(const SpatialDirection & direction) const
      -> Real {
    return upper_bounds(direction) - lower_bounds(direction);
  }

  [[nodiscard]] auto size() const -> Vector<Real> {
    Vector<Real> size(dim);
    for (auto s : arange(this->dim)) {
      size(s) = this->size(SpatialDirection(s));
    }
    return size;
  }

  inline operator bool() const { return not empty_; }

  /* ------------------------------------------------------------------------ */
  [[nodiscard]] auto allSum(const Communicator & communicator) const -> BBox {
    Matrix<Real> reduce_bounds(dim, 2);

    reduce_bounds(0) = lower_bounds;
    reduce_bounds(1) = Real(-1.) * upper_bounds;

    communicator.allReduce(reduce_bounds, SynchronizerOperation::_min);

    BBox global(dim);
    global.lower_bounds = reduce_bounds(0);
    global.upper_bounds = Real(-1.) * reduce_bounds(1);
    global.empty_ = false;
    return global;
  }

  [[nodiscard]] auto allGather(const Communicator & communicator) const
      -> std::vector<BBox> {
    auto prank = communicator.whoAmI();
    auto nb_proc = communicator.getNbProc();
    Array<Real> bboxes_data(nb_proc, dim * (2 + 1));

    auto it_bbox = make_view(bboxes_data, dim, 2 + 1).begin();
    auto && bounds = it_bbox[prank];
    bounds(0) = lower_bounds;
    bounds(1) = upper_bounds;
    bounds(2)(0) = empty_ ? 1. : 0.; // ugly trick

    communicator.allGather(bboxes_data);

    std::vector<BBox> bboxes;
    bboxes.reserve(nb_proc);

    for (auto p : arange(nb_proc)) {
      bboxes.emplace_back(dim);
      auto & bbox = bboxes.back();
      auto && bounds = it_bbox[p];

      bbox.lower_bounds = bounds(0);
      bbox.upper_bounds = bounds(1);
      bbox.empty_ = (bounds(2)(0) == 1.);
    }

    return bboxes;
  }

  [[nodiscard]] auto intersection(const BBox & other,
                                  const Communicator & communicator) const
      -> decltype(auto) {
    // todo: change for a custom reduction algorithm
    auto other_bboxes = other.allGather(communicator);
    std::map<decltype(other_bboxes.size()), BBox> intersections;
    for (auto && [bbox_i, bbox] : enumerate(other_bboxes)) {
      auto && tmp = this->intersection(bbox);
      if (tmp) {
        intersections[bbox_i] = tmp;
      }
    }
    return intersections;
  }

  void printself(std::ostream & stream) const {
    stream << "BBox[";
    if (not empty()) {
      stream << lower_bounds << " - " << upper_bounds;
    }
    stream << "]";
  }

  [[nodiscard]] auto empty() const -> bool { return empty_; }

protected:
  Int dim{0};
  bool empty_{true};
  Vector<Real> lower_bounds;
  Vector<Real> upper_bounds;
};

inline auto operator<<(std::ostream & stream, const BBox & bbox)
    -> std::ostream & {
  bbox.printself(stream);
  return stream;
}

} // namespace akantu

#endif /* AKANTU_AKA_BBOX_HH_ */
