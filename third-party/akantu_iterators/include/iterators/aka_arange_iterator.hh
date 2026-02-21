/**
 * @file   aka_arange_iterator.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  jeu déc 12 2019
 *
 * @brief implementation of arange
 *
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * akantu-iterators is free  software: you can redistribute it and/or  modify it
 * under the terms  of the  GNU Lesser  General Public  License as  published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * akantu-iterators is  distributed in the  hope that it  will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public
 * License  for more details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with akantu-iterators. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include "aka_iterator_tools.hh"
/* -------------------------------------------------------------------------- */
#include <cstddef>
#include <iterator>
#include <tuple>
#include <utility>
/* -------------------------------------------------------------------------- */

#ifndef AKA_ARANGE_ITERATOR_HH
#define AKA_ARANGE_ITERATOR_HH

namespace AKANTU_ITERATORS_NAMESPACE {

namespace containers AKA_ITERATOR_EXPORT_NAMESPACE {

  template <class Iterator> class Range {
  public:
    using size_type =
        std::size_t; // typename
                     // std::iterator_traits<Iterator>::difference_type;
    using iterator = Iterator;
    // ugly trick
    using const_iterator = Iterator;

    explicit Range(Iterator && it1, Iterator && it2)
        : iterators(std::forward<Iterator>(it1), std::forward<Iterator>(it2)) {}

    auto begin() const -> decltype(auto) { return iterators.first; }
    auto begin() -> decltype(auto) { return iterators.first; }

    auto end() const -> decltype(auto) { return iterators.second; }
    auto end() -> decltype(auto) { return iterators.second; }

  private:
    std::pair<Iterator, Iterator> iterators;
  };
} // namespace containers AKA_ITERATOR_EXPORT_NAMESPACE

template <class Iterator>
decltype(auto) range(Iterator && it1, Iterator && it2) {
  return containers::Range<Iterator>(std::forward<Iterator>(it1),
                                     std::forward<Iterator>(it2));
}
/* -------------------------------------------------------------------------- */
/* Arange                                                                     */
/* -------------------------------------------------------------------------- */
namespace iterators AKA_ITERATOR_EXPORT_NAMESPACE {
  template <class T> class ArangeIterator {
  public:
    using value_type = T;
    using pointer = T *;
    using reference = T &;
    using difference_type = size_t;
    using iterator_category = std::forward_iterator_tag;

    constexpr ArangeIterator(T value, T step) : value(value), step(step) {}
    constexpr ArangeIterator(const ArangeIterator &) = default;

    constexpr auto operator++() -> ArangeIterator & {
      value += step;
      return *this;
    }

    constexpr auto operator*() const -> T { return value; }

    constexpr auto operator==(const ArangeIterator & other) const -> bool {
      return (value == other.value) and (step == other.step);
    }

    constexpr auto operator!=(const ArangeIterator & other) const -> bool {
      return not operator==(other);
    }

  private:
    T value{0};
    const T step{1};
  };
} // namespace iterators AKA_ITERATOR_EXPORT_NAMESPACE

namespace containers AKA_ITERATOR_EXPORT_NAMESPACE {
  template <class T> class ArangeContainer {
  public:
    using iterator = iterators::ArangeIterator<T>;
    using const_iterator = iterators::ArangeIterator<T>;
    using size_type = T;

    constexpr ArangeContainer(T start, T stop, T step = 1)
        : start(start), stop((stop - start) % step == 0
                                 ? stop
                                 : start + (1 + (stop - start) / step) * step),
          step(step) {}
    explicit constexpr ArangeContainer(T stop) : ArangeContainer(0, stop, 1) {}

    constexpr auto operator[](std::size_t i) -> T {
      T val = start + i * step;
      assert(val < stop && "i is out of range");
      return val;
    }

    constexpr auto size() -> T { return (stop - start) / step; }

    constexpr auto begin() -> iterator { return iterator(start, step); }
    constexpr auto end() -> iterator { return iterator(stop, step); }

  private:
    const T start{0}, stop{0}, step{1};
  };
} // namespace containers AKA_ITERATOR_EXPORT_NAMESPACE

template <class T,
          typename = std::enable_if_t<std::is_integral<std::decay_t<T>>::value>>
inline decltype(auto) arange(const T & stop) {
  return containers::ArangeContainer<T>(stop);
}

template <class T1, class T2,
          typename = std::enable_if_t<
              std::is_integral<std::common_type_t<T1, T2>>::value>>
inline constexpr decltype(auto) arange(const T1 & start, const T2 & stop) {
  return containers::ArangeContainer<std::common_type_t<T1, T2>>(start, stop);
}

template <class T1, class T2, class T3,
          typename = std::enable_if_t<
              std::is_integral<std::common_type_t<T1, T2, T3>>::value>>
inline constexpr decltype(auto) arange(const T1 & start, const T2 & stop,
                                       const T3 & step) {
  return containers::ArangeContainer<std::common_type_t<T1, T2, T3>>(
      start, stop, step);
}
} // namespace AKANTU_ITERATORS_NAMESPACE

namespace std {
template <class T>
struct iterator_traits<
    ::AKANTU_ITERATORS_NAMESPACE::iterators::ArangeIterator<T>> {
private:
  using iterator_type =
      typename ::AKANTU_ITERATORS_NAMESPACE::iterators::ArangeIterator<T>;

public:
  using iterator_category = typename iterator_type::iterator_category;
  using value_type = typename iterator_type::value_type;
  using difference_type = typename iterator_type::difference_type;
  using pointer = typename iterator_type::pointer;
  using reference = typename iterator_type::reference;
};

} // namespace std

#endif /* AKA_ARANGE_ITERATOR_HH */
