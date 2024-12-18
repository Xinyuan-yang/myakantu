/**
 * @file   aka_zip_iterator.hh
 *
 * @author Nicolas Richart
 *
 * @date creation  jeu déc 12 2019
 *
 * @brief A Documented file.
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
#include "aka_compatibilty_with_cpp_standard.hh"
#include "aka_iterator_tools.hh"
#include "aka_tuple_tools.hh"
/* -------------------------------------------------------------------------- */
#include <iterator>
#include <utility>
/* -------------------------------------------------------------------------- */

#ifndef AKA_ZIP_ITERATOR_HH
#define AKA_ZIP_ITERATOR_HH

#ifndef AKANTU_ITERATORS_NAMESPACE
#define AKANTU_ITERATORS_NAMESPACE akantu
#endif

#define FWD(x) std::forward<decltype(x)>(x)

namespace AKANTU_ITERATORS_NAMESPACE {

/* -------------------------------------------------------------------------- */
namespace iterators AKA_ITERATOR_EXPORT_NAMESPACE {
  /* ------------------------------------------------------------------------ */
  template <template <class...> class Tuple, class... Iterators>
  class ZipIterator_
      : public details::CopyAssignmentEnabler<
            aka::conjunction<std::is_copy_assignable<Iterators>...,
                             std::is_copy_constructible<Iterators>...>::value>,
        public details::MoveAssignmentEnabler<
            aka::conjunction<std::is_move_assignable<Iterators>...,
                             std::is_move_constructible<Iterators>...>::value> {
  private:
    using tuple_t = Tuple<Iterators...>;

  public:
    using value_type =
        Tuple<typename std::iterator_traits<Iterators>::value_type...>;
    using difference_type = std::common_type_t<
        typename std::iterator_traits<Iterators>::difference_type...>;
    using pointer = Tuple<typename std::iterator_traits<Iterators>::pointer...>;
    using reference =
        Tuple<typename std::iterator_traits<Iterators>::reference...>;
    using iterator_category = // std::input_iterator_tag;
        std::common_type_t<
            typename std::iterator_traits<Iterators>::iterator_category...>;

    explicit ZipIterator_(tuple_t iterators)
        : iterators(std::move(iterators)) {}

    /* ---------------------------------------------------------------------- */
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::bidirectional_iterator_tag>::value> * = nullptr>
    auto operator--() -> ZipIterator_ & {
      tuple::foreach ([](auto && it) { --it; }, iterators);
      return *this;
    }

    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::bidirectional_iterator_tag>::value> * = nullptr>
    auto operator--(int) -> ZipIterator_ {
      auto cpy = *this;
      this->operator--();
      return cpy;
    }

    // input iterator ++it
    auto operator++() -> ZipIterator_ & {
      tuple::foreach ([](auto && it) { ++it; }, iterators);
      return *this;
    }

    // input iterator it++
    auto operator++(int) -> ZipIterator_ {
      auto cpy = *this;
      this->operator++();
      return cpy;
    }

    // input iterator it != other_it
    auto operator!=(const ZipIterator_ & other) const -> bool {
      // return tuple::are_not_equal(iterators, other.iterators);
      return std::get<0>(iterators) !=
             std::get<0>(other.iterators); // helps the compiler to optimize
    }

    // input iterator dereference *it
    auto operator*() -> decltype(auto) {
      return tuple::transform([](auto && it) -> decltype(auto) { return *it; },
                              iterators);
    }

    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    auto operator-(const ZipIterator_ & other) const -> difference_type {
      return std::get<0>(this->iterators) - std::get<0>(other.iterators);
    }

    // random iterator it[idx]
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    auto operator[](std::size_t idx) -> decltype(auto) {
      return tuple::transform(
          [idx](auto && it) -> decltype(auto) { return it[idx]; }, iterators);
    }

    // random iterator it + n
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    auto operator+(std::size_t n) const -> decltype(auto) {
      return ZipIterator_(std::forward<tuple_t>(tuple::transform(
          [n](auto && it) -> decltype(auto) { return it + n; }, iterators)));
    }

    // random iterator it - n
    template <class iterator_category_ = iterator_category,
              std::enable_if_t<aka::is_iterator_category_at_least<
                  iterator_category_,
                  std::random_access_iterator_tag>::value> * = nullptr>
    auto operator-(std::size_t n) const -> decltype(auto) {
      return ZipIterator_(std::forward<tuple_t>(tuple::transform(
          [n](auto && it) -> decltype(auto) { return it - n; }, iterators)));
    }

    template <
        class iterator_category_ = iterator_category,
        std::enable_if_t<aka::is_iterator_category_at_least<
            iterator_category_, std::forward_iterator_tag>::value> * = nullptr>
    auto operator==(const ZipIterator_ & other) const -> bool {
      return not tuple::are_not_equal(iterators, other.iterators);
    }

  private:
    tuple_t iterators;
  };

  template <class... Iterators>
  using ZipIterator = ZipIterator_<std::tuple, Iterators...>;

  template <class... Iterators>
  using NamedZipIterator = ZipIterator_<tuple::named_tuple, Iterators...>;
} // namespace AKA_ITERATOR_EXPORT_NAMESPACE

/* -------------------------------------------------------------------------- */
template <class... Iterators>
auto zip_iterator(std::tuple<Iterators...> && iterators_tuple)
    -> decltype(auto) {
  auto zip = iterators::ZipIterator<Iterators...>(FWD(iterators_tuple));
  return zip;
}

template <class... Iterators>
auto zip_iterator(tuple::named_tuple<Iterators...> && iterators_tuple)
    -> decltype(auto) {
  auto zip = iterators::NamedZipIterator<Iterators...>(FWD(iterators_tuple));
  return zip;
}

/* -------------------------------------------------------------------------- */
namespace containers AKA_ITERATOR_EXPORT_NAMESPACE {
  template <template <class...> class Tuple, class... Containers>
  class ZipContainer_ {
  public:
    using containers_t = Tuple<Containers...>;
    using size_type = std::common_type_t<aka::size_type_t<Containers>...>;

    explicit ZipContainer_(Containers &&... containers)
        : containers(std::forward<Containers>(containers)...) {}

    auto begin() const -> decltype(auto) {
      return zip_iterator(
          tuple::transform([](auto && c) { return c.begin(); },
                           std::forward<containers_t>(containers)));
    }

    auto end() const -> decltype(auto) {
      return zip_iterator(
          tuple::transform([](auto && c) { return c.end(); },
                           std::forward<containers_t>(containers)));
    }

    auto begin() -> decltype(auto) {
      return zip_iterator(
          tuple::transform([](auto && c) { return c.begin(); },
                           std::forward<containers_t>(containers)));
    }

    auto end() -> decltype(auto) {
      return zip_iterator(
          tuple::transform([](auto && c) { return c.end(); },
                           std::forward<containers_t>(containers)));
    }

    template <class... OtherContainers>
    static inline auto append(ZipContainer_ && zip_container,
                              OtherContainers &&... cont);

    template <typename Tag, class OtherContainers>
    static inline auto replace(ZipContainer_ && zip_container,
                               OtherContainers && cont);

    template <typename Tag>
    static inline auto remove(ZipContainer_ && zip_container);

  private:
    template <template <class...> class OtherTuple, class... OtherContainers>
    friend class ZipContainer_;

    containers_t containers;
  };

  template <class... Containers>
  using ZipContainer = ZipContainer_<std::tuple, Containers...>;

  template <class... Containers>
  using NamedZipContainer = ZipContainer_<tuple::named_tuple, Containers...>;

} // namespace AKA_ITERATOR_EXPORT_NAMESPACE

/* -------------------------------------------------------------------------- */
template <class... Containers,
          std::enable_if_t<not aka::conjunction<tuple::is_named_tag<
              std::decay_t<Containers>>...>::value> * = nullptr>
constexpr auto zip(Containers &&... conts)
    -> containers::ZipContainer<Containers...> {
  return containers::ZipContainer<Containers...>(
      std::forward<Containers>(conts)...);
}

template <class... Containers,
          std::enable_if_t<aka::conjunction<tuple::is_named_tag<
              std::decay_t<Containers>>...>::value> * = nullptr>
constexpr auto zip(Containers &&... conts)
    -> containers::NamedZipContainer<Containers...> {
  return containers::NamedZipContainer<Containers...>(
      std::forward<Containers>(conts)...);
}

/* -------------------------------------------------------------------------- */
template <class... ZipContainer> inline auto zip_cat(ZipContainer &&... cont) {
  return make_transform_adaptor(
      zip(std::forward<ZipContainer>(cont)...),
      [](auto && value) { return tuple::flatten(value); });
}

template <class ZipContainer, class... OtherContainers,
          std::enable_if_t<std::is_rvalue_reference<ZipContainer &&>::value> * =
              nullptr>
constexpr inline auto zip_append(ZipContainer && zip_container,
                                 OtherContainers &&... cont) {
  return std::decay_t<ZipContainer>::append(FWD(zip_container), FWD(cont)...);
}

template <typename Tag, class ZipContainer, class OtherContainer,
          std::enable_if_t<std::is_rvalue_reference<ZipContainer &&>::value> * =
              nullptr>
constexpr inline auto zip_replace(ZipContainer && zip_container,
                                  OtherContainer && cont) {
  return std::decay_t<ZipContainer>::template replace<Tag>(FWD(zip_container),
                                                           FWD(cont));
}

template <class ZipContainer, class OtherContainer,
          std::enable_if_t<std::is_rvalue_reference_v<ZipContainer &&> and
                           tuple::is_named_tag_v<OtherContainer>> * = nullptr>
constexpr inline auto zip_replace(ZipContainer && zip_container,
                                  OtherContainer && cont) {
  using Tag = decltype(tuple::make_named_tag<
                       typename std::decay_t<OtherContainer>::_tag>());
  return std::decay_t<ZipContainer>::template replace<Tag>(FWD(zip_container),
                                                           FWD(cont._value));
}

template <typename Tag, class ZipContainer,
          std::enable_if_t<std::is_rvalue_reference<ZipContainer &&>::value> * =
              nullptr>
constexpr inline auto zip_remove(ZipContainer && zip_container) {
  return std::decay_t<ZipContainer>::template remove<Tag>(FWD(zip_container));
}

/* -------------------------------------------------------------------------- */
namespace details AKA_ITERATOR_EXPORT_NAMESPACE {
  template <class Tuple, std::size_t... Is,
            std::enable_if_t<not tuple::is_named_tuple<
                std::decay_t<Tuple>>::value> * = nullptr>
  constexpr inline auto
  make_zip_from_tuple_impl(std::index_sequence<Is...> && /*unused*/,
                           Tuple && tuple) {
    return zip(FWD(std::get<Is>(FWD(tuple)))...);
  }

  template <class Tuple, std::size_t... Is,
            std::enable_if_t<
                tuple::is_named_tuple<std::decay_t<Tuple>>::value> * = nullptr>
  constexpr inline auto
  make_zip_from_tuple_impl(std::index_sequence<Is...> && /*unused*/,
                           Tuple && tuple) {
    using namespace tuple;
    return zip(make_named_tag<tuple_name_tag_t<Is, std::decay_t<Tuple>>>() =
                   FWD(std::get<Is>(FWD(tuple)))...);
  }
} // namespace AKA_ITERATOR_EXPORT_NAMESPACE

template <class Tuple>
constexpr inline auto make_zip_from_tuple(Tuple && tuple) {
  return details::make_zip_from_tuple_impl(
      std::make_index_sequence<
          std::tuple_size<typename std::decay_t<Tuple>>::value>{},
      std::forward<Tuple>(tuple));
}

namespace containers AKA_ITERATOR_EXPORT_NAMESPACE {
  template <template <class...> class Tuple, class... Containers>
  template <class... OtherContainers>
  inline auto
  ZipContainer_<Tuple, Containers...>::append(ZipContainer_ && zip_container,
                                              OtherContainers &&... cont) {
    return make_zip_from_tuple(
        tuple::append(FWD(zip_container.containers), FWD(cont)...));
  }

  template <template <class...> class Tuple, class... Containers>
  template <typename Tag, class OtherContainer>
  inline auto
  ZipContainer_<Tuple, Containers...>::replace(ZipContainer_ && zip_container,
                                               OtherContainer && cont) {
    return make_zip_from_tuple(
        tuple::replace<Tag>(FWD(zip_container.containers), FWD(cont)));
  }

  template <template <class...> class Tuple, class... Containers>
  template <typename Tag>
  inline auto
  ZipContainer_<Tuple, Containers...>::remove(ZipContainer_ && zip_container) {
    return make_zip_from_tuple(
        tuple::remove<Tag>(FWD(zip_container.containers)));
  }
} // namespace AKA_ITERATOR_EXPORT_NAMESPACE

} // namespace AKANTU_ITERATORS_NAMESPACE

#undef FWD

namespace std {
template <template <class...> class Tuple, typename... Its>
struct iterator_traits<
    ::AKANTU_ITERATORS_NAMESPACE::iterators::ZipIterator_<Tuple, Its...>> {
private:
  using iterator_type =
      typename ::AKANTU_ITERATORS_NAMESPACE::iterators::ZipIterator_<Tuple,
                                                                     Its...>;

public:
  using iterator_category = typename iterator_type::iterator_category;
  using value_type = typename iterator_type::value_type;
  using difference_type = typename iterator_type::difference_type;
  using pointer = typename iterator_type::pointer;
  using reference = typename iterator_type::reference;
};
} // namespace std

#endif /* AKA_ZIP_ITERATOR_HH */
