/**
 * @file   aka_tuple_tools.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 *
 * @date creation: Fri Aug 11 2017
 * @date last modification: Mon Jan 29 2018
 *
 * @brief  iterator interfaces
 *
 *
 * Copyright 2019 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */

/* -------------------------------------------------------------------------- */
#include "aka_compatibilty_with_cpp_standard.hh"
#include "aka_str_hash.hh"
/* -------------------------------------------------------------------------- */
#include <iostream>
#include <tuple>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_NAMED_TUPLE_HH
#define AKANTU_AKA_NAMED_TUPLE_HH

#ifndef AKANTU_ITERATORS_NAMESPACE
#define AKANTU_ITERATORS_NAMESPACE akantu
#endif

namespace AKANTU_ITERATORS_NAMESPACE {

namespace tuple {
  using hash_type = uint64_t;

  /* ------------------------------------------------------------------------ */
  template <typename tag, typename type> struct named_tag {
    using _tag = tag;
    using _type = type;

    template <typename T,
              std::enable_if_t<not std::is_same_v<named_tag, T>> * = nullptr>
    named_tag(T && value) : _value(std::forward<T>(value)) {}

    type _value;
  };

#if (defined(__GNUC__) || defined(__GNUG__))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Weffc++"
#endif
  namespace details {
    template <typename tag> struct named_tag_proxy {
      using _tag = tag;

      template <typename T> auto operator=(T && value) -> decltype(auto) {
        return named_tag<tag, T>{std::forward<T>(value)};
      }
    };
  } // namespace details
#if (defined(__GNUC__) || defined(__GNUG__))
#pragma GCC diagnostic pop
#endif

  template <typename CharT, CharT... chars> struct string_literal {
    static constexpr tuple::hash_type hash =
        hash::details::fnv1a<tuple::hash_type, CharT, chars...>();
    using hash_type = std::integral_constant<tuple::hash_type, hash>;
  };

  /* ------------------------------------------------------------------------ */
  template <typename T> struct is_named_tag : public std::false_type {};
  template <typename tag>
  struct is_named_tag<details::named_tag_proxy<tag>> : public std::true_type {};
  template <typename tag, typename type>
  struct is_named_tag<named_tag<tag, type>> : public std::true_type {};

  template <typename T> constexpr bool is_named_tag_v = is_named_tag<T>::value;
  /* ------------------------------------------------------------------------ */

  template <class... Params>
  struct named_tuple : public std::tuple<typename Params::_type...> {
    using Index = int;
    using Names_t = std::tuple<typename Params::_tag...>;
    using parent = std::tuple<typename Params::_type...>;

    named_tuple(Params &&... params)
        : parent(std::forward<typename Params::_type>(params._value)...) {}

    named_tuple(typename Params::_type &&... args)
        : parent(std::forward<typename Params::_type>(args)...) {}

    template <typename tag, Index Idx>
    struct get_element_index_helper
        : public std::integral_constant<
              Index, (std::is_same_v<std::tuple_element_t<Idx, Names_t>, tag>)
                         ? Idx
                         : get_element_index_helper<tag, Idx + 1>()> {};

    template <typename tag>
    struct get_element_index_helper<tag, sizeof...(Params)>
        : public std::integral_constant<Index, -1> {};

    template <typename NT>
    static constexpr auto get_element_index() noexcept -> Index {
      return get_element_index_helper<typename NT::_tag, 0>::value;
    }

    template <typename NT>
    static constexpr auto get_element_index(NT && /*unused*/) noexcept
        -> Index {
      return get_element_index_helper<typename NT::_tag, 0>::value;
    }

    template <typename NT, std::enable_if_t<is_named_tag_v<NT>> * = nullptr>
    inline constexpr auto get(NT && /*unused*/) noexcept -> decltype(auto) {
      const auto index = get_element_index<NT>();
      static_assert((index != Index(-1)), "wrong named_tag");
      return (std::get<index>(*this));
    }

    template <typename NT, std::enable_if_t<is_named_tag_v<NT>> * = nullptr>
    inline constexpr auto get(NT && /*unused*/) const noexcept
        -> decltype(auto) {
      const auto index = get_element_index<NT>();
      static_assert((index != Index(-1) and is_named_tag_v<NT>),
                    "wrong named_tag");
      return (std::get<index>(*this));
    }

    template <typename NT,
              std::enable_if_t<is_named_tag<NT>::value> * = nullptr>
    inline constexpr auto operator[](NT && name_tag) const noexcept
        -> decltype(auto) {
      return get(std::forward<NT>(name_tag));
    }

    template <typename NT,
              std::enable_if_t<is_named_tag<NT>::value> * = nullptr>
    inline static constexpr auto has(NT && /*unused*/) noexcept -> bool {
      const auto index = get_element_index<NT>();
      return (index != Index(-1));
    }

    template <typename NT,
              std::enable_if_t<is_named_tag<NT>::value> * = nullptr>
    static constexpr auto has() noexcept -> bool {
      constexpr auto index = get_element_index<NT>();
      return (index != Index(-1));
    }

    // template <hash_type HashCode> static constexpr auto has() noexcept ->
    // bool {
    //   constexpr auto index =
    //       get_element_index_helper<std::integral_constant<hash_type,
    //       HashCode>,
    //                                0>::value;
    //   return (index != Index(-1));
    // }
  };

  /* ---------------------------------------------------------------------- */
  template <typename T> struct is_named_tuple : public std::false_type {};
  template <typename... Params>
  struct is_named_tuple<named_tuple<Params...>> : public std::true_type {};

  template <typename T>
  constexpr bool is_named_tuple_v = is_named_tuple<T>::value;
  /* ---------------------------------------------------------------------- */

  template <typename... Params>
  constexpr auto make_named_tuple(Params &&... params) noexcept
      -> named_tuple<Params...> {
    return named_tuple<Params...>(std::forward<Params>(params)...);
  }

  template <typename tag>
  constexpr auto make_named_tag() noexcept -> decltype(auto) {
    return details::named_tag_proxy<tag>{};
  }

  template <typename Tag, class Tuple> constexpr auto has() -> bool {
    return std::decay_t<Tuple>::template has<Tag>();
  }

  template <typename Tag, class Tuple>
  constexpr auto has(Tag && /**/, Tuple && /**/) -> bool {
    return has<std::decay_t<Tag>, std::decay_t<Tuple>>();
  }

  template <typename Tag, class Tuple> constexpr auto has(Tag && /**/) -> bool {
    return has<std::decay_t<Tag>, std::decay_t<Tuple>>();
  }

  template <typename Param, typename Tuple,
            std::enable_if_t<is_named_tag_v<Param>> * = nullptr>
  constexpr auto get(Tuple && tuple) noexcept -> decltype(auto) {
    return tuple.template get<typename Param::hash>();
  }

  template <std::size_t I, class Tuple> struct tuple_name_tag {
    using type = std::decay_t<
        std::tuple_element_t<I, typename std::decay_t<Tuple>::Names_t>>;
  };

  template <std::size_t I, class Tuple>
  using tuple_name_tag_t = typename tuple_name_tag<I, Tuple>::type;

  template <std::size_t I, class Tuple> struct tuple_element {
    using type = std::tuple_element_t<I, typename Tuple::parent>;
  };

  template <std::size_t I, class Tuple>
  using tuple_element_t = typename tuple_element<I, Tuple>::type;

  namespace details {
    template <class Tuple, std::size_t... Is>
    void printTuple(std::ostream & stream, Tuple && /*tuple*/,
                    std::index_sequence<Is...> && /*unused*/) {
      std::initializer_list<int>{
          ((stream << "{"
                   << std::tuple_element_t<
                          Is, typename std::decay_t<Tuple>::Names_t>::value
                   << "}"),
           0)...};
    }
  } // namespace details

  template <class... Params>
  std::ostream & operator<<(std::ostream & stream,
                            const named_tuple<Params...> & tuple) {
    details::printTuple(stream, tuple,
                        std::make_index_sequence<sizeof...(Params)>{});
    return stream;
  }

} // namespace tuple

#if defined(__INTEL_COMPILER)
// intel warnings here
#elif defined(__clang__)
// clang warnings here
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wgnu-string-literal-operator-template"
#elif (defined(__GNUC__) || defined(__GNUG__))
// gcc warnings here
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
#endif

/// this is a GNU exstension
/// http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3599.html
template <typename C, C... chars>
constexpr auto operator""_n() -> decltype(auto) {
  return tuple::make_named_tag<tuple::string_literal<C, chars...>>();
}

#if defined(__clang__)
#pragma clang diagnostic pop
#elif (defined(__GNUC__) || defined(__GNUG__))
#pragma GCC diagnostic pop
#endif

template <typename T> struct named_tuple_type {
  static_assert(tuple::is_named_tuple_v<std::decay_t<T>>,
                "T is not named_tuple type");
  using type = std::decay_t<T>;
};

// just for readability
template <typename T> using named_tuple_t = typename named_tuple_type<T>::type;

} // namespace AKANTU_ITERATORS_NAMESPACE

namespace aka {
template <typename tag, typename type_>
struct size_type<AKANTU_ITERATORS_NAMESPACE::tuple::named_tag<tag, type_>> {
  using type = typename std::decay_t<type_>::size_type;
};
} // namespace aka

/* -------------------------------------------------------------------------- */
#include <iterator>
/* -------------------------------------------------------------------------- */

namespace std {
template <typename tag, typename type>
struct iterator_traits<
    ::AKANTU_ITERATORS_NAMESPACE::tuple::named_tag<tag, type>>
    : public iterator_traits<type> {};

template <class... Types>
struct tuple_size<AKANTU_ITERATORS_NAMESPACE::tuple::named_tuple<Types...>>
    : std::integral_constant<std::size_t, sizeof...(Types)> {};

} // namespace std

#endif /* AKANTU_AKA_NAMED_TUPLE_HH */
