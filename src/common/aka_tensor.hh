/**
 * Copyright (©) 2022-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef AKA_TENSOR_HH
#define AKA_TENSOR_HH

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <cstddef>
#include <iostream>
#include <ostream>
#include <type_traits>
#include <utility>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
/* TensorBase                                                                 */
/* -------------------------------------------------------------------------- */
template <typename T, Int ndim> class TensorBase : public TensorTrait<ndim> {
  using RetType = TensorBase<T, ndim>;
  static_assert(ndim > 2, "TensorBase cannot by used for dimensions < 3");

protected:
  using idx_t = Idx;

  template <typename... Args>
  using valid_args_t = typename std::enable_if<
      aka::conjunction<aka::disjunction<std::is_integral<Args>,
                                        std::is_enum<Args>>...>::value and
          ndim == sizeof...(Args),
      int>::type;

public:
  using proxy = TensorBase<T, ndim>;
  using size_type = Idx;
  template <Int _ndim = ndim,
            std::enable_if_t<_ndim == 1 or _ndim == 2, int> = 0>
  TensorBase() {
    n.fill(0);
  }

  TensorBase() { n.fill(0); }

  template <typename... Args, valid_args_t<Args...> = 0>
  constexpr TensorBase(Args... args)
      : n{idx_t(args)...}, _size(detail::product_all(args...)) {}

  constexpr TensorBase(const TensorBase & other)
      : n(other.n), _size(other._size), values(other.values) {}

  constexpr TensorBase(TensorBase && other) noexcept
      : n(std::move(other.n)), _size(std::exchange(other._size, 0)),
        values(std::exchange(other.values, nullptr)) {}

protected:
  template <typename Array, idx_t... I>
  constexpr auto check_indices(
      const Array & idx,
      std::integer_sequence<idx_t, I...> /* for_template_deduction */) const {
    bool result = true;
    (void)std::initializer_list<int>{(result &= idx[I] < n[I], 0)...};
    return result;
  }

  template <typename... Args> constexpr auto compute_index(Args... args) const {
    std::array<idx_t, sizeof...(Args)> idx{idx_t(args)...};
    assert(check_indices(
               idx, std::make_integer_sequence<idx_t, sizeof...(Args)>{}) &&
           "The indexes are out of bound");
    idx_t index = 0, i = (sizeof...(Args)) - 1;
    for (; i > 0; i--) {
      index += idx[i];
      if (i > 0) {
        index *= n[i - 1];
      }
    }
    return index + idx[0];
  }

  template <typename S, int... I>
  constexpr auto get_slice(idx_t s, ::std::index_sequence<I...>) {
    return S(this->values + s * detail::product_all(n[I]...), n[I]...);
  }

  template <typename S, ::std::size_t... I>
  constexpr auto get_slice(idx_t s, ::std::index_sequence<I...>) const {
    return S(this->values + s * detail::product_all(n[I]...), n[I]...);
  }

public:
  template <typename... Args, valid_args_t<Args...> = 0>
  inline auto operator()(Args... args) -> T & {
    return *(this->values + compute_index(std::move(args)...));
  }

  template <typename... Args, valid_args_t<Args...> = 0>
  inline auto operator()(Args... args) const -> const T & {
    return *(this->values + compute_index(std::move(args)...));
  }

  template <class R = T, idx_t _ndim = ndim,
            std::enable_if_t<(_ndim > 3) and not ::std::is_const<R>::value> * =
                nullptr>
  inline auto operator()(idx_t s) {
    return get_slice<TensorProxy<T, ndim - 1>>(
        s, std::make_index_sequence<ndim - 1>());
  }

  template <idx_t _ndim = ndim, ::std::enable_if_t<(_ndim > 3)> * = nullptr>
  inline auto operator()(idx_t s) const {
    return get_slice<TensorProxy<T, ndim - 1>>(
        s, std::make_index_sequence<ndim - 1>());
  }

  template <class R = T, idx_t _ndim = ndim,
            std::enable_if_t<(_ndim == 3) and not std::is_const<R>::value> * =
                nullptr>
  inline auto operator()(idx_t s) {
    return get_slice<
        Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>>>(
        s, std::make_index_sequence<ndim - 1>());
  }

  template <idx_t _ndim = ndim, std::enable_if_t<_ndim == 3> * = nullptr>
  inline auto operator()(idx_t s) const {
    return get_slice<Eigen::Map<const Eigen::Matrix<
        std::remove_const_t<T>, Eigen::Dynamic, Eigen::Dynamic>>>(
        s, std::make_index_sequence<ndim - 1>());
  }

protected:
  template <class Operator> auto transform(Operator && op) -> RetType & {
    std::transform(this->values, this->values + this->_size, this->values,
                   std::forward<Operator>(op));
    return *(static_cast<RetType *>(this));
  }

  template <class Other, class Operator>
  auto transform(Other && other, Operator && op) -> RetType & {
    AKANTU_DEBUG_ASSERT(_size == other.size(),
                        "The two tensors do not have the same size "
                            << this->_size << " != " << other._size);

    std::transform(this->values, this->values + this->_size, other.values,
                   this->values, std::forward<Operator>(op));
    return *(static_cast<RetType *>(this));
  }

  template <class Operator> auto accumulate(T init, Operator && op) -> T {
    return std::accumulate(this->values, this->values + this->_size,
                           std::move(init), std::forward<Operator>(op));
  }

  template <class Other, class Init, class Accumulate, class Operator>
  auto transform_reduce(Other && other, T init, Accumulate && acc,
                        Operator && op) -> T {
    return std::inner_product(
        this->values, this->values + this->_size, other.data(), std::move(init),
        std::forward<Accumulate>(acc), std::forward<Operator>(op));
  }

  // element wise arithmetic operators -----------------------------------------
public:
  inline decltype(auto) operator+=(const TensorBase & other) {
    return transform(other, [](auto && a, auto && b) { return a + b; });
  }

  /* ------------------------------------------------------------------------ */
  inline auto operator-=(const TensorBase & other) -> TensorBase & {
    return transform(other, [](auto && a, auto && b) { return a - b; });
  }

  /* ------------------------------------------------------------------------ */
  inline auto operator+=(const T & x) -> TensorBase & {
    return transform([&x](auto && a) { return a + x; });
  }

  /* ------------------------------------------------------------------------ */
  inline auto operator-=(const T & x) -> TensorBase & {
    return transform([&x](auto && a) { return a - x; });
  }

  /* ------------------------------------------------------------------------ */
  inline auto operator*=(const T & x) -> TensorBase & {
    return transform([&x](auto && a) { return a * x; });
  }

  /* ---------------------------------------------------------------------- */
  inline auto operator/=(const T & x) -> TensorBase & {
    return transform([&x](auto && a) { return a / x; });
  }

  /// Y = \alpha X + Y
  inline auto aXplusY(const TensorBase & other, const T alpha = 1.)
      -> TensorBase & {
    return transform(other,
                     [&alpha](auto && a, auto && b) { return alpha * a + b; });
  }

  /* ------------------------------------------------------------------------ */
  auto data() -> T * { return values; }
  auto data() const -> const T * { return values; }

  // clang-format off
  [[deprecated("use data instead to be stl compatible")]]
  auto storage() -> T*{
    return values;
  }

  [[deprecated("use data instead to be stl compatible")]]
  auto storage() const -> const T * {
    return values;
  }
  // clang-format on

  auto size() const { return _size; }
  auto size(idx_t i) const {
    AKANTU_DEBUG_ASSERT(i < ndim, "This tensor has only " << ndim
                                                          << " dimensions, not "
                                                          << (i + 1));
    return n[i];
  };

  inline void set(const T & t) { std::fill_n(values, _size, t); };
  inline void zero() { set(T()); };

public:
  /// "Entrywise" norm norm<L_p> @f[ \|\boldsymbol{T}\|_p = \left(
  /// \sum_i^{n[0]}\sum_j^{n[1]}\sum_k^{n[2]} |T_{ijk}|^p \right)^{\frac{1}{p}}
  /// @f]
  template <Int norm_type,
            std::enable_if_t<norm_type == Eigen::Infinity> * = nullptr>
  auto lpNorm() const -> T {
    return accumulate(
        T(), [](auto && init, auto && a) { return init + std::abs(a); });
  }

  template <Int norm_type, std::enable_if_t<norm_type == 1> * = nullptr>
  auto lpNorm() const -> T {
    return accumulate(T(), [](auto && init, auto && a) {
      return std::max(init, std::abs(a));
    });
  }

  template <Int norm_type, std::enable_if_t<norm_type == 2> * = nullptr>
  auto norm() const -> T {
    return std::sqrt(
        accumulate(T(), [](auto && init, auto && a) { return init + a * a; }));
  }

  template <Int norm_type, std::enable_if_t<(norm_type > 2)> * = nullptr>
  auto norm() const -> T {
    return std::pow(accumulate(T(),
                               [](auto && init, auto && a) {
                                 return init + std::pow(a, norm_type);
                               }),
                    1. / norm_type);
  }

  auto norm() const -> T { return lpNorm<2>(); }

protected:
  template <Int N, typename... Args,
            std::enable_if_t<(sizeof...(Args) == ndim), int> = 0>
  void serialize(std::ostream & stream, Args... args) const {
    stream << this->operator()(std::move(args)...);
  }

  template <Int N, typename... Args,
            std::enable_if_t<(sizeof...(Args) < ndim), int> = 0>
  void serialize(std::ostream & stream, Args... args) const {
    stream << "[";
    for (idx_t i = 0; i < n[N]; ++i) {
      if (i != 0) {
        stream << ",";
      }
      serialize<N + 1>(stream, std::move(args)..., i);
    }
    stream << "]";
  }

public:
  void printself(std::ostream & stream) const { serialize<0>(stream); };

protected:
  template <std::size_t... I>
  constexpr decltype(auto) begin(std::index_sequence<I...>) {
    return view_iterator<ViewIteratorHelper_t<sizeof...(I), T>>(values,
                                                                n[I]...);
  }

  template <std::size_t... I>
  constexpr decltype(auto) end(std::index_sequence<I...>) {
    return view_iterator<ViewIteratorHelper_t<sizeof...(I), T>>(values + _size,
                                                                n[I]...);
  }

  template <std::size_t... I>
  constexpr decltype(auto) begin(std::index_sequence<I...>) const {
    return const_view_iterator<ViewIteratorHelper_t<sizeof...(I), T>>(values,
                                                                      n[I]...);
  }

  template <std::size_t... I>
  constexpr decltype(auto) end(std::index_sequence<I...>) const {
    return const_view_iterator<ViewIteratorHelper_t<sizeof...(I), T>>(
        values + _size, n[I]...);
  }

public:
  decltype(auto) begin() { return begin(std::make_index_sequence<ndim - 1>{}); }
  decltype(auto) end() { return end(std::make_index_sequence<ndim - 1>{}); }

  decltype(auto) begin() const {
    return begin(std::make_index_sequence<ndim - 1>{});
  }
  decltype(auto) end() const {
    return end(std::make_index_sequence<ndim - 1>{});
  }

protected:
  // size per dimension
  std::array<idx_t, ndim> n;

  // total storage size
  idx_t _size{0};

  // actual data location
  T * values{nullptr};
};

/* -------------------------------------------------------------------------- */
/* TensorProxy                                                                */
/* -------------------------------------------------------------------------- */
template <typename T, Int ndim> class TensorProxy : public TensorBase<T, ndim> {
private:
  using parent = TensorBase<T, ndim>;

public:
  // proxy constructor
  template <typename... Args>
  constexpr TensorProxy(T * data = reinterpret_cast<T *>(0xdeadbeef),
                        Args... args)
      : parent(args...) {
    this->values = data;
  }

  constexpr TensorProxy(const TensorProxy<T, ndim> & other) : parent(other) {
    this->values = other.values;
  }

  constexpr TensorProxy(const Tensor<T, ndim> & other) : parent(other) {
    this->values = other.values;
  }

  // move constructors ---------------------------------------------------------
  // proxy -> proxy
  TensorProxy(TensorProxy && other) noexcept : parent(other) {}

  auto operator=(const TensorBase<T, ndim> & other) -> TensorProxy & {
    AKANTU_DEBUG_ASSERT(
        other.size() == this->size(),
        "You are trying to copy too a tensors proxy with the wrong size "
            << this->_size << " != " << other._size);

    static_assert(std::is_trivially_copyable<T>{},
                  "Cannot copy a tensor on non trivial types");

    std::copy(other.values, other.values + this->_size, this->values);
    return *this;
  }
};

/* -------------------------------------------------------------------------- */
/* Tensor                                                                     */
/* -------------------------------------------------------------------------- */
template <typename T, Int ndim> class Tensor : public TensorBase<T, ndim> {
private:
  using parent = TensorBase<T, ndim>;

public:
  template <typename... Args> constexpr Tensor(Args... args) : parent(args...) {
    static_assert(
        std::is_trivially_constructible<T>{},
        "Cannot create a tensor on non trivially constructible types");
    this->values = new T[this->_size];
  }

  /* ------------------------------------------------------------------------ */
  virtual ~Tensor() { delete[] this->values; }

  // copy constructors ---------------------------------------------------------
  constexpr Tensor(const Tensor & other) : parent(other) {
    this->values = new T[this->_size];
    std::copy(other.values, other.values + this->_size, this->values);
  }

  constexpr explicit Tensor(const TensorProxy<T, ndim> & other)
      : parent(other) {
    //    static_assert(false, "Copying data are you sure");
    this->values = new T[this->_size];
    std::copy(other.values, other.values + this->_size, this->values);
  }

  // move constructors ---------------------------------------------------------
  // proxy -> proxy, non proxy -> non proxy
  Tensor(Tensor && other) noexcept : parent(other) {}

  // copy operator -------------------------------------------------------------
  /// operator= copy-and-swap
  auto operator=(const TensorBase<T, ndim> & other) -> Tensor & {
    if (&other == this)
      return *this;

    std::cout << "Warning: operator= delete data" << std::endl;
    delete[] this->values;
    this->n = other.n;
    this->_size = other._size;

    static_assert(
        std::is_trivially_constructible<T>{},
        "Cannot create a tensor on non trivially constructible types");

    this->values = new T[this->_size];

    static_assert(std::is_trivially_copyable<T>{},
                  "Cannot copy a tensor on non trivial types");

    std::copy(other.values, other.values + this->_size, this->values);

    return *this;
  }
};

/* -------------------------------------------------------------------------- */
template <typename T> using Tensor3 = Tensor<T, 3>;
template <typename T> using Tensor3Proxy = TensorProxy<T, 3>;
template <typename T> using Tensor3Base = TensorBase<T, 3>;

} // namespace akantu

#endif // AKA_TENSOR_HH
