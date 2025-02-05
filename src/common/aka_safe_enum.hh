/**
 * Copyright (©) 2013-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_AKA_SAFE_ENUM_HH_
#define AKANTU_AKA_SAFE_ENUM_HH_

namespace akantu {

/// Safe enumerated type
template <typename def, typename inner = typename def::type>
class safe_enum : public def {
  using type = typename def::type;

public:
  constexpr explicit safe_enum(type v = def::_end_) : val(v) {}

  constexpr auto underlying() const -> inner { return val; }

  constexpr auto operator==(const safe_enum & s) const -> bool {
    return this->val == s.val;
  }
  constexpr auto operator!=(const safe_enum & s) const -> bool {
    return this->val != s.val;
  }
  constexpr auto operator<(const safe_enum & s) const -> bool {
    return this->val < s.val;
  }
  constexpr auto operator<=(const safe_enum & s) const -> bool {
    return this->val <= s.val;
  }
  constexpr auto operator>(const safe_enum & s) const -> bool {
    return this->val > s.val;
  }
  constexpr auto operator>=(const safe_enum & s) const -> bool {
    return this->val >= s.val;
  }

  constexpr operator inner() { return val; };

  // Works only if enumerations are contiguous.
  class const_iterator {
  public:
    constexpr explicit const_iterator(type v) : it(v) {}
    constexpr auto operator++() -> const_iterator & {
      ++it;
      return *this;
    }
    constexpr auto operator*() -> safe_enum {
      return safe_enum(static_cast<type>(it));
    }
    constexpr auto operator!=(const_iterator const & it) -> bool {
      return it.it != this->it;
    }

  private:
    int it;
  };

  constexpr auto begin() const { return const_iterator(def::_begin_); }
  constexpr auto end() const { return const_iterator(def::_end_); }

private:
  inner val;
};

} // namespace akantu

#endif /* AKANTU_AKA_SAFE_ENUM_HH_ */
