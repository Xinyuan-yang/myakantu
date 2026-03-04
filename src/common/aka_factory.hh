/**
 * Copyright (©) 2017-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
/* -------------------------------------------------------------------------- */
#include <functional>
#include <map>
#include <memory>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_AKA_FACTORY_HH
#define AKANTU_AKA_FACTORY_HH

namespace akantu {

template <class Base, class T = ID, class... Args> class Factory {
  using allocator_t = std::function<std::unique_ptr<Base>(Args...)>;

private:
  Factory() = default;

public:
  Factory(Factory &&) = delete;
  auto operator=(Factory &&) -> Factory & = delete;
  Factory(const Factory &) = delete;
  auto operator=(const Factory &) -> Factory & = delete;

  static auto getInstance() -> Factory & {
    static Factory instance;
    return instance;
  }
  /* ------------------------------------------------------------------------ */
  auto registerAllocator(const T & id, const allocator_t & allocator) -> bool {
    if (allocators.find(id) != allocators.end()) {
      AKANTU_EXCEPTION("The id \"" << id << "\" is already registered in the "
                                   << debug::demangle(typeid(Base).name())
                                   << " factory");
    }
    allocators[id] = allocator;
    return true;
  }

  template <typename... AArgs>
  auto allocate(const T & id, AArgs &&... args) const -> std::unique_ptr<Base> {
    if (allocators.find(id) == allocators.end()) {
      AKANTU_EXCEPTION("The id \"" << id << "\" is not registered in the "
                                   << debug::demangle(typeid(Base).name())
                                   << " factory.");
    }
    return std::forward<std::unique_ptr<Base>>(
        allocators.at(id)(std::forward<AArgs>(args)...));
  }

  auto knows(const T & id) { return allocators.contains(id); }

  auto getPossibleAllocators() -> std::vector<T> {
    std::vector<T> keys;
    for (auto & e : allocators) {
      keys.push_back(e.first);
    }
    return keys;
  }

private:
  std::map<T, allocator_t> allocators;
};

} // namespace akantu

#endif /* AKANTU_AKA_FACTORY_HH */
