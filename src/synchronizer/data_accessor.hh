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
#include "aka_common.hh"
#include "aka_iterators.hh"
#include "communication_buffer.hh"
#include "element.hh"
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_DATA_ACCESSOR_HH_
#define AKANTU_DATA_ACCESSOR_HH_

namespace akantu {
class FEEngine;
} // namespace akantu

namespace akantu {

class DataAccessorBase {
public:
  DataAccessorBase() = default;
  virtual ~DataAccessorBase() = default;

  DataAccessorBase(const DataAccessorBase & other) = default;
  DataAccessorBase(DataAccessorBase && other) = default;

  DataAccessorBase & operator=(const DataAccessorBase & other) = default;
  DataAccessorBase & operator=(DataAccessorBase && other) = default;
};

template <class T> class DataAccessor : public virtual DataAccessorBase {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  DataAccessor() = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /**
   * @brief get  the number of  data to exchange  for a given array of T
   * (elements or dofs) and a given akantu::SynchronizationTag
   */
  virtual Int getNbData(const Array<T> & elements,
                        const SynchronizationTag & tag) const = 0;

  /**
   * @brief pack the data for a given array of T (elements or dofs) and a given
   * akantu::SynchronizationTag
   */
  virtual void packData(CommunicationBuffer & buffer, const Array<T> & element,
                        const SynchronizationTag & tag) const = 0;

  /**
   * @brief unpack the data for a given array of T (elements or dofs) and a
   * given akantu::SynchronizationTag
   */
  virtual void unpackData(CommunicationBuffer & buffer,
                          const Array<T> & element,
                          const SynchronizationTag & tag) = 0;
};

/* -------------------------------------------------------------------------- */
/* Specialization                                                             */
/* -------------------------------------------------------------------------- */
template <> class DataAccessor<Element> : public virtual DataAccessorBase {
public:
  DataAccessor() = default;

  [[nodiscard]] virtual Int
  getNbData(const Array<Element> & /*elements*/,
            const SynchronizationTag & /*tag*/) const {
    return 0;
  };
  virtual void packData(CommunicationBuffer & /*buffer*/,
                        const Array<Element> & /*element*/,
                        const SynchronizationTag & /*tag*/) const {};
  virtual void unpackData(CommunicationBuffer & /*buffer*/,
                          const Array<Element> & /*element*/,
                          const SynchronizationTag & /*tag*/) {};

private:
  template <
      typename T, bool pack_helper, class Func,
      std::enable_if_t<not std::is_base_of_v<FEEngine, std::decay_t<Func>>> * =
          nullptr>
  static void packUnpackElementalDataHelper(
      std::conditional_t<pack_helper, const ElementTypeMapArray<T>,
                         ElementTypeMapArray<T>> & data_to_pack,
      CommunicationBuffer & buffer, const Array<Element> & element,
      Func && data_per_element);

  /* ------------------------------------------------------------------------ */
public:
  template <typename T, bool pack_helper>
  static void packUnpackNodalDataHelper(
      std::conditional_t<pack_helper, const Array<T>, Array<T>> & data,
      CommunicationBuffer & buffer, const Array<Element> & elements,
      const Mesh & mesh);

  /* ------------------------------------------------------------------------ */
  template <typename T, bool pack_helper>
  static void packUnpackElementalDataHelper(
      std::conditional_t<pack_helper, const ElementTypeMapArray<T>,
                         ElementTypeMapArray<T>> & data_to_pack,
      CommunicationBuffer & buffer, const Array<Element> & element);

  /* ------------------------------------------------------------------------ */
  template <typename T, bool pack_helper>
  static void packUnpackElementalDataHelper(
      std::conditional_t<pack_helper, const ElementTypeMapArray<T>,
                         ElementTypeMapArray<T>> & data_to_pack,
      CommunicationBuffer & buffer, const Array<Element> & element,
      const FEEngine & fem);

  /* ------------------------------------------------------------------------ */
  template <typename T>
  static void
  packNodalDataHelper(const Array<T> & data, CommunicationBuffer & buffer,
                      const Array<Element> & elements, const Mesh & mesh) {
    packUnpackNodalDataHelper<T, true>(data, buffer, elements, mesh);
  }

  template <typename T>
  static inline void
  unpackNodalDataHelper(Array<T> & data, CommunicationBuffer & buffer,
                        const Array<Element> & elements, const Mesh & mesh) {
    packUnpackNodalDataHelper<T, false>(data, buffer, elements, mesh);
  }

  /* ------------------------------------------------------------------------ */
  template <typename T>
  static inline void
  packElementalDataHelper(const ElementTypeMapArray<T> & data_to_pack,
                          CommunicationBuffer & buffer,
                          const Array<Element> & elements) {
    packUnpackElementalDataHelper<T, true>(data_to_pack, buffer, elements);
  }

  template <typename T>
  static inline void
  unpackElementalDataHelper(ElementTypeMapArray<T> & data_to_unpack,
                            CommunicationBuffer & buffer,
                            const Array<Element> & elements) {
    packUnpackElementalDataHelper<T, false>(data_to_unpack, buffer, elements);
  }
  /* ------------------------------------------------------------------------ */
  template <typename T>
  static inline void packElementalDataHelper(
      const ElementTypeMapArray<T> & data_to_pack, CommunicationBuffer & buffer,
      const Array<Element> & elements, const FEEngine & fem) {
    packUnpackElementalDataHelper<T, true>(data_to_pack, buffer, elements, fem);
  }

  template <typename T>
  static inline void unpackElementalDataHelper(
      ElementTypeMapArray<T> & data_to_unpack, CommunicationBuffer & buffer,
      const Array<Element> & elements, const FEEngine & fem) {
    packUnpackElementalDataHelper<T, false>(data_to_unpack, buffer, elements,
                                            fem);
  }
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
template <> class DataAccessor<Idx> : public virtual DataAccessorBase {
public:
  DataAccessor() = default;

  [[nodiscard]] virtual Int
  getNbData(const Array<Idx> & /*elements*/,
            const SynchronizationTag & /*tag*/) const {
    return 0;
  }
  virtual void packData(CommunicationBuffer & /*buffer*/,
                        const Array<Idx> & /*element*/,
                        const SynchronizationTag & /*tag*/) const {}
  virtual void unpackData(CommunicationBuffer & /*buffer*/,
                          const Array<Idx> & /*element*/,
                          const SynchronizationTag & /*tag*/) {}
  /* ------------------------------------------------------------------------ */
public:
  template <typename T, bool pack_helper>
  static void packUnpackDOFDataHelper(
      std::conditional_t<pack_helper, const Array<T>, Array<T>> & data,
      CommunicationBuffer & buffer, const Array<Idx> & dofs);

  template <typename T>
  static inline void packDOFDataHelper(const Array<T> & data_to_pack,
                                       CommunicationBuffer & buffer,
                                       const Array<Idx> & dofs) {
    packUnpackDOFDataHelper<T, true>(data_to_pack, buffer, dofs);
  }

  template <typename T>
  static inline void unpackDOFDataHelper(Array<T> & data_to_unpack,
                                         CommunicationBuffer & buffer,
                                         const Array<Idx> & dofs) {
    packUnpackDOFDataHelper<T, false>(data_to_unpack, buffer, dofs);
  }
};

/* -------------------------------------------------------------------------- */
template <typename T> struct AddOperation {
  inline T operator()(const T & a, const T & b) { return a + b; }
};

template <typename T> struct IdentityOperation {
  inline T operator()(const T & /*unused*/, const T & b) { return b; }
};

template <typename T> struct OrOperation {
  inline T operator()(const T & /*a*/, const T & /*b*/) {
    AKANTU_TO_IMPLEMENT();
  }
};

template <> struct AddOperation<Vector<bool>> {
  inline Vector<bool> operator()(const Vector<bool> & /*a*/,
                                 const Vector<bool> & /*b*/) {
    AKANTU_TO_IMPLEMENT();
  }
};

template <> struct OrOperation<Vector<bool>> {
  inline Vector<bool> operator()(const Vector<bool> & a,
                                 const Vector<bool> & b) {
    Vector<bool> c(a.size());
    for (auto && [c_, a_, b_] : zip(c, a, b)) {
      c_ = a_ | b_;
    }

    return c;
  }
};

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
template <class Entity, template <class> class Op, class T>
class ReduceDataAccessor : public virtual DataAccessor<Entity> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  ReduceDataAccessor(Array<T> & data, const SynchronizationTag & tag)
      : data(data), tag(tag) {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  [[nodiscard]] Int getNbData(const Array<Entity> & entities,
                              const SynchronizationTag & tag) const override {
    if (tag != this->tag) {
      return 0;
    }

    Vector<T> tmp(data.getNbComponent());
    return entities.size() * CommunicationBuffer::sizeInBuffer(tmp);
  }

  /* ------------------------------------------------------------------------ */
  void packData(CommunicationBuffer & buffer, const Array<Entity> & entities,
                const SynchronizationTag & tag) const override {
    if (tag != this->tag) {
      return;
    }

    for (auto && vect :
         filter(entities, make_view(data, data.getNbComponent()))) {
      buffer << vect;
    }
  }

  /* ------------------------------------------------------------------------ */
  void unpackData(CommunicationBuffer & buffer, const Array<Entity> & entities,
                  const SynchronizationTag & tag) override {
    if (tag != this->tag) {
      return;
    }

    Vector<T> unpacked(data.getNbComponent());
    auto data_it = data.begin(data.getNbComponent());
    for (auto && vect :
         filter(entities, make_view(data, data.getNbComponent()))) {
      buffer >> unpacked;
      vect = oper(vect, unpacked);
    }
  }

protected:
  /// data to (un)pack
  Array<T> & data;

  /// Tag to consider
  SynchronizationTag tag;

  /// reduction operator
  Op<Vector<T>> oper;
};

/* -------------------------------------------------------------------------- */
template <class T>
using SimpleIdxDataAccessor = ReduceDataAccessor<Idx, IdentityOperation, T>;

/* -------------------------------------------------------------------------- */
template <class T>
class SimpleElementDataAccessor : public virtual DataAccessor<Element> {
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  SimpleElementDataAccessor(ElementTypeMapArray<T> & data,
                            const SynchronizationTag & tag)
      : data(data), tag(tag) {}

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /* ------------------------------------------------------------------------ */
  [[nodiscard]] Int getNbData(const Array<Element> & elements,
                              const SynchronizationTag & tag) const override {
    if (tag != this->tag) {
      return 0;
    }

    Int size = 0;

    for (const auto & el : elements) {
      auto && data_type = data(el.type, el.ghost_type);
      size += sizeof(T) * data_type.getNbComponent();
    }

    return size;
  }

  /* ------------------------------------------------------------------------ */
  void packData(CommunicationBuffer & buffer, const Array<Element> & elements,
                const SynchronizationTag & tag) const override {
    if (tag != this->tag) {
      return;
    }

    for (const auto & el : elements) {
      auto && data_type = data(el.type, el.ghost_type);
      for (auto c : arange(data_type.getNbComponent())) {
        const auto & data_per_element = data_type(el.element, c);
        buffer << data_per_element;
      }
    }
  }

  /* ------------------------------------------------------------------------ */
  void unpackData(CommunicationBuffer & buffer, const Array<Element> & elements,
                  const SynchronizationTag & tag) override {
    if (tag != this->tag) {
      return;
    }

    for (const auto & el : elements) {
      auto && data_type = data(el.type, el.ghost_type);
      for (auto c : arange(data_type.getNbComponent())) {
        auto & data_per_element = data_type(el.element, c);
        buffer >> data_per_element;
      }
    }
  }

protected:
  /// data to (un)pack
  ElementTypeMapArray<T> & data;

  /// Tag to consider
  SynchronizationTag tag;
};

} // namespace akantu

#include "data_accessor_tmpl.hh"

#endif /* AKANTU_DATA_ACCESSOR_HH_ */
