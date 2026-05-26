/**
 * Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "integration_point.hh"
#include "mesh.hh"
/* -------------------------------------------------------------------------- */
#include "element_type_conversion.hh"
/* -------------------------------------------------------------------------- */
#include <functional>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_ELEMENT_TYPE_MAP_TMPL_HH_
#define AKANTU_ELEMENT_TYPE_MAP_TMPL_HH_

namespace akantu {

/* -------------------------------------------------------------------------- */
/* ElementTypeMap                                                             */
/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline std::string
ElementTypeMap<Stored, SupportType>::printType(SupportType type,
                                               GhostType ghost_type) {
  std::stringstream sstr;
  sstr << "(" << ghost_type << ":" << type << ")";
  return sstr.str();
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline bool
ElementTypeMap<Stored, SupportType>::exists(SupportType type,
                                            GhostType ghost_type) const {
  return this->getData(ghost_type).find(type) !=
         this->getData(ghost_type).end();
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline const Stored &
ElementTypeMap<Stored, SupportType>::operator()(SupportType type,
                                                GhostType ghost_type) const {
  auto it = this->getData(ghost_type).find(type);

  if (it == this->getData(ghost_type).end()) {
    AKANTU_SILENT_EXCEPTION("No element of type "
                            << ElementTypeMap::printType(type, ghost_type)
                            << " in this ElementTypeMap<"
                            << debug::demangle(typeid(Stored).name())
                            << "> class");
  }
  return it->second;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline Stored &
ElementTypeMap<Stored, SupportType>::operator()(SupportType type,
                                                GhostType ghost_type) {
  return this->getData(ghost_type)[type];
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
template <typename U>
inline Stored &
ElementTypeMap<Stored, SupportType>::operator()(U && insertee, SupportType type,
                                                GhostType ghost_type) {
  auto it = this->getData(ghost_type).find(type);

  if (it != this->getData(ghost_type).end()) {
    AKANTU_SILENT_EXCEPTION("Element of type "
                            << ElementTypeMap::printType(type, ghost_type)
                            << " already in this ElementTypeMap<"
                            << debug::demangle(typeid(Stored).name())
                            << "> class");
  } else {
    auto & data = this->getData(ghost_type);
    const auto & res =
        data.insert(std::make_pair(type, std::forward<U>(insertee)));
    it = res.first;
  }

  return it->second;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline typename ElementTypeMap<Stored, SupportType>::DataMap &
ElementTypeMap<Stored, SupportType>::getData(GhostType ghost_type) {
  if (ghost_type == _not_ghost) {
    return data;
  }

  return ghost_data;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline const typename ElementTypeMap<Stored, SupportType>::DataMap &
ElementTypeMap<Stored, SupportType>::getData(GhostType ghost_type) const {
  if (ghost_type == _not_ghost) {
    return data;
  }

  return ghost_data;
}

/* -------------------------------------------------------------------------- */
/// Works only if stored is a pointer to a class with a printself method
template <class Stored, typename SupportType>
void ElementTypeMap<Stored, SupportType>::printself(std::ostream & stream,
                                                    int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "ElementTypeMap<" << debug::demangle(typeid(Stored).name())
         << "> ["
         << "\n";
  for (auto && gt : ghost_types) {
    const DataMap & data = getData(gt);
    for (auto & pair : data) {
      stream << space << space << ElementTypeMap::printType(pair.first, gt)
             << "\n";
    }
  }

  stream << space << "]"
         << "\n";
}

/* -------------------------------------------------------------------------- */
/* ElementTypeMapArray                                                        */
/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
void ElementTypeMapArray<T, SupportType>::copy(
    const ElementTypeMapArray & other) {
  for (auto ghost_type : ghost_types) {
    for (auto type :
         this->elementTypes(_all_dimensions, ghost_type, _ek_not_defined)) {
      const auto & array_to_copy = other(type, ghost_type);
      auto & array =
          this->alloc(0, array_to_copy.getNbComponent(), type, ghost_type);
      array.copy(array_to_copy);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
ElementTypeMapArray<T, SupportType>::ElementTypeMapArray(
    const ElementTypeMapArray & other)
    : parent(), id(other.id + "_copy"), name(other.name + "_copy") {
  this->copy(other);
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
auto ElementTypeMapArray<T, SupportType>::operator=(
    const ElementTypeMapArray & other) -> ElementTypeMapArray & {
  if (this != &other) {
    AKANTU_DEBUG_WARNING("You are copying the ElementTypeMapArray "
                         << this->id << " are you sure it is on purpose");

    this->id = other.id + "_copy";
    this->name = other.name + "_copy";
    this->copy(other);
  }
  return *this;
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline Array<T> & ElementTypeMapArray<T, SupportType>::alloc(
    Int size, Int nb_component, SupportType type, GhostType ghost_type,
    const T & default_value) {
  std::string ghost_id;
  if (ghost_type == _ghost) {
    ghost_id = ":ghost";
  }

  auto it = this->getData(ghost_type).find(type);

  if (it == this->getData(ghost_type).end()) {
    auto id = this->id + ":" + std::to_string(type) + ghost_id;

    this->getData(ghost_type)[type] =
        std::make_unique<Array<T>>(size, nb_component, default_value, id);
    return *(this->getData(ghost_type)[type]);
  }

  AKANTU_DEBUG_INFO("The vector "
                    << this->id << this->printType(type, ghost_type)
                    << " already exists, it is resized instead of allocated.");
  auto && array = *(it->second);
  array.resize(size);
  return array;
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline void ElementTypeMapArray<T, SupportType>::alloc(
    Int size, Int nb_component, SupportType type, const T & default_value) {
  this->alloc(size, nb_component, type, _not_ghost, default_value);
  this->alloc(size, nb_component, type, _ghost, default_value);
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline void ElementTypeMapArray<T, SupportType>::free() {
  AKANTU_DEBUG_IN();

  for (auto gt : ghost_types) {
    auto & data = this->getData(gt);
    data.clear();
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline void ElementTypeMapArray<T, SupportType>::clear() {
  for (auto gt : ghost_types) {
    auto & data = this->getData(gt);
    for (auto & vect : data) {
      vect.second->clear();
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline bool ElementTypeMapArray<T, SupportType>::empty() const {
  bool is_empty = true;
  for (auto gt : ghost_types) {
    auto & data = this->getData(gt);
    for (auto & vect : data) {
      is_empty &= vect.second->empty();
      if (not is_empty) {
        return false;
      }
    }
  }
  return is_empty;
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
template <typename ST>
inline void ElementTypeMapArray<T, SupportType>::set(const ST & value) {
  for (auto gt : ghost_types) {
    auto & data = this->getData(gt);
    for (auto & vect : data) {
      vect.second->set(value);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline const Array<T> &
ElementTypeMapArray<T, SupportType>::operator()(SupportType type,
                                                GhostType ghost_type) const {
  auto it = this->getData(ghost_type).find(type);

  if (it == this->getData(ghost_type).end()) {
    AKANTU_SILENT_EXCEPTION("No element of type "
                            << ElementTypeMapArray::printType(type, ghost_type)
                            << " in this const ElementTypeMapArray<"
                            << debug::demangle(typeid(T).name()) << "> class(\""
                            << this->id << "\")");
  }
  return *(it->second);
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline Array<T> &
ElementTypeMapArray<T, SupportType>::operator()(SupportType type,
                                                GhostType ghost_type) {
  auto it = this->getData(ghost_type).find(type);

  if (it == this->getData(ghost_type).end()) {
    AKANTU_SILENT_EXCEPTION("No element of type "
                            << ElementTypeMapArray::printType(type, ghost_type)
                            << " in this ElementTypeMapArray<"
                            << debug::demangle(typeid(T).name())
                            << "> class (\"" << this->id << "\")");
  }

  return *(it->second);
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline void ElementTypeMapArray<T, SupportType>::setArray(SupportType type,
                                                          GhostType ghost_type,
                                                          Array<T> & vect) {
  auto it = this->getData(ghost_type).find(type);

  if (AKANTU_DEBUG_TEST(dblWarning) && it != this->getData(ghost_type).end() &&
      it->second != &vect) {
    AKANTU_DEBUG_WARNING(
        "The Array "
        << this->printType(type, ghost_type)
        << " is already registred, this call can lead to a memory leak.");
  }

  this->getData(ghost_type)[type] = &vect;
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
inline void ElementTypeMapArray<T, SupportType>::onElementsRemoved(
    const ElementTypeMapArray<Idx> & new_numbering) {
  for (auto gt : ghost_types) {
    for (auto && type :
         new_numbering.elementTypes(_all_dimensions, gt, _ek_not_defined)) {
      auto support_type = convertType<ElementType, SupportType>(type);
      if (this->exists(support_type, gt)) {
        const auto & renumbering = new_numbering(type, gt);
        if (renumbering.empty()) {
          continue;
        }

        auto & vect = this->operator()(support_type, gt);
        auto nb_entry_per_element = vect.size() / renumbering.size();
        auto nb_component = vect.getNbComponent();

        Array<T> tmp(renumbering.size() * nb_entry_per_element, nb_component);
        auto tmp_it =
            make_view(tmp, nb_entry_per_element * nb_component).begin();
        Int new_size = 0;

        for (auto && [new_i, vect] :
             zip(renumbering,
                 make_view(vect, nb_entry_per_element * nb_component))) {
          if (new_i != -1) {
            tmp_it[new_i] = vect;
            ++new_size;
          }
        }
        tmp.resize(new_size);
        vect.copy(tmp);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
void ElementTypeMapArray<T, SupportType>::printself(std::ostream & stream,
                                                    int indent) const {
  std::string space(indent, AKANTU_INDENT);

  stream << space << "ElementTypeMapArray<" << debug::demangle(typeid(T).name())
         << "> ["
         << "\n";
  for (UInt g = _not_ghost; g <= _ghost; ++g) {
    auto gt = (GhostType)g;

    const DataMap & data = this->getData(gt);
    typename DataMap::const_iterator it;
    for (it = data.begin(); it != data.end(); ++it) {
      stream << space << space << ElementTypeMapArray::printType(it->first, gt)
             << " ["
             << "\n";
      it->second->printself(stream, indent + 3);
      stream << space << space << " ]"
             << "\n";
    }
  }
  stream << space << "]"
         << "\n";
}

/* -------------------------------------------------------------------------- */
/* SupportType Iterator                                                       */
/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
ElementTypeMap<Stored, SupportType>::type_iterator::type_iterator(
    DataMapIterator & list_begin, DataMapIterator & list_end, Int dim,
    ElementKind ek)
    : list_begin(list_begin), list_end(list_end), dim(dim), kind(ek) {}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline auto ElementTypeMap<Stored, SupportType>::type_iterator::operator*()
    -> reference {
  return list_begin->first;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline auto
ElementTypeMap<Stored, SupportType>::type_iterator::operator*() const
    -> reference {
  return list_begin->first;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline auto ElementTypeMap<Stored, SupportType>::type_iterator::operator++()
    -> type_iterator & {
  ++list_begin;
  while ((list_begin != list_end) &&
         (((dim != _all_dimensions) &&
           (dim != Mesh::getSpatialDimension(list_begin->first))) ||
          ((kind != _ek_not_defined) &&
           (kind != Mesh::getKind(list_begin->first))))) {
    ++list_begin;
  }
  return *this;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
auto ElementTypeMap<Stored, SupportType>::type_iterator::operator++(int)
    -> type_iterator {
  type_iterator tmp(*this);
  operator++();
  return tmp;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline bool ElementTypeMap<Stored, SupportType>::type_iterator::operator==(
    const type_iterator & other) const {
  return this->list_begin == other.list_begin;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
inline bool ElementTypeMap<Stored, SupportType>::type_iterator::operator!=(
    const type_iterator & other) const {
  return this->list_begin != other.list_begin;
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
auto ElementTypeMap<Stored, SupportType>::ElementTypesIteratorHelper::begin()
    -> type_iterator {
  auto && data = container.get().getData(ghost_type);
  auto b = data.begin();
  auto e = data.end();

  // loop until the first valid type
  while ((b != e) &&
         (((dim != _all_dimensions) &&
           (dim != Mesh::getSpatialDimension(b->first))) ||
          ((kind != _ek_not_defined) && (kind != Mesh::getKind(b->first))))) {
    ++b;
  }

  return type_iterator(b, e, dim, kind);
}

template <class Stored, typename SupportType>
auto ElementTypeMap<Stored, SupportType>::ElementTypesIteratorHelper::end()
    -> type_iterator {
  auto && data = container.get().getData(ghost_type);
  auto e = data.end();
  return type_iterator(e, e, dim, kind);
}

/* -------------------------------------------------------------------------- */
template <class Stored, typename SupportType>
auto ElementTypeMap<Stored, SupportType>::elementTypesImpl(
    Int dim, GhostType ghost_type, ElementKind kind) const
    -> ElementTypesIteratorHelper {
  return ElementTypesIteratorHelper(*this, dim, ghost_type, kind);
}

/* -------------------------------------------------------------------------- */
/// standard output stream operator
template <class Stored, typename SupportType>
inline std::ostream &
operator<<(std::ostream & stream,
           const ElementTypeMap<Stored, SupportType> & _this) {
  _this.printself(stream);
  return stream;
}

/* -------------------------------------------------------------------------- */
class ElementTypeMapArrayInitializer { // NOLINT
protected:
  using CompFunc = std::function<Int(ElementType, GhostType)>;

public:
  ElementTypeMapArrayInitializer(const CompFunc & comp_func,
                                 Int spatial_dimension = _all_dimensions,
                                 GhostType ghost_type = _not_ghost,
                                 ElementKind element_kind = _ek_not_defined)
      : comp_func(comp_func), spatial_dimension(spatial_dimension),
        ghost_type(ghost_type), element_kind(element_kind) {}

  [[nodiscard]] GhostType ghostType() const { return ghost_type; }

  [[nodiscard]] virtual Int nbComponent(ElementType type) const {
    return comp_func(type, ghostType());
  }

  [[nodiscard]] virtual bool isNodal() const { return false; }

protected:
  CompFunc comp_func;
  Int spatial_dimension;
  GhostType ghost_type;
  ElementKind element_kind;
};

/* -------------------------------------------------------------------------- */
class MeshElementTypeMapArrayInitializer // NOLINT
    : public ElementTypeMapArrayInitializer {
  using CompFunc = ElementTypeMapArrayInitializer::CompFunc;

public:
  MeshElementTypeMapArrayInitializer(
      const Mesh & mesh, Int nb_component = 1,
      Int spatial_dimension = _all_dimensions,
      GhostType ghost_type = _not_ghost,
      ElementKind element_kind = _ek_not_defined, bool with_nb_element = false,
      bool with_nb_nodes_per_element = false,
      const ElementTypeMapArray<Idx> * filter = nullptr)
      : MeshElementTypeMapArrayInitializer(
            mesh,
            [nb_component](ElementType, GhostType) -> Int {
              return nb_component;
            },
            spatial_dimension, ghost_type, element_kind, with_nb_element,
            with_nb_nodes_per_element, filter) {}

  MeshElementTypeMapArrayInitializer(
      const Mesh & mesh, const CompFunc & comp_func,
      Int spatial_dimension = _all_dimensions,
      GhostType ghost_type = _not_ghost,
      ElementKind element_kind = _ek_not_defined, bool with_nb_element = false,
      bool with_nb_nodes_per_element = false,
      const ElementTypeMapArray<Idx> * filter = nullptr)
      : ElementTypeMapArrayInitializer(comp_func, spatial_dimension, ghost_type,
                                       element_kind),
        mesh(mesh), with_nb_element(with_nb_element),
        with_nb_nodes_per_element(with_nb_nodes_per_element), filter(filter) {}

  [[nodiscard]] decltype(auto) elementTypes() const {
    if (filter != nullptr) {
      return filter->elementTypes(this->spatial_dimension, this->ghost_type,
                                  this->element_kind);
    }
    return mesh.elementTypes(this->spatial_dimension, this->ghost_type,
                             this->element_kind);
  }

  [[nodiscard]] virtual Idx size(ElementType type) const {
    if (with_nb_element) {
      if (filter != nullptr) {
        return (*filter)(type, this->ghost_type).size();
      }
      return mesh.getNbElement(type, this->ghost_type);
    }

    return 0;
  }

  [[nodiscard]] Int nbComponent(ElementType type) const override {
    auto res = ElementTypeMapArrayInitializer::nbComponent(type);
    if (with_nb_nodes_per_element) {
      return (res * Mesh::getNbNodesPerElement(type));
    }

    return res;
  }

  [[nodiscard]] bool isNodal() const override {
    return with_nb_nodes_per_element;
  }

protected:
  const Mesh & mesh;
  bool with_nb_element{false};
  bool with_nb_nodes_per_element{false};
  const ElementTypeMapArray<Idx> * filter{nullptr};
};

/* -------------------------------------------------------------------------- */
class FEEngineElementTypeMapArrayInitializer // NOLINT
    : public MeshElementTypeMapArrayInitializer {
public:
  FEEngineElementTypeMapArrayInitializer(
      const FEEngine & fe_engine, Int nb_component = 1,
      Int spatial_dimension = _all_dimensions,
      GhostType ghost_type = _not_ghost,
      ElementKind element_kind = _ek_not_defined,
      const ElementTypeMapArray<Idx> * filter = nullptr);

  FEEngineElementTypeMapArrayInitializer(
      const FEEngine & fe_engine,
      const ElementTypeMapArrayInitializer::CompFunc & nb_component,
      Int spatial_dimension = _all_dimensions,
      GhostType ghost_type = _not_ghost,
      ElementKind element_kind = _ek_not_defined,
      const ElementTypeMapArray<Idx> * filter = nullptr);

  [[nodiscard]] Int size(ElementType type) const override;

  using ElementTypesIteratorHelper =
      ElementTypeMapArray<Real, ElementType>::ElementTypesIteratorHelper;

  [[nodiscard]] auto elementTypes() const -> std::vector<ElementType>;

protected:
  const FEEngine & fe_engine;
};

/* -------------------------------------------------------------------------- */
template <typename T, typename SupportType>
template <class Func>
void ElementTypeMapArray<T, SupportType>::initialize(const Func & f,
                                                     const T & default_value,
                                                     bool do_not_default) {
  this->is_nodal = f.isNodal();
  auto ghost_type = f.ghostType();
  for (const auto & type : f.elementTypes()) {
    if (not this->exists(type, ghost_type)) {
      if (do_not_default) {
        auto & array = this->alloc(0, f.nbComponent(type), type, ghost_type);
        array.resize(f.size(type));
      } else {
        this->alloc(f.size(type), f.nbComponent(type), type, ghost_type,
                    default_value);
      }
    } else {
      auto & array = this->operator()(type, ghost_type);
      if (not do_not_default) {
        array.resize(f.size(type), default_value);
      } else {
        array.resize(f.size(type));
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
/**
 * All parameters are named optionals
 *  \param _nb_component a functor giving the number of components per
 *  (ElementType, GhostType) pair or a scalar giving a unique number of
 * components
 *  regardless of type
 *  \param _spatial_dimension a filter for the elements of a specific dimension
 *  \param _element_kind filter with element kind (_ek_regular, _ek_structural,
 * ...)
 *  \param _with_nb_element allocate the arrays with the number of elements for
 * each
 *  type in the mesh
 *  \param _with_nb_nodes_per_element multiply the number of components by the
 *  number of nodes per element
 *  \param _default_value default inital value
 *  \param _do_not_default do not initialize the allocated arrays
 *  \param _ghost_type filter a type of ghost
 */
template <typename T, typename SupportType>
template <typename... pack>
void ElementTypeMapArray<T, SupportType>::initialize(const Mesh & mesh,
                                                     pack &&... _pack) {
  GhostType requested_ghost_type = OPTIONAL_NAMED_ARG(ghost_type, _casper);
  bool all_ghost_types =
      OPTIONAL_NAMED_ARG(all_ghost_types, requested_ghost_type == _casper);

  for (GhostType ghost_type : ghost_types) {
    if ((not(ghost_type == requested_ghost_type)) and (not all_ghost_types)) {
      continue;
    }

    auto functor = MeshElementTypeMapArrayInitializer(
        mesh, OPTIONAL_NAMED_ARG(nb_component, 1),
        OPTIONAL_NAMED_ARG(spatial_dimension, mesh.getSpatialDimension()),
        ghost_type, OPTIONAL_NAMED_ARG(element_kind, _ek_not_defined),
        OPTIONAL_NAMED_ARG(with_nb_element, false),
        OPTIONAL_NAMED_ARG(with_nb_nodes_per_element, false),
        OPTIONAL_NAMED_ARG(element_filter, nullptr));

    this->initialize(functor, OPTIONAL_NAMED_ARG(default_value, T()),
                     OPTIONAL_NAMED_ARG(do_not_default, false));
  }
}

/* -------------------------------------------------------------------------- */
/**
 * All parameters are named optionals
 *  \param _nb_component a functor giving the number of components per
 *  (ElementType, GhostType) pair or a scalar giving a unique number of
 * components
 *  regardless of type
 *  \param _spatial_dimension a filter for the elements of a specific dimension
 *  \param _element_kind filter with element kind (_ek_regular, _ek_structural,
 * ...)
 *  \param _default_value default inital value
 *  \param _do_not_default do not initialize the allocated arrays
 *  \param _ghost_type filter a specific ghost type
 *  \param _all_ghost_types get all ghost types
 */
template <typename T, typename SupportType>
template <typename... pack>
void ElementTypeMapArray<T, SupportType>::initialize(const FEEngine & fe_engine,
                                                     pack &&... _pack) {
  GhostType requested_ghost_type = OPTIONAL_NAMED_ARG(ghost_type, _casper);
  bool all_ghost_types =
      OPTIONAL_NAMED_ARG(all_ghost_types, requested_ghost_type == _casper);

  for (auto ghost_type : ghost_types) {
    if ((not(ghost_type == requested_ghost_type)) and (not all_ghost_types)) {
      continue;
    }

    auto functor = FEEngineElementTypeMapArrayInitializer(
        fe_engine, OPTIONAL_NAMED_ARG(nb_component, 1),
        OPTIONAL_NAMED_ARG(spatial_dimension, Int(-2)), ghost_type,
        OPTIONAL_NAMED_ARG(element_kind, _ek_not_defined),
        OPTIONAL_NAMED_ARG(element_filter, nullptr));

    this->initialize(functor, OPTIONAL_NAMED_ARG(default_value, T()),
                     OPTIONAL_NAMED_ARG(do_not_default, false));
  }
}

/* -------------------------------------------------------------------------- */
template <class T, typename SupportType>
inline T &
ElementTypeMapArray<T, SupportType>::operator()(const Element & element,
                                                Int component) {
  return this->operator()(element.type, element.ghost_type)(element.element,
                                                            component);
}

/* -------------------------------------------------------------------------- */
template <class T, typename SupportType>
inline const T &
ElementTypeMapArray<T, SupportType>::operator()(const Element & element,
                                                Int component) const {
  return this->operator()(element.type, element.ghost_type)(element.element,
                                                            component);
}

/* -------------------------------------------------------------------------- */
template <class T, typename SupportType>
inline T &
ElementTypeMapArray<T, SupportType>::operator()(const IntegrationPoint & point,
                                                Int component) {
  return this->operator()(point.type, point.ghost_type)(point.global_num,
                                                        component);
}

/* -------------------------------------------------------------------------- */
template <class T, typename SupportType>
inline const T &
ElementTypeMapArray<T, SupportType>::operator()(const IntegrationPoint & point,
                                                Int component) const {
  return this->operator()(point.type, point.ghost_type)(point.global_num,
                                                        component);
}

/* -------------------------------------------------------------------------- */
template <class T, typename SupportType>
inline decltype(auto)
ElementTypeMapArray<T, SupportType>::get(const Element & element) {
  auto & array = operator()(element.type, element.ghost_type);
  auto it = array.begin(array.getNbComponent());
  return it[element.element];
}

/* -------------------------------------------------------------------------- */
template <class T, typename SupportType>
inline decltype(auto)
ElementTypeMapArray<T, SupportType>::get(const Element & element) const {
  const auto & array = operator()(element.type, element.ghost_type);
  auto it = array.begin(array.getNbComponent());
  return it[element.element];
}

/* -------------------------------------------------------------------------- */
template <class T, typename SupportType>
inline decltype(auto)
ElementTypeMapArray<T, SupportType>::get(const IntegrationPoint & point) {
  auto & array = operator()(point.type, point.ghost_type);
  auto it = array.begin(array.getNbComponent());
  return it[point.global_num];
}

/* -------------------------------------------------------------------------- */
template <class T, typename SupportType>
inline decltype(auto)
ElementTypeMapArray<T, SupportType>::get(const IntegrationPoint & point) const {
  const auto & array = operator()(point.type, point.ghost_type);
  auto it = array.begin(array.getNbComponent());
  return it[point.global_num];
}

/* -------------------------------------------------------------------------- */
template <class T, typename SupportType>
template <typename... Ns,
          std::enable_if_t<
              std::conjunction_v<std::is_integral<std::decay_t<Ns>>...> and
              sizeof...(Ns) >= 1> *>
inline decltype(auto)
ElementTypeMapArray<T, SupportType>::get(const Element & element, Ns &&... ns) {
  auto & array = this->operator()(element.type, element.ghost_type);
  auto it = make_view(array, std::forward<Ns>(ns)...).begin();
  return it[element.element];
}

/* -------------------------------------------------------------------------- */
template <class T, typename SupportType>
template <typename... Ns,
          std::enable_if_t<
              std::conjunction_v<std::is_integral<std::decay_t<Ns>>...> and
              sizeof...(Ns) >= 1> *>
inline decltype(auto)
ElementTypeMapArray<T, SupportType>::get(const Element & element,
                                         Ns &&... ns) const {
  const auto & array = this->operator()(element.type, element.ghost_type);
  auto it = make_view(array, std::forward<Ns>(ns)...).begin();
  return it[element.element];
}

/* -------------------------------------------------------------------------- */
template <class T, typename SupportType>
Int ElementTypeMapArray<T, SupportType>::sizeImpl(Int spatial_dimension,
                                                  GhostType ghost_type,
                                                  ElementKind kind) const {
  Int size = 0;
  for (auto && type : this->elementTypes(spatial_dimension, ghost_type, kind)) {
    size += this->operator()(type, ghost_type).size();
  }
  return size;
}

/* -------------------------------------------------------------------------- */
template <class T, typename SupportType>
template <typename... pack>
Int ElementTypeMapArray<T, SupportType>::size(pack &&... _pack) const {
  Int size = 0;
  GhostType requested_ghost_type = OPTIONAL_NAMED_ARG(ghost_type, _casper);
  bool all_ghost_types =
      OPTIONAL_NAMED_ARG(all_ghost_types, requested_ghost_type == _casper);

  for (auto ghost_type : ghost_types) {
    if ((not(ghost_type == requested_ghost_type)) and (not all_ghost_types)) {
      continue;
    }

    size +=
        sizeImpl(OPTIONAL_NAMED_ARG(spatial_dimension, _all_dimensions),
                 ghost_type, OPTIONAL_NAMED_ARG(element_kind, _ek_not_defined));
  }
  return size;
}

} // namespace akantu

#endif /* AKANTU_ELEMENT_TYPE_MAP_TMPL_HH_ */
