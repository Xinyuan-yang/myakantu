/**
 * @file   constitutive_laws_handler_tmpl.hh
 *
 * @author Nicolas Richart <nicolas.richart@epfl.ch>
 * @author Mohit Pundir <mohit.pundir@epfl.ch>
 *
 * @date creation  ven mar 26 2021
 *
 * @brief Implementation of the consistutive laws handler.
 *
 * @section LICENSE
 *
 * Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
 *
 * Akantu is free  software: you can redistribute it and/or  modify it under the
 * terms  of the  GNU Lesser  General Public  License as  published by  the Free
 * Software Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
 * details.
 *
 * You should  have received  a copy  of the GNU  Lesser General  Public License
 * along with Akantu. If not, see <http://www.gnu.org/licenses/>.
 *
 */
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "constitutive_law_non_local_interface.hh"
#include "constitutive_law_selector_tmpl.hh"
#include "constitutive_laws_handler.hh"
#include "mesh_iterators.hh"
#include "non_local_manager.hh"
#include "parsable.hh"
#include <sstream>
#include <string>
#include <utility>
#include <vector>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_CONSTITUTIVE_LAWS_HANDLER_TMPL_HH
#define AKANTU_CONSTITUTIVE_LAWS_HANDLER_TMPL_HH

namespace akantu {

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
inline decltype(auto)
ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::getConstitutiveLaws() {
  return make_transform_adaptor(
      constitutive_laws, [](auto && value) -> decltype(auto) {
        return aka::as_type<ConstitutiveLawType>(*value);
      });
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
inline decltype(auto)
ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::getConstitutiveLaws()
    const {
  return make_transform_adaptor(
      constitutive_laws, [](auto && value) -> decltype(auto) {
        return aka::as_type<ConstitutiveLawType>(*value);
      });
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
inline Idx
ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::getConstitutiveLawIndex(
    const std::string & name) const {
  auto it = constitutive_laws_names_to_id.find(name);
  if (it == constitutive_laws_names_to_id.end()) {
    AKANTU_SILENT_EXCEPTION(
        "The model " << this->id << " has no constitutive_law named " << name);
  }

  return it->second;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
inline auto &
ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::getConstitutiveLaw(
    Idx cl_index) {
  AKANTU_DEBUG_ASSERT(cl_index < Idx(constitutive_laws.size()),
                      "The model " << this->id
                                   << " has no constitutive law number "
                                   << cl_index);
  return *(constitutive_laws.at(cl_index));
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
inline const auto &
ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::getConstitutiveLaw(
    Idx cl_index) const {
  AKANTU_DEBUG_ASSERT(cl_index < Idx(constitutive_laws.size()),
                      "The model " << this->id
                                   << " has no constitutive law number "
                                   << cl_index);
  return *(constitutive_laws.at(cl_index));
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
inline auto &
ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::getConstitutiveLaw(
    const std::string & name) {
  auto it = constitutive_laws_names_to_id.find(name);
  if (it == constitutive_laws_names_to_id.end()) {
    AKANTU_SILENT_EXCEPTION(
        "The model " << this->id << " has no constitutive_law named " << name);
  }

  return *(constitutive_laws[it->second]);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
inline const auto &
ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::getConstitutiveLaw(
    const std::string & name) const {
  auto it = constitutive_laws_names_to_id.find(name);
  if (it == constitutive_laws_names_to_id.end()) {
    AKANTU_SILENT_EXCEPTION(
        "The model " << this->id << " has no constitutive_law named " << name);
  }
  return *(constitutive_laws[it->second]);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
inline auto &
ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::getConstitutiveLaw(
    const Element & element) {
  auto index = constitutive_law_index(element);
  return *(constitutive_laws[index]);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
inline const auto &
ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::getConstitutiveLaw(
    const Element & element) const {
  auto index = constitutive_law_index(element);
  return *(constitutive_laws[index]);
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
#define AKA_FWD_CL(...) ::std::forward<decltype(__VA_ARGS__)>(__VA_ARGS__)

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
template <typename Operation>
void ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::
    splitByConstitutiveLaw(const Array<Element> & elements,
                           Operation && op) const {
  std::vector<Array<Element>> elements_per_cl(constitutive_laws.size());
  this->splitElementByConstitutiveLaw(elements, elements_per_cl);

  for (auto && cl : zip(constitutive_laws, elements_per_cl)) {
    AKA_FWD_CL(op)(AKA_FWD_CL(*std::get<0>(cl)), AKA_FWD_CL(std::get<1>(cl)));
  }
}

#undef AKA_FWD_CL

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
auto & ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::
    registerNewConstitutiveLaw(const ParserSection & section) {
  std::string cl_name;
  std::string cl_type = section.getName();
  std::string opt_param = section.getOption();

  try {
    std::string tmp = section.getParameter("name");
    cl_name = tmp; /** this can seam weird, but there is an ambiguous operator
                    * overload that i couldn't solve. @todo remove the
                    * weirdness of this code
                    */
  } catch (debug::Exception &) {
    AKANTU_ERROR("A constitutive_law of type \'"
                 << cl_type
                 << "\' in the input file has been defined without a name!");
  }

  auto && cl = this->registerNewConstitutiveLaw(cl_name, cl_type, opt_param);

  cl.parseSection(section);

  return cl;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
auto & ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::
    registerNewConstitutiveLaw(const ID & cl_name, const ID & cl_type,
                               const ID & opt_param) {
  AKANTU_DEBUG_ASSERT(constitutive_laws_names_to_id.find(cl_name) ==
                          constitutive_laws_names_to_id.end(),
                      "A constitutive_law with this name '"
                          << cl_name << "' has already been registered. "
                          << "Please use unique names for constitutive_laws");

  UInt cl_count = constitutive_laws.size();
  constitutive_laws_names_to_id[cl_name] = cl_count;

  std::stringstream sstr_cl;
  sstr_cl << this->id << ":" << cl_count << ":" << cl_type;
  ID cl_id = sstr_cl.str();

  auto constitutive_law = ConstitutiveLawType::getFactory().allocate(
      cl_type, this->spatial_dimension, opt_param,
      aka::as_type<typename ConstitutiveLawType::ConstitutiveLawsHandler>(
          *this),
      cl_id);

  constitutive_laws.push_back(std::move(constitutive_law));

  return aka::as_type<ConstitutiveLawType>(*(constitutive_laws.back()));
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
void ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::
    instantiateConstitutiveLaws(ParserSection & parser_section) {
  auto model_constitutive_laws =
      parser_section.getSubSections(this->parser_type);
  for (const auto & section : model_constitutive_laws) {
    this->registerNewConstitutiveLaw(section);
  }

#ifdef AKANTU_DAMAGE_NON_LOCAL
  for (auto & constitutive_law : constitutive_laws) {
    if (dynamic_cast<ConstitutiveLawNonLocalInterface *>(
            constitutive_law.get()) == nullptr) {
      continue;
    }

    this->non_local_manager = std::make_unique<NonLocalManager>(
        *this, *this, this->id + ":non_local_manager");
    break;
  }
#endif

  if (constitutive_laws.empty()) {
    AKANTU_DEBUG_WARNING(
        "No constitutive_laws where instantiated for the model " << this->id);
  }
  are_constitutive_laws_instantiated = true;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
void ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::
    assignConstitutiveLawToElements(const ElementTypeMapArray<Idx> * filter) {

  for_each_element(
      this->mesh,
      [&](auto && element) {
        Idx cl_index = (*constitutive_law_selector)(element);
        AKANTU_DEBUG_ASSERT(cl_index < Idx(constitutive_laws.size()),
                            "The constitutive_law selector returned an index "
                            "that does not exists");
        constitutive_law_index(element) = cl_index;
      },
      _element_filter = filter, _ghost_type = _not_ghost);

  if (non_local_manager) {
    non_local_manager->synchronize(*this,
                                   SynchronizationTag::_constitutive_law_id);
  }

  for_each_element(
      this->mesh,
      [&](auto && element) {
        auto cl_index = constitutive_law_index(element);
        auto index = constitutive_laws[cl_index]->addElement(element);
        constitutive_law_local_numbering(element) = index;
      },
      _element_filter = filter, _ghost_type = _not_ghost);

  // synchronize the element constitutive_law arrays
  this->synchronize(SynchronizationTag::_constitutive_law_id);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
void ConstitutiveLawsHandler<ConstitutiveLawType,
                             Model_>::initConstitutiveLaws() {
  AKANTU_DEBUG_ASSERT(not constitutive_laws.empty(),
                      "No constitutive_law to initialize !");

  this->assignConstitutiveLawToElements();

  for (auto & constitutive_law : constitutive_laws) {
    /// init internals properties
    constitutive_law->initConstitutiveLaw();
  }

  this->synchronize(SynchronizationTag::_clh_init_cl);

  if (this->non_local_manager) {
    this->non_local_manager->initialize();
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
Int ConstitutiveLawsHandler<
    ConstitutiveLawType, Model_>::getInternalIndexFromID(const ID & id) const {
  AKANTU_DEBUG_IN();

  auto it = constitutive_laws.begin();
  auto end = constitutive_laws.end();

  for (; it != end; ++it) {
    if ((*it)->getID() == id) {
      AKANTU_DEBUG_OUT();
      return (it - constitutive_laws.begin());
    }
  }

  AKANTU_DEBUG_OUT();
  return -1;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
void ConstitutiveLawsHandler<ConstitutiveLawType,
                             Model_>::reassignConstitutiveLaw() {
  AKANTU_DEBUG_IN();

  std::vector<Array<Element>> element_to_add(constitutive_laws.size());
  std::vector<Array<Element>> element_to_remove(constitutive_laws.size());

  for_each_element(this->mesh, [&](auto && element) {
    auto old_constitutive_law = constitutive_law_index(element);
    auto new_constitutive_law = (*constitutive_law_selector)(element);
    if (old_constitutive_law != new_constitutive_law) {
      element_to_add[new_constitutive_law].push_back(element);
      element_to_remove[old_constitutive_law].push_back(element);
    }
  });

  for (auto && data : enumerate(constitutive_laws)) {
    auto cl_index = std::get<0>(data);
    auto & cl = *std::get<1>(data);

    cl.removeElements(element_to_remove[cl_index]);
    cl.addElements(element_to_add[cl_index]);
  }

  AKANTU_DEBUG_OUT();
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
void ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::
    updateDataForNonLocalCriterion(ElementTypeMapReal & criterion) {
  const ID field_name = criterion.getName();
  for (auto & constitutive_law : constitutive_laws) {
    if (not constitutive_law->template isInternal<Real>(field_name,
                                                        _ek_regular)) {
      continue;
    }

    for (auto ghost_type : ghost_types) {
      constitutive_law->flattenInternal(field_name, criterion, ghost_type,
                                        _ek_regular);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
void ConstitutiveLawsHandler<ConstitutiveLawType,
                             Model_>::initializeNonLocal() {
  this->non_local_manager->synchronize(
      *this, SynchronizationTag::_constitutive_law_id);
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
void ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::
    insertIntegrationPointsInNeighborhoods(GhostType ghost_type) {
  for (auto & cl : constitutive_laws) {
    ConstitutiveLawNonLocalInterface * cl_non_local;
    if ((cl_non_local = dynamic_cast<ConstitutiveLawNonLocalInterface *>(
             cl.get())) == nullptr) {
      continue;
    }

    auto & fe_engine = cl->getFEEngine();

    ElementTypeMapArray<Real> quadrature_points_coordinates(
        "quadrature_points_coordinates_tmp_nl", this->id);
    quadrature_points_coordinates.initialize(
        fe_engine, _nb_component = this->spatial_dimension,
        _ghost_type = ghost_type);

    for (const auto & type : quadrature_points_coordinates.elementTypes(
             this->spatial_dimension, ghost_type)) {
      fe_engine.computeIntegrationPointsCoordinates(
          quadrature_points_coordinates(type, ghost_type), type, ghost_type);
    }

    cl_non_local->initConstitutiveLawNonLocal();

    cl_non_local->insertIntegrationPointsInNeighborhoods(
        ghost_type, quadrature_points_coordinates);
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
void ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::updateLocalInternal(
    ElementTypeMapReal & internal_flat, GhostType ghost_type,
    ElementKind kind) {
  const ID field_name = internal_flat.getName();
  for (auto & constitutive_law : constitutive_laws) {
    if (constitutive_law->template isInternal<Real>(field_name, kind)) {
      constitutive_law->flattenInternal(field_name, internal_flat, ghost_type,
                                        kind);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
void ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::
    updateNonLocalInternal(ElementTypeMapReal & internal_flat,
                           GhostType ghost_type, ElementKind kind) {

  const ID field_name = internal_flat.getName();

  for (auto & mat : constitutive_laws) {
    if (not aka::is_of_type<ConstitutiveLawNonLocalInterface>(*mat)) {
      continue;
    }

    auto & mat_non_local =
        dynamic_cast<ConstitutiveLawNonLocalInterface &>(*mat);
    mat_non_local.updateNonLocalInternals(internal_flat, field_name, ghost_type,
                                          kind);
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
void ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::
    splitElementByConstitutiveLaw(
        const Array<Element> & elements,
        std::vector<Array<Element>> & elements_per_cl) const {
  for (const auto & el : elements) {
    Element cl_el = el;
    cl_el.element = this->constitutive_law_local_numbering(el);
    auto cl_id = this->constitutive_law_index(el);
    elements_per_cl[cl_id].push_back(cl_el);
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
Int ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::getNbData(
    const Array<Element> & elements, const SynchronizationTag & tag) const {
  Int size{0};

  if (tag == SynchronizationTag::_constitutive_law_id) {
    size += Int(elements.size()) * sizeof(Idx);
  }

  if (tag != SynchronizationTag::_constitutive_law_id) {
    splitByConstitutiveLaw(elements, [&](auto && cl, auto && elements) {
      size += cl.getNbData(elements, tag);
    });
  }
  return size;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
void ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::packData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & tag) const {
  if (tag == SynchronizationTag::_constitutive_law_id) {
    DataAccessor<Element>::packElementalDataHelper(constitutive_law_index,
                                                   buffer, elements);
  }

  if (tag != SynchronizationTag::_constitutive_law_id) {
    splitByConstitutiveLaw(elements, [&](auto && cl, auto && elements) {
      cl.packData(buffer, elements, tag);
    });
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
void ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::unpackData(
    CommunicationBuffer & buffer, const Array<Element> & elements,
    const SynchronizationTag & tag) {
  if (tag == SynchronizationTag::_constitutive_law_id) {
    for (auto && element : elements) {
      Idx recv_cl_index{};
      buffer >> recv_cl_index;
      Idx & cl_index = constitutive_law_index(element);
      if (cl_index != -1) {
        continue;
      }

      // add ghosts element to the correct constitutive_law
      cl_index = recv_cl_index;
      auto index = constitutive_laws[cl_index]->addElement(element);
      constitutive_law_local_numbering(element) = index;
    }
  }

  if (tag != SynchronizationTag::_constitutive_law_id) {
    splitByConstitutiveLaw(elements, [&](auto && cl, auto && elements) {
      cl.unpackData(buffer, elements, tag);
    });
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
void ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::onElementsAdded(
    const Array<Element> & element_list, const NewElementsEvent & event) {
  this->constitutive_law_index.initialize(
      this->mesh, _element_kind = _ek_not_defined, _with_nb_element = true,
      _default_value = -1);
  this->constitutive_law_local_numbering.initialize(
      this->mesh, _element_kind = _ek_not_defined, _with_nb_element = true,
      _default_value = -1);

  ElementTypeMapArray<Idx> filter("new_element_filter", this->id);

  for (const auto & elem : element_list) {
    if (this->mesh.getSpatialDimension(elem.type) != this->spatial_dimension) {
      continue;
    }

    if (!filter.exists(elem.type, elem.ghost_type)) {
      filter.alloc(0, 1, elem.type, elem.ghost_type);
    }
    filter(elem.type, elem.ghost_type).push_back(elem.element);
  }

  // this fails in parallel if the event is sent on facet between constructor
  // and initFull \todo: to debug...
  this->assignConstitutiveLawToElements(&filter);

  for (auto & constitutive_law : constitutive_laws) {
    constitutive_law->onElementsAdded(element_list, event);
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
void ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::onElementsRemoved(
    const Array<Element> & element_list,
    const ElementTypeMapArray<Idx> & new_numbering,
    const RemovedElementsEvent & event) {
  for (auto & constitutive_law : constitutive_laws) {
    constitutive_law->onElementsRemoved(element_list, new_numbering, event);
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
bool ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::isInternal(
    const std::string & field_name, ElementKind element_kind) {
  /// check if at least one constitutive_law contains field_id as an internal
  for (auto & constitutive_law : constitutive_laws) {
    auto is_internal =
        constitutive_law->template isInternal<Real>(field_name, element_kind);
    if (is_internal) {
      return true;
    }
  }

  return false;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
ElementTypeMap<Int>
ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::getInternalDataPerElem(
    const std::string & field_name, ElementKind element_kind) {

  if (!(this->isInternal(field_name, element_kind))) {
    AKANTU_EXCEPTION("unknown internal " << field_name);
  }

  for (auto & constitutive_law : constitutive_laws) {
    if (constitutive_law->template isInternal<Real>(field_name, element_kind)) {
      return constitutive_law->template getInternalDataPerElem<Real>(
          field_name, element_kind);
    }
  }

  return ElementTypeMap<Int>();
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
template <class T>
ElementTypeMapArray<T> &
ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::flattenInternal(
    const std::string & field_name, ElementKind kind,
    const GhostType ghost_type) {
  auto key = std::make_pair(field_name, kind);

  ElementTypeMapArray<T> * internal_flat;

  auto it = this->registered_internals.find(key);
  if (it == this->registered_internals.end()) {
    auto internal =
        std::make_unique<ElementTypeMapArray<T>>(field_name, this->id);

    internal_flat = internal.get();
    this->registered_internals[key] = std::move(internal);
  } else {
    internal_flat = aka::as_type<ElementTypeMapArray<T>>(it->second.get());
  }

  for (auto type :
       this->mesh.elementTypes(this->spatial_dimension, ghost_type, kind)) {
    if (internal_flat->exists(type, ghost_type)) {
      auto & internal = (*internal_flat)(type, ghost_type);
      internal.resize(0);
    }
  }

  for_each_constitutive_law([&](auto && constitutive_law) {
    if (constitutive_law.template isInternal<T>(field_name, kind)) {
      constitutive_law.flattenInternal(field_name, *internal_flat, ghost_type,
                                       kind);
    }
  });

  return *internal_flat;
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
void ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::
    flattenAllRegisteredInternals(ElementKind kind) {
  for (auto & internal : this->registered_internals) {
    auto && [_id, _kind] = internal.first;
    if (kind == _kind) {
      this->flattenInternal(_id, kind);
    }
  }
}

/* -------------------------------------------------------------------------- */
template <class ConstitutiveLawType, class Model_>
template <class T>
void ConstitutiveLawsHandler<ConstitutiveLawType, Model_>::inflateInternal(
    const std::string & field_name, const ElementTypeMapArray<T> & field,
    GhostType ghost_type, ElementKind kind) {
  for_each_constitutive_law([&](auto && constitutive_law) {
    if (constitutive_law.template isInternal<T>(field_name, kind)) {
      constitutive_law.inflateInternal(field_name, field, ghost_type, kind);
    } else {
      AKANTU_ERROR("A internal of name \'"
                   << field_name
                   << "\' has not been defined in the constitutive law");
    }
  });
}

} // namespace akantu

#endif /* AKANTU_CONSTITUTIVE_LAWS_HANDLER_TMPL_HH */
