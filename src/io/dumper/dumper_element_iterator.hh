/**
 * Copyright (©) 2014-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_DUMPER_ELEMENT_ITERATOR_HH_
#define AKANTU_DUMPER_ELEMENT_ITERATOR_HH_
/* -------------------------------------------------------------------------- */
#include "element.hh"
/* -------------------------------------------------------------------------- */

namespace akantu::dumpers {
/* ------------------------------------------------------------------------ */

template <class types, template <class> class final_iterator>
class element_iterator {
  /* ---------------------------------------------------------------------- */
  /* Typedefs                                                               */
  /* ---------------------------------------------------------------------- */
public:
  using it_type = typename types::it_type;
  using field_type = typename types::field_type;
  using array_type = typename types::array_type;
  using array_iterator = typename types::array_iterator;
  using iterator = final_iterator<types>;

  /* ---------------------------------------------------------------------- */
  /* Constructors/Destructors                                               */
  /* ---------------------------------------------------------------------- */

  element_iterator(const field_type & field,
                   const typename field_type::type_iterator & t_it,
                   const typename field_type::type_iterator & t_it_end,
                   const array_iterator & array_it,
                   const array_iterator & array_it_end,
                   const GhostType ghost_type = _not_ghost)
      : field(field), tit(t_it), tit_end(t_it_end), array_it(array_it),
        array_it_end(array_it_end), ghost_type(ghost_type) {}

  /* ---------------------------------------------------------------------- */
  /* Methods                                                                */
  /* ---------------------------------------------------------------------- */
  auto operator!=(const iterator & it) const -> bool {
    return (ghost_type != it.ghost_type) ||
           (tit != it.tit || (array_it != it.array_it));
  }

  auto operator++() -> iterator & {
    ++array_it;
    while (array_it == array_it_end && tit != tit_end) {
      ++tit;
      if (tit != tit_end) {

        const array_type & vect = field(*tit, ghost_type);
        auto nb_data_per_elem = getNbDataPerElem(*tit);

        auto view = make_view(vect, nb_data_per_elem);
        array_it = view.begin();
        array_it_end = view.end();
      }
    }
    return *(static_cast<iterator *>(this));
  };

  [[nodiscard]] auto getType() const -> ElementType { return *tit; }

  [[nodiscard]] auto element_type() const -> Int {
    return getIOHelperType(*tit);
  }

  [[nodiscard]] auto getCurrentElement() const -> Element {
    return Element{*tit, array_it.getCurrentIndex(), _not_ghost};
  }

  [[nodiscard]] auto getNbDataPerElem(ElementType type) const -> Int {
    if (!nb_data_per_elem.exists(type, ghost_type)) {
      return field(type, ghost_type).getNbComponent();
    }

    return nb_data_per_elem(type, ghost_type);
  }

  void setNbDataPerElem(const ElementTypeMap<Int> & nb_data) {
    this->nb_data_per_elem = nb_data;
  }

  /* ---------------------------------------------------------------------- */
  /* Class Members                                                          */
  /* ---------------------------------------------------------------------- */
protected:
  /// the field to iterate on
  const field_type & field;
  /// field iterator
  typename field_type::type_iterator tit;
  /// field iterator end
  typename field_type::type_iterator tit_end;
  /// array iterator
  array_iterator array_it;
  /// internal iterator end
  array_iterator array_it_end;
  /// ghost type identification
  const GhostType ghost_type;
  /// number of data per element
  ElementTypeMap<Int> nb_data_per_elem;
};

/* ------------------------------------------------------------------------ */
template <typename types>
class elemental_field_iterator
    : public element_iterator<types, elemental_field_iterator> {
public:
  /* ---------------------------------------------------------------------- */
  /* Typedefs                                                               */
  /* ---------------------------------------------------------------------- */

  using parent =
      element_iterator<types, ::akantu::dumpers::elemental_field_iterator>;
  using it_type = typename types::it_type;
  using return_type = typename types::return_type;
  using field_type = typename types::field_type;
  using array_iterator = typename types::array_iterator;

  /* ---------------------------------------------------------------------- */
  /* Constructors/Destructors                                               */
  /* ---------------------------------------------------------------------- */

  elemental_field_iterator(const field_type & field,
                           const typename field_type::type_iterator & t_it,
                           const typename field_type::type_iterator & t_it_end,
                           const array_iterator & array_it,
                           const array_iterator & array_it_end,
                           const GhostType ghost_type = _not_ghost)
      : parent(field, t_it, t_it_end, array_it, array_it_end, ghost_type) {}

  /* ---------------------------------------------------------------------- */
  /* Methods                                                                */
  /* ---------------------------------------------------------------------- */
  return_type operator*() { return *this->array_it; }

private:
};

/* ------------------------------------------------------------------------ */
} // namespace akantu::dumpers

/* -------------------------------------------------------------------------- */

#endif /* AKANTU_DUMPER_ELEMENT_ITERATOR_HH_ */
