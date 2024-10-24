/**
 * Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef AST_SYNCHRONIZED_ARRAY_HH_
#define AST_SYNCHRONIZED_ARRAY_HH_

/* -------------------------------------------------------------------------- */
// std
#include <vector>

// akantu
#include "aka_array.hh"

namespace akantu {

/* -------------------------------------------------------------------------- */
enum SyncChoice { _added, _deleted };

/* -------------------------------------------------------------------------- */
class SynchronizedArrayBase {
public:
  SynchronizedArrayBase() = default;
  ~SynchronizedArrayBase() = default;

  virtual ID getID() const { return "call should be virtual"; }

  virtual Int syncDeletedElements(std::vector<Idx> & delete_el) = 0;
  virtual Int syncAddedElements(Int nb_added_el) = 0;
};

/* -------------------------------------------------------------------------- */
template <class T>
class SynchronizedArray : public SynchronizedArrayBase, protected Array<T> {
  using parent = Array<T>;
  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  using value_type = typename parent::value_type;
  using size_type = typename parent::size_type;
  using reference = typename parent::reference;
  using pointer_type = typename parent::pointer_type;
  using const_reference = typename parent::const_reference;
  using array_type = parent;

  SynchronizedArray(Int size, Int nb_component, const_reference value,
                    const ID & id, const_reference default_value,
                    const std::string & restart_name);
  ~SynchronizedArray() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
public:
  /// push_back
  template <typename P> inline void push_back(P && value);

  /// erase
  inline void erase(Idx i);
  // template<typename R>
  // inline void erase(const iterator<R> & it);

  /// synchronize elements
  void syncElements(SyncChoice sync_choice);

  /// dump restart file
  void dumpRestartFile(const std::string & file_name) const;

  /// read restart file
  void readRestartFile(const std::string & file_name);

  /// register depending array
  void registerDependingArray(SynchronizedArrayBase & array);

  /// function to print the contain of the class
  void printself(std::ostream & stream, int indent = 0) const override;

  /// find position of element
  using parent::find;

  /// set values to zero
  using parent::zero;

  /// set all entries of the array to the value t
  using parent::set;

  /// set all entries of the array to value t and set default value
  inline void setAndChangeDefault(T t) {
    this->set(t);
    this->default_value = t;
  }

  /// copy the content of an other array
  using parent::copy;

  /// give the address of the memory allocated for this array
  using parent::data;
  using parent::storage;

  // get nb component
  using parent::getNbComponent;

protected:
  Int syncDeletedElements(std::vector<Idx> & del_elements) override;
  Int syncAddedElements(Int nb_add_elements) override;

  /* ------------------------------------------------------------------------ */
  /* Operators                                                                */
  /* ------------------------------------------------------------------------ */
public:
  using parent::operator();

  auto begin() { return parent::begin(); }
  auto begin() const { return parent::begin(); }
  auto end() { return parent::end(); }
  auto end() const { return parent::end(); }

  /* ------------------------------------------------------------------------ */
  /* Accessors                                                                */
  /* ------------------------------------------------------------------------ */
public:
  AKANTU_SET_MACRO(DefaultValue, default_value, T);

  Int size() const { return this->size_; }

  ID getID() const override { return parent::getID(); }

  const Array<T> & getArray() const {
    const Array<T> & a = *(dynamic_cast<const Array<T> *>(this));
    return a;
  }

  Array<T> & getArray() {
    Array<T> & a = *(dynamic_cast<Array<T> *>(this));
    return a;
  }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// init value when new elements added
  T default_value;

  /// restart file_name
  const std::string restart_name;

  /// elements that have been deleted
  std::vector<Idx> deleted_elements;

  /// number of elements to add
  UInt nb_added_elements;

  /// pointers to arrays to be updated
  std::vector<SynchronizedArrayBase *> depending_arrays;
};

/// standard output stream operator
template <typename T>
inline std::ostream & operator<<(std::ostream & stream,
                                 const SynchronizedArray<T> & _this) {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "synchronized_array_inline_impl.hh"

#endif /* AST_SYNCHRONIZED_ARRAY_HH_ */
