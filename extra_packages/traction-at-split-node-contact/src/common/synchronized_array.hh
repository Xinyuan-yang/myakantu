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
#include "aka_array.hh"
#include "aka_common.hh"
/* -------------------------------------------------------------------------- */
#include <ostream>
#include <string>
#include <vector>
/* -------------------------------------------------------------------------- */

namespace akantu {

/* -------------------------------------------------------------------------- */
enum SyncChoice { _added, _deleted };

/* -------------------------------------------------------------------------- */
class SynchronizedArrayBase {
public:
  SynchronizedArrayBase() = default;
  SynchronizedArrayBase(const SynchronizedArrayBase &) = default;
  SynchronizedArrayBase(SynchronizedArrayBase &&) noexcept = default;
  auto operator=(const SynchronizedArrayBase &)
      -> SynchronizedArrayBase & = default;
  auto operator=(SynchronizedArrayBase &&) noexcept
      -> SynchronizedArrayBase & = default;

  virtual ~SynchronizedArrayBase() = default;

  [[nodiscard]] virtual auto getID() const -> ID { AKANTU_TO_IMPLEMENT(); }

  virtual auto syncDeletedElements(std::vector<Idx> & delete_el) -> Int = 0;
  virtual auto syncAddedElements(Int nb_added_el) -> Int = 0;
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

  SynchronizedArray(const SynchronizedArray &) = delete;
  SynchronizedArray(SynchronizedArray &&) = delete;
  auto operator=(const SynchronizedArray &) -> SynchronizedArray & = delete;
  auto operator=(SynchronizedArray &&) -> SynchronizedArray & = delete;

  SynchronizedArray(Int size, Int nb_component, const_reference value,
                    const ID & id, const_reference default_value,
                    ID restart_name);

  ~SynchronizedArray() override = default;

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
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
  auto syncDeletedElements(std::vector<Idx> & del_elements) -> Int override;
  auto syncAddedElements(Int nb_add_elements) -> Int override;

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
  AKANTU_SET_MACRO(DefaultValue, default_value, T);

  [[nodiscard]] auto size() const -> Int { return this->size_; }

  [[nodiscard]] auto getID() const -> ID override { return parent::getID(); }

  auto getArray() const -> const Array<T> & {
    return dynamic_cast<const Array<T> &>(*this);
  }

  auto getArray() -> Array<T> & { return dynamic_cast<Array<T> &>(*this); }

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  /// init value when new elements added
  T default_value{0};

  /// restart file_name
  const std::string restart_name;

  /// elements that have been deleted
  std::vector<Idx> deleted_elements;

  /// number of elements to add
  UInt nb_added_elements{0};

  /// pointers to arrays to be updated
  std::vector<SynchronizedArrayBase *> depending_arrays;
};

/// standard output stream operator
template <typename T>
inline auto operator<<(std::ostream & stream,
                       const SynchronizedArray<T> & _this) -> std::ostream & {
  _this.printself(stream);
  return stream;
}

} // namespace akantu

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */
#include "synchronized_array_inline_impl.hh"

#endif /* AST_SYNCHRONIZED_ARRAY_HH_ */
