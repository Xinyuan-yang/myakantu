/**
 * Copyright (©) 2012-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#ifndef AKANTU_DUMPER_NODAL_FIELD_HH_
#define AKANTU_DUMPER_NODAL_FIELD_HH_

/* -------------------------------------------------------------------------- */
#include "dumper_compute.hh"
#include "dumper_field.hh"
#include <io_helper.hh>
/* -------------------------------------------------------------------------- */

namespace akantu {
namespace dumpers {

  /* ------------------------------------------------------------------------ */
  // This represents a iohelper compatible field
  template <typename T, bool filtered = false, class Container = Array<T>,
            class Filter = Array<Idx>>
  class NodalField : public dumpers::Field {
    /* ---------------------------------------------------------------------- */
    /* Typedefs                                                               */
    /* ---------------------------------------------------------------------- */
  public:
    using support_type = Idx;
    using types = TypeTraits<T, Vector<T>, Container>;

    class iterator
        : public iohelper::iterator<T, iterator, VectorProxy<const T>> {
    public:
      iterator(const T * vect, Int _offset, Int _n, Int _stride,
               const Int * filter)
          : internal_it(vect), offset(_offset), n(_n), stride(_stride),
            filter(filter) {}

      bool operator!=(const iterator & it) const override {
        if (filter != nullptr) {
          return filter != it.filter;
        }
        return internal_it != it.internal_it;
      }

      iterator & operator++() override {
        if (filter != nullptr) {
          ++filter;
        } else {
          internal_it += offset;
        }
        return *this;
      }

      VectorProxy<const T> operator*() override {
        if (filter != nullptr) {
          return VectorProxy<const T>(internal_it + *(filter)*offset + stride,
                                      n);
        }
        return VectorProxy<const T>(internal_it + stride, n);
      }

    private:
      const T * internal_it;
      Int offset, n, stride;
      const Int * filter{nullptr};
    };

    /* ---------------------------------------------------------------------- */
    /* Constructors/Destructors                                               */
    /* ---------------------------------------------------------------------- */
  public:
    NodalField(const Container & _field, Int _n = 0, Int _stride = 0,
               const Filter * filter = nullptr)
        : field(_field), n(_n), stride(_stride), filter(filter), padding(0) {
      AKANTU_DEBUG_ASSERT(((not filtered) and filter == nullptr) or filtered,
                          "Filter passed to unfiltered NodalField!");
      AKANTU_DEBUG_ASSERT((filtered and this->filter != nullptr) or
                              (not filtered),
                          "No filter passed to filtered NodalField!");
      AKANTU_DEBUG_ASSERT(
          (filter != nullptr and this->filter->getNbComponent() == 1) or
              (filter == nullptr),
          "Multi-component filter given to NodalField ("
              << this->filter->getNbComponent()
              << " components detected, sould be 1");
      if (n == 0) {
        this->n = field.getNbComponent() - stride;
      }
    }

    /* ---------------------------------------------------------------------- */
    /* Methods                                                                */
    /* ---------------------------------------------------------------------- */
  public:
    void registerToDumper(const std::string & id,
                          iohelper::Dumper & dumper) override {
      dumper.addNodeDataField(id, *this);
    }

    inline iterator begin() {
      return iterator(field.data(), field.getNbComponent(), n, stride,
                      filter == nullptr ? nullptr : filter->data());
    }

    inline iterator end() {
      return iterator(field.data() + field.getNbComponent() * field.size(),
                      field.getNbComponent(), n, stride,
                      filter == nullptr ? nullptr
                                        : filter->data() + filter->size());
    }

    bool isHomogeneous() override { return true; }
    void checkHomogeneity() override { this->homogeneous = true; }

    virtual Int getDim() {
      if (this->padding) {
        return this->padding;
      }
      return n;
    }

    void setPadding(Int padding) { this->padding = padding; }

    Int size() {
      if (filter != nullptr) {
        return filter->size();
      }
      return field.size();
    }

    inline std::shared_ptr<Field> connect(FieldComputeProxy & proxy) override {
      return proxy.connectToField(this);
    }

    /// for connection to a Homogenizer
    inline std::unique_ptr<ComputeFunctorInterface>
    connect(HomogenizerProxy & /*proxy*/) override {
      throw;
    }

    template <class T1 = T,
              std::enable_if_t<std::is_enum<T1>::value> * = nullptr>
    iohelper::DataType getDataType() {
      return iohelper::getDataType<std::underlying_type_t<T1>>();
    }

    template <class T1 = T,
              std::enable_if_t<not std::is_enum<T1>::value> * = nullptr>
    iohelper::DataType getDataType() {
      return iohelper::getDataType<T>();
    }

    /* ---------------------------------------------------------------------- */
    /* Class Members */
    /* ---------------------------------------------------------------------- */
  private:
    const Container & field;
    Int n, stride;
    const Filter * filter{nullptr};
    Int padding;
  };

} // namespace dumpers
} // namespace akantu
/* -------------------------------------------------------------------------- */
#endif /* AKANTU_DUMPER_NODAL_FIELD_HH_ */
