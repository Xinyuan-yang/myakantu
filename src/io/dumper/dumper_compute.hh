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

#ifndef AKANTU_DUMPER_COMPUTE_HH_
#define AKANTU_DUMPER_COMPUTE_HH_
/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "dumper_field.hh"
#include "dumper_iohelper.hh"
#include "dumper_type_traits.hh"
/* -------------------------------------------------------------------------- */
#include <aka_iterators.hh>
/* -------------------------------------------------------------------------- */
#include <io_helper.hh>
/* -------------------------------------------------------------------------- */
#include <type_traits>
/* -------------------------------------------------------------------------- */

namespace akantu {
namespace dumpers {
  /* ------------------------------------------------------------------------ */
  class ComputeFunctorInterface {
  public:
    virtual ~ComputeFunctorInterface() = default;

    virtual Int getDim() = 0;
    virtual Int getNbComponent(Int old_nb_comp) = 0;
  };

  /* ------------------------------------------------------------------------ */

  template <typename return_type>
  class ComputeFunctorOutput : public ComputeFunctorInterface {
  public:
    ComputeFunctorOutput() = default;
    ~ComputeFunctorOutput() override = default;
  };

  /* ------------------------------------------------------------------------ */
  template <typename input_type, typename return_type>
  class ComputeFunctor : public ComputeFunctorOutput<return_type> {
  public:
    ComputeFunctor() = default;
    ~ComputeFunctor() override = default;

    virtual return_type func(const input_type & /*d*/,
                             Element /*global_index*/) {
      AKANTU_TO_IMPLEMENT();
    }
    virtual return_type func(const input_type & /*d*/) {
      AKANTU_TO_IMPLEMENT();
    }
  };

  /* ------------------------------------------------------------------------ */
  template <class EnumType>
  class ComputeIntFromEnum
      : public ComputeFunctor<Vector<EnumType>, Vector<Int>> {
  public:
    ComputeIntFromEnum() = default;

    inline Vector<Int> func(const Vector<EnumType> & in) override {
      Vector<Int> out(in.size());
      for (auto && data : zip(in, out)) {
        std::get<1>(data) =
            static_cast<std::underlying_type_t<EnumType>>(std::get<0>(data));
      }

      return out;
    }

    Int getDim() override { return 1; };
    Int getNbComponent(Int old_nb_comp) override { return old_nb_comp; };
  };

  /* ------------------------------------------------------------------------ */
  template <typename SubFieldCompute, typename _return_type,
            class support_type_ = typename SubFieldCompute::support_type>
  class FieldCompute : public Field {
    /* ---------------------------------------------------------------------- */
    /* Typedefs                                                               */
    /* ---------------------------------------------------------------------- */
  public:
    using return_type = _return_type;
    using support_type = support_type_;

    using sub_iterator = typename SubFieldCompute::iterator;
    using sub_types = typename SubFieldCompute::types;
    using sub_return_type = typename sub_types::return_type;
    using data_type = typename return_type::value_type;
    using functor_type = ComputeFunctor<sub_return_type, return_type>;

    using types = TypeTraits<data_type, return_type, Array<data_type>>;

  public:
    class iterator {
    public:
      iterator(const sub_iterator & it, functor_type & func)
          : it(it), func(func) {}

      bool operator!=(const iterator & it) const { return it.it != this->it; }

      iterator operator++() {
        ++this->it;
        return *this;
      }

      return_type operator*() { return func.func(*it); }

      /// Do to IOHelper the needs it...
      Int element_type() { return this->it.element_type(); }

    protected:
      sub_iterator it;
      functor_type & func;
    };

    /* ---------------------------------------------------------------------- */
    /* Constructors/Destructors                                               */
    /* ---------------------------------------------------------------------- */
  public:
    FieldCompute(SubFieldCompute & cont,
                 std::unique_ptr<ComputeFunctorInterface> func)
        : sub_field(aka::as_type<SubFieldCompute>(cont.shared_from_this())),
          func(aka::as_type<functor_type>(func.release())) {
      this->checkHomogeneity();
    };

    void registerToDumper(const std::string & id,
                          iohelper::Dumper & dumper) override {
      dumper.addNodeDataField(id, *this);
    }

    /* ---------------------------------------------------------------------- */
    /* Class Members                                                          */
    /* ---------------------------------------------------------------------- */
  public:
    iterator begin() { return iterator(sub_field->begin(), *func); }
    iterator end() { return iterator(sub_field->end(), *func); }

    Int getDim() { return func->getDim(); }

    Int size() {
      throw;
      // return Functor::size();
      return 0;
    }

    void checkHomogeneity() override { this->homogeneous = true; };

    iohelper::DataType getDataType() {
      return iohelper::getDataType<data_type>();
    }

    /// for connection to a FieldCompute
    inline std::shared_ptr<Field> connect(FieldComputeProxy & proxy) override;

    /// for connection to a FieldCompute
    std::unique_ptr<ComputeFunctorInterface>
    connect(HomogenizerProxy & proxy) override;

    /* ---------------------------------------------------------------------- */
    /* Class Members                                                          */
    /* ---------------------------------------------------------------------- */
  public:
    std::shared_ptr<SubFieldCompute> sub_field;
    std::unique_ptr<functor_type> func;
  };

  /* ------------------------------------------------------------------------ */
  template <typename SubFieldCompute, typename _return_type>
  class FieldCompute<SubFieldCompute, _return_type, Element> : public Field {
    /* ---------------------------------------------------------------------- */
    /* Typedefs                                                               */
    /* ---------------------------------------------------------------------- */
  public:
    using return_type = _return_type;
    using support_type = Element;

    using sub_iterator = typename SubFieldCompute::iterator;
    using sub_types = typename SubFieldCompute::types;
    using sub_return_type = typename sub_types::return_type;
    using data_type = typename sub_types::data_type;
    using functor_type = ComputeFunctor<sub_return_type, return_type>;

    using types =
        TypeTraits<data_type, return_type, ElementTypeMapArray<data_type>>;

  public:
    class iterator {
    public:
      iterator(const sub_iterator & it, functor_type & func)
          : it(it), func(func) {}

      bool operator!=(const iterator & it) const { return it.it != this->it; }
      iterator operator++() {
        ++this->it;
        return *this;
      }

      Idx currentGlobalIndex() { return this->it.currentGlobalIndex(); }

      return_type operator*() { return func.func(*it, it.getCurrentElement()); }

      Element getCurrentElement() { return this->it.getCurrentElement(); }

      Int element_type() { return this->it.element_type(); }

    protected:
      sub_iterator it;
      functor_type & func;
    };

    /* ---------------------------------------------------------------------- */
    /* Constructors/Destructors                                               */
    /* ---------------------------------------------------------------------- */
  public:
    FieldCompute(SubFieldCompute & cont,
                 std::unique_ptr<ComputeFunctorInterface> func)
        : sub_field(aka::as_type<SubFieldCompute>(cont.shared_from_this())),
          func(aka::as_type<functor_type>(func.release())) {
      this->checkHomogeneity();
    };

    ~FieldCompute() override = default;

    void registerToDumper(const std::string & id,
                          iohelper::Dumper & dumper) override {
      dumper.addElemDataField(id, *this);
    }

    /* ---------------------------------------------------------------------- */
    /* Class Members                                                          */
    /* ---------------------------------------------------------------------- */
  public:
    iterator begin() { return iterator(sub_field->begin(), *func); }
    iterator end() { return iterator(sub_field->end(), *func); }

    Int getDim() { return func->getDim(); }

    Int size() {
      throw;
      // return Functor::size();
      return 0;
    }

    void checkHomogeneity() override { this->homogeneous = true; };

    template <class T1 = data_type,
              std::enable_if_t<std::is_enum<T1>::value> * = nullptr>
    iohelper::DataType getDataType() {
      return iohelper::getDataType<Int>();
    }

    template <class T1 = data_type,
              std::enable_if_t<not std::is_enum<T1>::value> * = nullptr>
    iohelper::DataType getDataType() {
      return iohelper::getDataType<data_type>();
    }

    /// get the number of components of the hosted field
    ElementTypeMap<Int>
    getNbComponents(Int dim = _all_dimensions,
                    GhostType ghost_type = _not_ghost,
                    ElementKind kind = _ek_not_defined) override {
      ElementTypeMap<Int> nb_components;
      const auto & old_nb_components =
          this->sub_field->getNbComponents(dim, ghost_type, kind);

      for (auto type : old_nb_components.elementTypes(dim, ghost_type, kind)) {
        auto nb_comp = old_nb_components(type, ghost_type);
        nb_components(type, ghost_type) = func->getNbComponent(nb_comp);
      }
      return nb_components;
    };

    /// for connection to a FieldCompute
    inline std::shared_ptr<Field> connect(FieldComputeProxy & proxy) override;

    /// for connection to a FieldCompute
    std::unique_ptr<ComputeFunctorInterface>
    connect(HomogenizerProxy & proxy) override;

    /* ---------------------------------------------------------------------- */
    /* Class Members                                                          */
    /* ---------------------------------------------------------------------- */
  public:
    std::shared_ptr<SubFieldCompute> sub_field;
    std::unique_ptr<functor_type> func;
  };

  /* ------------------------------------------------------------------------ */
  class FieldComputeProxy {
    /* ---------------------------------------------------------------------- */
    /* Constructors/Destructors                                               */
    /* ---------------------------------------------------------------------- */
  public:
    FieldComputeProxy(std::unique_ptr<ComputeFunctorInterface> func)
        : func(std::move(func)) {};

    inline static std::shared_ptr<Field>
    createFieldCompute(std::shared_ptr<Field> & field,
                       std::unique_ptr<ComputeFunctorInterface> func) {
      FieldComputeProxy compute_proxy(std::move(func));
      return field->connect(compute_proxy);
    }

    template <typename T> std::shared_ptr<Field> connectToField(T * ptr) {
      if (aka::is_of_type<ComputeFunctorOutput<Vector<Real>>>(func)) {
        return this->connectToFunctor<Vector<Real>>(ptr);
      }

      if (aka::is_of_type<ComputeFunctorOutput<Vector<Int>>>(func)) {
        return this->connectToFunctor<Vector<Int>>(ptr);
      }

      if (aka::is_of_type<ComputeFunctorOutput<Matrix<Int>>>(func)) {
        return this->connectToFunctor<Matrix<Int>>(ptr);
      }

      if (aka::is_of_type<ComputeFunctorOutput<Matrix<Real>>>(func)) {
        return this->connectToFunctor<Matrix<Real>>(ptr);
      }
      throw;
    }

    template <typename output, typename T>
    std::shared_ptr<Field> connectToFunctor(T * ptr) {
      return std::make_shared<FieldCompute<T, output>>(*ptr, std::move(func));
    }

    template <typename output, typename SubFieldCompute, typename return_type1,
              typename return_type2>
    std::shared_ptr<Field>
    connectToFunctor(FieldCompute<FieldCompute<SubFieldCompute, return_type1>,
                                  return_type2> * /*ptr*/) {
      throw; //    return new FieldCompute<T,output>(*ptr,func);
      return nullptr;
    }

    template <typename output, typename SubFieldCompute, typename return_type1,
              typename return_type2, typename return_type3,
              typename return_type4>
    std::shared_ptr<Field> connectToFunctor(
        FieldCompute<FieldCompute<FieldCompute<FieldCompute<SubFieldCompute,
                                                            return_type1>,
                                               return_type2>,
                                  return_type3>,
                     return_type4> * /*ptr*/) {
      throw; //    return new FieldCompute<T,output>(*ptr,func);
      return nullptr;
    }

    /* ---------------------------------------------------------------------- */
    /* Class Members                                                          */
    /* ---------------------------------------------------------------------- */
  public:
    std::unique_ptr<ComputeFunctorInterface> func;
  };

  /* ------------------------------------------------------------------------ */
  /// for connection to a FieldCompute
  template <typename SubFieldCompute, typename return_type,
            typename support_type_>
  inline std::shared_ptr<Field>
  FieldCompute<SubFieldCompute, return_type, support_type_>::connect(
      FieldComputeProxy & proxy) {
    return proxy.connectToField(this);
  }

  template <typename SubFieldCompute, typename return_type>
  inline std::shared_ptr<Field>
  FieldCompute<SubFieldCompute, return_type, Element>::connect(
      FieldComputeProxy & proxy) {
    return proxy.connectToField(this);
  }

  /* ------------------------------------------------------------------------ */
} // namespace dumpers
} // namespace akantu

#endif /* AKANTU_DUMPER_COMPUTE_HH_ */
