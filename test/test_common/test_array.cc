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
#include "test_gtest_utils.hh"
/* -------------------------------------------------------------------------- */
#include <aka_array.hh>
#include <aka_types.hh>
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <memory>
#include <typeindex>
#include <typeinfo>
/* -------------------------------------------------------------------------- */

using namespace akantu;

namespace {

class NonTrivial {
public:
  NonTrivial() = default;
  NonTrivial(int a) : a(a) {};

  bool operator==(const NonTrivial & rhs) { return a == rhs.a; }
  int a{0};
};

bool operator==(const int & a, const NonTrivial & rhs) { return a == rhs.a; }

std::ostream & operator<<(std::ostream & stream, const NonTrivial & _this) {
  stream << _this.a;
  return stream;
}

/* -------------------------------------------------------------------------- */
using TestTypes = ::testing::Types<Real, UInt, NonTrivial>;
/* -------------------------------------------------------------------------- */

::testing::AssertionResult AssertType(const char * /*a_expr*/,
                                      const char * /*b_expr*/,
                                      const std::type_info & a,
                                      const std::type_info & b) {
  if (std::type_index(a) == std::type_index(b))
    return ::testing::AssertionSuccess();

  return ::testing::AssertionFailure()
         << debug::demangle(a.name()) << " != " << debug::demangle(b.name())
         << ") are different";
}

/* -------------------------------------------------------------------------- */

template <typename T> class ArrayConstructor : public ::testing::Test {
protected:
  using type = T;

  void SetUp() override { type_str = debug::demangle(typeid(T).name()); }

  template <typename... P> decltype(auto) construct(P &&... params) {
    return std::make_unique<Array<T>>(std::forward<P>(params)...);
  }

protected:
  std::string type_str;
};

TYPED_TEST_SUITE(ArrayConstructor, TestTypes, );

TYPED_TEST(ArrayConstructor, ConstructDefault1) {
  auto array = this->construct();
  EXPECT_EQ(0, array->size());
  EXPECT_EQ(1, array->getNbComponent());
  EXPECT_STREQ("", array->getID().c_str());
}

TYPED_TEST(ArrayConstructor, ConstructDefault2) {
  auto array = this->construct(1000);
  EXPECT_EQ(1000, array->size());
  EXPECT_EQ(1, array->getNbComponent());
  EXPECT_STREQ("", array->getID().c_str());
}

TYPED_TEST(ArrayConstructor, ConstructDefault3) {
  auto array = this->construct(1000, 10);
  EXPECT_EQ(1000, array->size());
  EXPECT_EQ(10, array->getNbComponent());
  EXPECT_STREQ("", array->getID().c_str());
}

TYPED_TEST(ArrayConstructor, ConstructDefault4) {
  auto array = this->construct(1000, 10, "test");
  EXPECT_EQ(1000, array->size());
  EXPECT_EQ(10, array->getNbComponent());
  EXPECT_STREQ("test", array->getID().c_str());
}

TYPED_TEST(ArrayConstructor, ConstructDefault5) {
  auto array = this->construct(1000, 10, 1);
  EXPECT_EQ(1000, array->size());
  EXPECT_EQ(10, array->getNbComponent());
  EXPECT_EQ(1, array->operator()(10, 6));
  EXPECT_STREQ("", array->getID().c_str());
}

/* -------------------------------------------------------------------------- */
template <typename T> class ArrayFixture : public ArrayConstructor<T> {
public:
  void SetUp() override {
    ArrayConstructor<T>::SetUp();
    array = this->construct(1000, 10);
  }

  void TearDown() override { array.reset(nullptr); }

protected:
  std::unique_ptr<Array<T>> array;
};

TYPED_TEST_SUITE(ArrayFixture, TestTypes, );

TYPED_TEST(ArrayFixture, Copy) {
  Array<typename TestFixture::type> copy(*this->array);

  EXPECT_EQ(1000, copy.size());
  EXPECT_EQ(10, copy.getNbComponent());
  EXPECT_NE(this->array->data(), copy.data());
}

TYPED_TEST(ArrayFixture, Set) {
  auto & arr = *(this->array);
  arr.set(12);
  EXPECT_EQ(12, arr(156, 5));
  EXPECT_EQ(12, arr(520, 7));
  EXPECT_EQ(12, arr(999, 9));
}

TYPED_TEST(ArrayFixture, Resize) {
  auto & arr = *(this->array);

  auto * ptr = arr.data();

  arr.resize(0);
  EXPECT_EQ(0, arr.size());
  EXPECT_TRUE(arr.data() == nullptr or arr.data() == ptr);
  EXPECT_LE(0, arr.getAllocatedSize());

  arr.resize(3000);
  EXPECT_EQ(3000, arr.size());
  EXPECT_LE(3000, arr.getAllocatedSize());

  ptr = arr.data();

  arr.resize(0);
  EXPECT_EQ(0, arr.size());
  EXPECT_TRUE(arr.data() == nullptr or arr.data() == ptr);
  EXPECT_LE(0, arr.getAllocatedSize());
}

TYPED_TEST(ArrayFixture, PushBack) {
  auto & arr = *(this->array);

  auto * ptr = arr.data();

  arr.resize(0);
  EXPECT_EQ(0, arr.size());
  EXPECT_TRUE(arr.data() == nullptr or arr.data() == ptr);
  EXPECT_LE(0, arr.getAllocatedSize());

  arr.resize(3000);
  EXPECT_EQ(3000, arr.size());
  EXPECT_LE(3000, arr.getAllocatedSize());

  ptr = arr.data();

  arr.resize(0);
  EXPECT_EQ(0, arr.size());
  EXPECT_TRUE(arr.data() == nullptr or arr.data() == ptr);
  EXPECT_LE(0, arr.getAllocatedSize());
}

TYPED_TEST(ArrayFixture, ViewVectorDynamic) {
  auto && view = make_view(*this->array, 10);
  EXPECT_NO_THROW(view.begin());
  {
    auto it = view.begin();
    EXPECT_EQ(10, it->size());
    EXPECT_PRED_FORMAT2(AssertType, typeid(*it),
                        typeid(VectorProxy<typename TestFixture::type>));
    EXPECT_PRED_FORMAT2(AssertType, typeid(it[0]),
                        typeid(VectorProxy<typename TestFixture::type>));
  }
}

TYPED_TEST(ArrayFixture, ViewVectorStatic) {
  auto && view = make_view<10>(*this->array);
  EXPECT_NO_THROW(view.begin());
  {
    auto it = view.begin();
    EXPECT_EQ(10, it->size());
    EXPECT_PRED_FORMAT2(AssertType, typeid(*it),
                        typeid(VectorProxy<typename TestFixture::type, 10>));
    EXPECT_PRED_FORMAT2(AssertType, typeid(it[0]),
                        typeid(VectorProxy<typename TestFixture::type, 10>));
  }
}

TYPED_TEST(ArrayFixture, ViewMatrixStatic) {
  auto && view = make_view(*this->array, 2, 5);

  EXPECT_NO_THROW(view.begin());
  {
    auto it = view.begin();
    EXPECT_EQ(10, it->size());
    EXPECT_EQ(2, it->size(0));
    EXPECT_EQ(5, it->size(1));

    EXPECT_PRED_FORMAT2(AssertType, typeid(*it),
                        typeid(MatrixProxy<typename TestFixture::type>));
    EXPECT_PRED_FORMAT2(AssertType, typeid(it[0]),
                        typeid(MatrixProxy<typename TestFixture::type>));
  }
}

TYPED_TEST(ArrayFixture, ViewMatrixDynamic) {
  auto && view = make_view<2, 5>(*this->array);

  EXPECT_NO_THROW(view.begin());
  {
    auto it = view.begin();
    EXPECT_EQ(10, it->size());
    EXPECT_EQ(2, it->size(0));
    EXPECT_EQ(5, it->size(1));
    EXPECT_EQ(2, it->rows());
    EXPECT_EQ(5, it->cols());

    EXPECT_PRED_FORMAT2(AssertType, typeid(*it),
                        typeid(MatrixProxy<typename TestFixture::type, 2, 5>));
    EXPECT_PRED_FORMAT2(AssertType, typeid(it[0]),
                        typeid(MatrixProxy<typename TestFixture::type, 2, 5>));
  }
}

TYPED_TEST(ArrayFixture, ViewVectorWrong) {
  auto && view = make_view(*this->array, 11);
  EXPECT_THROW(view.begin(), debug::ArrayException);
}

TYPED_TEST(ArrayFixture, ViewMatrixWrong) {
  auto && view = make_view(*this->array, 3, 7);
  EXPECT_THROW(view.begin(), debug::ArrayException);
}

TYPED_TEST(ArrayFixture, ViewMatrixIter) {
  std::size_t count = 0;
  for (auto && mat : make_view(*this->array, 10, 10)) {
    EXPECT_EQ(100, mat.size());
    EXPECT_EQ(10, mat.size(0));
    EXPECT_EQ(10, mat.size(1));
    EXPECT_PRED_FORMAT2(AssertType, typeid(mat),
                        typeid(MatrixProxy<typename TestFixture::type>));

    ++count;
  }

  EXPECT_EQ(100, count);
}

TYPED_TEST(ArrayFixture, ConstViewVector) {
  const auto & carray = *this->array;
  auto && view = make_view(carray, 10);
  EXPECT_NO_THROW(view.begin());
  {
    auto it = view.begin();
    EXPECT_EQ(10, it->size());
    EXPECT_PRED_FORMAT2(AssertType, typeid(*it),
                        typeid(VectorProxy<const typename TestFixture::type>));
    EXPECT_PRED_FORMAT2(AssertType, typeid(it[0]),
                        typeid(VectorProxy<const typename TestFixture::type>));
  }
}

TYPED_TEST(ArrayFixture, EnumerateArray) {
  this->array->set(12);
  const auto & carray = *this->array;
  auto && view = enumerate(make_view(carray, 2));
  int i = 0;
  for (auto && data : view) {
    EXPECT_EQ(i, std::get<0>(data));
    EXPECT_EQ(12, std::get<1>(data)[0]);
    EXPECT_EQ(12, std::get<1>(data)[1]);
    ++i;
  }
}

TYPED_TEST(ArrayFixture, ZipArray) {
  this->array->set(12);
  const auto & carray = *this->array;
  auto && view = zip(arange(carray.size() * carray.getNbComponent() / 2),
                     make_view(carray, 2));
  int i = 0;
  for (auto && data : view) {
    EXPECT_EQ(i, std::get<0>(data));
    EXPECT_EQ(12, std::get<1>(data)[0]);
    ++i;
  }
}

TYPED_TEST(ArrayFixture, IteratorIncrement) {
  this->array->set(12);

  auto it = make_view(*this->array, this->array->getNbComponent()).begin() + 10;

  EXPECT_EQ(12, (*it)[0]);
}

TYPED_TEST(ArrayFixture, IteratorFixSize) {
  this->array->set(12);

  auto it = make_view<10>(*this->array).begin();
  EXPECT_EQ(12, (*it)[0]);

  auto && vect = *it;
  vect(0) = 120;
  EXPECT_EQ(120, (*this->array)(0, 0));
}

TYPED_TEST(ArrayFixture, IteratorBracket) {
  this->array->set(12);

  auto && vect =
      make_view(*this->array, this->array->getNbComponent()).begin()[10];

  EXPECT_EQ(12, vect[0]);
}

TYPED_TEST(ArrayFixture, IteratorSimple) {
  this->array->set(12);

  auto it = this->array->begin(this->array->getNbComponent());

  EXPECT_EQ(12, (*it)[0]);
}

TYPED_TEST(ArrayFixture, IteratorThrow) {
  this->array->set(12);

  EXPECT_THROW(this->array->begin(2 * this->array->getNbComponent()),
               debug::Exception);
}

TYPED_TEST(ArrayFixture, IteratorRange) {
  this->array->set(12);

  auto && view = make_view(*this->array, this->array->getNbComponent());

  auto begin = view.begin();
  auto end = view.end();

  for (auto && data : range(begin, end)) {
    EXPECT_EQ(12, data[0]);
  }
}

TYPED_TEST(ArrayFixture, DynamicSizeIteratorFilter) {
  this->array->set(12);
  (*this->array)(3) = 13;
  (*this->array)(50) = 13;

  std::vector<Idx> list_filter{3, 50};
  auto && view = make_view(*this->array, 10);

  for (auto && data : filter(list_filter, view)) {
    EXPECT_EQ(13, data[0]);
    EXPECT_EQ(12, data[1]);
  }
}

TYPED_TEST(ArrayFixture, IteratorFilter) {
  this->array->set(12);
  (*this->array)(3) = 13;
  (*this->array)(50) = 13;

  std::vector<Idx> list_filter{3, 50};
  auto && view = make_view<10>(*this->array);

  for (auto && data : filter(list_filter, view)) {
    EXPECT_EQ(13, data[0]);
    EXPECT_EQ(12, data[1]);
  }
}

} // namespace
