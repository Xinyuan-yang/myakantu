/**
 * Copyright (©) 2010-2024 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "dof_manager_default.hh"
#include "mesh.hh"
#include "mesh_accessor.hh"
#include "test_mesh_periodic_fixture.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <ranges>
#include <set>
#include <tuple>
/* -------------------------------------------------------------------------- */

using TestDIMs =
    ::testing::Types<std::tuple<dim_1_t, dim_1_t>, std::tuple<dim_2_t, dim_1_t>,
                     std::tuple<dim_2_t, dim_2_t>, std::tuple<dim_3_t, dim_1_t>,
                     std::tuple<dim_3_t, dim_2_t>,
                     std::tuple<dim_3_t, dim_3_t>>;

TYPED_TEST_SUITE(TestMeshPeriodic, TestDIMs);

TYPED_TEST(TestMeshPeriodic, SimpleMakePeriodic) {
  this->checkAll();

  this->dump("simple");
}

TYPED_TEST(TestMeshPeriodic, MakeConsistentForPeriodicityIdentity) {
  Array<Real> test_val_array(this->mesh->getNbNodes(), 1);
  test_val_array.set(10);

  for (auto && [n, val] : enumerate(test_val_array)) {
    if (this->mesh->isPeriodicSlave(n)) {
      val = -1.;
    }
  }

  this->mesh->makeConsistentForPeriodicity(test_val_array,
                                           ReduceOperation::_id);

  MeshAccessor accessor(*this->mesh);
  const auto & node_flags = accessor.getNodesFlags();

  for (auto && [n, val, flag] : enumerate(test_val_array, node_flags)) {
    EXPECT_EQ(10., val) << "Error node " << n << " with flag " << std::hex
                        << flag;
  }
}

TYPED_TEST(TestMeshPeriodic, MakeConsistentForPeriodicitySum) {
  Array<Real> test_val_array(this->mesh->getNbNodes(), 1);
  test_val_array.set(100.);

  auto && periodic_master_slaves = this->mesh->getPeriodicMasterSlaves();

  auto && keys = std::views::keys(periodic_master_slaves);
  std::set<Idx> masters(keys.begin(), keys.end());

  for (auto && master : masters) {
    test_val_array(master) = 100. / (1. + periodic_master_slaves.count(master));
  }

  for (auto && [slave, master] : this->mesh->getPeriodicSlaveMaster()) {
    test_val_array(slave) = test_val_array(master);
  }

  this->mesh->makeConsistentForPeriodicity(test_val_array,
                                           ReduceOperation::_sum);

  MeshAccessor accessor(*this->mesh);
  const auto & node_flags = accessor.getNodesFlags();

  for (auto && [n, val, flag] : enumerate(test_val_array, node_flags)) {
    EXPECT_NEAR(100., val, 1e-7)
        << "Error node " << n << " with flag " << std::hex << flag;
  }
}
