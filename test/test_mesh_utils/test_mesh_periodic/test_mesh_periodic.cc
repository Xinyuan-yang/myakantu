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
#include "dof_manager_default.hh"
#include "test_mesh_periodic_fixture.hh"
/* -------------------------------------------------------------------------- */

using TestDIMs =
    ::testing::Types<std::tuple<dim_1_t, dim_1_t>, std::tuple<dim_2_t, dim_1_t>,
                     std::tuple<dim_2_t, dim_2_t>, std::tuple<dim_3_t, dim_1_t>,
                     std::tuple<dim_3_t, dim_2_t>,
                     std::tuple<dim_3_t, dim_3_t>>;

TYPED_TEST_SUITE(TestMeshPeriodic, TestDIMs);

TYPED_TEST(TestMeshPeriodic, SimpleMakePeriodic) {
  std::string periodic_string;
  for (auto sp : spatial_directions<this->dim_periodic>) {
    this->mesh->makePeriodic(sp);
    periodic_string += std::to_string(sp);
  }
  this->checkAll();
  this->dump(std::to_string(this->dim) + "d-" + periodic_string);
}

TYPED_TEST(TestMeshPeriodic, MakeConsistentForPeriodicity) {
  std::string periodic_string;
  for (auto sp : spatial_directions<this->dim_periodic>) {
    this->mesh->makePeriodic(sp);
    periodic_string += std::to_string(sp);
  }

  DOFManagerDefault dof_manager;

  this->dump("assemble-" + std::to_string(this->dim) + "d-" + periodic_string);
}
