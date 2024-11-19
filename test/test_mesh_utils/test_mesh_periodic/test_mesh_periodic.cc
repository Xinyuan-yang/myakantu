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

#include "test_mesh_periodic_fixture.hh"


TEST_F(TestMeshPeriodic2D, X) {
  mesh->makePeriodic(_x);
  checkAll();
  dump("_2d_x");
}

TEST_F(TestMeshPeriodic2D, XY) {
  mesh->makePeriodic(_x);
  mesh->makePeriodic(_y);
  checkAll();
  dump("_2d_xy");
}

TEST_F(TestMeshPeriodic3D, X) {
  mesh->makePeriodic(_x);
  checkAll();
  dump("_3d_x");
}

TEST_F(TestMeshPeriodic3D, XY) {
  mesh->makePeriodic(_x);
  mesh->makePeriodic(_y);
  checkAll();
  dump("_3d_xy");
}

TEST_F(TestMeshPeriodic3D, XYZ) {
  mesh->makePeriodic(_x);
  mesh->makePeriodic(_y);
  mesh->makePeriodic(_z);
  checkAll();
  dump("_3d_xyz");
}
