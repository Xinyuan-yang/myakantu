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

#include "test_periodic_static_fixture.hh"
#include "dof_synchronizer.hh"

TEST_F(TestPeriodicStatic, None) {
  const double DISPLACEMENT_THRESHOLD = 1e-6;
  const double STRESS_THRESHOLD = 2e-9;

  initAll("none");
  solve();

  Int prank = mesh->getCommunicator().whoAmI();

  if (prank == 0) {
    // Check stress
    Real s_zz = traction(_z);

    for (auto && s : make_view(stress, dim, dim)) {
      EXPECT_NEAR(s(_x, _x) / s_zz, 0.0, STRESS_THRESHOLD);
      EXPECT_NEAR(s(_x, _y) / s_zz, 0.0, STRESS_THRESHOLD);
      EXPECT_NEAR(s(_x, _z) / s_zz, 0.0, STRESS_THRESHOLD);
      EXPECT_NEAR(s(_y, _y) / s_zz, 0.0, STRESS_THRESHOLD);
      EXPECT_NEAR(s(_y, _z) / s_zz, 0.0, STRESS_THRESHOLD);
      EXPECT_NEAR(s(_z, _z) / s_zz, 1.0, STRESS_THRESHOLD);
    }

    // Check displacement
    Real u_z = s_zz / E;
    Real u_xy = -nu * u_z;

    for (auto && [coord, disp] : zip(
      make_view(coordinates, dim),
      make_view(displacements, dim)
    )) {
      // std::cout << "Coordinate: " << coord(_x) << " " << coord(_y) << " " << coord(_z) << std::endl;
      // std::cout << "Displacement: " << disp(_x) << " " << disp(_y) << " " << disp(_z) << std::endl;
      // std::cout << "Expected dis: " << coord(_x) * u_xy << " " << coord(_y) * u_xy << " " << coord(_z) * u_z << std::endl << std::endl;
      EXPECT_NEAR(disp(_x) / u_xy, coord(_x), DISPLACEMENT_THRESHOLD);
      EXPECT_NEAR(disp(_y) / u_xy, coord(_y), DISPLACEMENT_THRESHOLD);
      EXPECT_NEAR(disp(_z) / u_z, coord(_z), DISPLACEMENT_THRESHOLD);
    }
  }
}

TEST_F(TestPeriodicStatic, X) {
  const double DISPLACEMENT_THRESHOLD = 1e-6;
  const double STRESS_THRESHOLD = 2e-9;

  mesh->makePeriodic(_x);
  initAll("x");
  solve();

  Int prank = mesh->getCommunicator().whoAmI();

  if (prank == 0) {
    // Check stress
    Real s_zz = traction(_z);
    Real s_xx = nu * s_zz;

    for (auto && s : make_view(stress, dim, dim)) {
      EXPECT_NEAR(s(_x, _x) / s_xx, 1.0, STRESS_THRESHOLD);
      EXPECT_NEAR(s(_x, _y) / s_zz, 0.0, STRESS_THRESHOLD);
      EXPECT_NEAR(s(_x, _z) / s_zz, 0.0, STRESS_THRESHOLD);
      EXPECT_NEAR(s(_y, _y) / s_zz, 0.0, STRESS_THRESHOLD);
      EXPECT_NEAR(s(_y, _z) / s_zz, 0.0, STRESS_THRESHOLD);
      EXPECT_NEAR(s(_z, _z) / s_zz, 1.0, STRESS_THRESHOLD);
    }

    // Check displacement
    Real u_y = -nu * (s_xx + s_zz) / E;
    Real u_z = (s_zz - nu * s_xx) / E;

    for (auto && [coord, disp] : zip(
      make_view(coordinates, dim),
      make_view(displacements, dim)
    )) {
      EXPECT_NEAR(disp(_x) / u_z, 0.0, DISPLACEMENT_THRESHOLD);
      EXPECT_NEAR(disp(_y) / u_y, coord(_y), DISPLACEMENT_THRESHOLD);
      EXPECT_NEAR(disp(_z) / u_z, coord(_z), DISPLACEMENT_THRESHOLD);
    }
  }
}

TEST_F(TestPeriodicStatic, XY) {
  const double DISPLACEMENT_THRESHOLD = 1e-6;
  const double STRESS_THRESHOLD = 2e-9;

  mesh->makePeriodic(_x);
  mesh->makePeriodic(_y);
  initAll("xy");
  solve();

  Int prank = mesh->getCommunicator().whoAmI();

  if (prank == 0) {
    // Check stress
    Real s_zz = traction(_z);
    Real s_xx = nu / (1.0 - nu) * s_zz;
    Real s_yy = s_xx;

    for (auto && s : make_view(stress, dim, dim)) {
      EXPECT_NEAR(s(_x, _x) / s_xx, 1.0, STRESS_THRESHOLD);
      EXPECT_NEAR(s(_x, _y) / s_zz, 0.0, STRESS_THRESHOLD);
      EXPECT_NEAR(s(_x, _z) / s_zz, 0.0, STRESS_THRESHOLD);
      EXPECT_NEAR(s(_y, _y) / s_yy, 1.0, STRESS_THRESHOLD);
      EXPECT_NEAR(s(_y, _z) / s_zz, 0.0, STRESS_THRESHOLD);
      EXPECT_NEAR(s(_z, _z) / s_zz, 1.0, STRESS_THRESHOLD);
    }

    // Check displacement
    Real u_z = (s_zz - nu * s_xx - nu * s_yy) / E;

    for (auto && [coord, disp] : zip(
      make_view(coordinates, dim),
      make_view(displacements, dim)
    )) {
      EXPECT_NEAR(disp(_x) / u_z, 0.0, DISPLACEMENT_THRESHOLD);
      EXPECT_NEAR(disp(_y) / u_z, 0.0, DISPLACEMENT_THRESHOLD);
      EXPECT_NEAR(disp(_z) / u_z, coord(_z), DISPLACEMENT_THRESHOLD);
    }
  }
}
