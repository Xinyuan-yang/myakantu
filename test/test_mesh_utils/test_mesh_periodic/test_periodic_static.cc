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
#include "dof_synchronizer.hh"
#include "test_periodic_static_fixture.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
/* -------------------------------------------------------------------------- */

using TestDIMs =
    ::testing::Types<std::tuple<dim1, dim1>, std::tuple<dim2, dim1>,
                     std::tuple<dim3, dim0>, std::tuple<dim3, dim1>,
                     std::tuple<dim3, dim2>>;

TYPED_TEST_SUITE(TestPeriodicStatic, TestDIMs);

const Real displacement_threshold = 1e-6;
const Real stress_threshold = 2e-9;

/* -------------------------------------------------------------------------- */
template <Int dim, Int dim_periodic> struct ValidationHelper {
  static auto stressRef(const Vector<Real> & traction, Real /*nu*/)
      -> Matrix<Real> {
    Matrix<Real> stress(dim, dim);
    stress(_x, _x) = traction(_x) / 2;

    return stress;
  }

  static auto strainRef(const Vector<Real> & traction, Real E, Real /*nu*/)
      -> Vector<Real> {
    Vector<Real> strain(dim);
    strain(_x) = traction(_x) / 2. / E;

    return strain;
  }
};

/* -------------------------------------------------------------------------- */
template <> struct ValidationHelper<2, 1> {
  static auto stressRef(const Vector<Real> & traction, Real nu)
      -> Matrix<Real> {
    Matrix<Real> stress(2, 2);
    stress.fill(0.);

    Real s_yy = traction(_y);
    Real s_xx = nu / (1.0 - nu) * s_yy;

    stress(_x, _x) = s_xx;
    stress(_y, _y) = s_yy;
    return stress;
  }

  static auto strainRef(const Vector<Real> & traction, Real E, Real nu)
      -> Vector<Real> {
    Vector<Real> strain(2);

    Real s_yy = traction(_y);
    Real s_xx = nu / (1.0 - nu) * s_yy;
    Real u_y = (1 + nu) * (-nu * s_xx + (1 - nu) * s_yy) / E;

    strain(_x) = 0.;
    strain(_y) = u_y;

    return strain;
  }
};

/* -------------------------------------------------------------------------- */
template <> struct ValidationHelper<3, 0> {
  static auto stressRef(const Vector<Real> & traction, Real /*nu*/)
      -> Matrix<Real> {
    Matrix<Real> stress(3, 3);
    stress.fill(0.);

    stress(_z, _z) = traction(_z);
    return stress;
  }

  static auto strainRef(const Vector<Real> & traction, Real E, Real nu)
      -> Vector<Real> {
    Vector<Real> strain(3);

    Real u_z = traction(_z) / E;
    Real u_xy = -nu * u_z;

    strain(_x) = strain(_y) = u_xy;
    strain(_z) = u_z;

    return strain;
  }
};

/* -------------------------------------------------------------------------- */
template <> struct ValidationHelper<3, 1> {
  static auto stressRef(const Vector<Real> & traction, Real nu)
      -> Matrix<Real> {
    Matrix<Real> stress(3, 3);
    stress.fill(0.);

    Real s_zz = traction(_z);
    Real s_xx = nu * s_zz;

    stress(_x, _x) = s_xx;
    stress(_z, _z) = s_zz;
    return stress;
  }

  static auto strainRef(const Vector<Real> & traction, Real E, Real nu)
      -> Vector<Real> {
    Vector<Real> strain(3);

    Real s_zz = traction(_z);
    Real s_xx = nu * s_zz;

    Real u_y = -nu * (s_xx + s_zz) / E;
    Real u_z = (s_zz - nu * s_xx) / E;

    strain(_x) = 0.;
    strain(_y) = u_y;
    strain(_z) = u_z;

    return strain;
  }
};

/* -------------------------------------------------------------------------- */
template <> struct ValidationHelper<3, 2> {
  static auto stressRef(const Vector<Real> & traction, Real nu)
      -> Matrix<Real> {
    Matrix<Real> stress(3, 3);
    stress.fill(0.);

    Real s_zz = traction(_z);
    Real s_xx = nu / (1.0 - nu) * s_zz;
    Real s_yy = s_xx;

    stress(_x, _x) = s_xx;
    stress(_y, _y) = s_yy;
    stress(_z, _z) = s_zz;
    return stress;
  }

  static auto strainRef(const Vector<Real> & traction, Real E, Real nu)
      -> Vector<Real> {
    Vector<Real> strain(3);

    Real s_zz = traction(_z);
    Real s_xx = nu / (1.0 - nu) * s_zz;
    Real s_yy = s_xx;

    Real u_z = (s_zz - nu * s_xx - nu * s_yy) / E;

    strain(_x) = strain(_y) = 0.;
    strain(_z) = u_z;

    return strain;
  }
};

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

TYPED_TEST(TestPeriodicStatic, Traction) {
  std::cout << "Seed: " << RandomGenerator<Idx>::seed() << std::endl;

  Int psize = this->mesh->getCommunicator().getNbProc();
  Int prank = this->mesh->getCommunicator().whoAmI();
  this->initAll(std::to_string(this->dim) +
                "d-pbc:" + std::to_string(this->dim_periodic) +
                "-proc:" + std::to_string(psize));
  this->solve();

  const auto & coordinates = this->mesh->getNodes();
  const auto & displacements = this->model->getDisplacement();
  const auto & stress =
      this->model->getMaterial(0).getStress(this->element_type);

  // Check stress
  auto && s_ref = ValidationHelper<this->dim, this->dim_periodic>::stressRef(
      this->traction, this->nu);

  for (auto && [e, s] : enumerate(make_view(stress, this->dim, this->dim))) {
    auto error = (s - s_ref).template lpNorm<Eigen::Infinity>() / s_ref.norm();

    if (this->dim == 1) {
      Matrix<Real> diff = s.array().abs() - s_ref.array();
      error = diff.template lpNorm<Eigen::Infinity>() / s_ref.norm();
    }

    EXPECT_NEAR(error, 0., stress_threshold)
        << "On proc " << prank << " for element " << e << " -  " << s
        << " != " << s_ref;
  }

  auto lower = this->mesh->getLowerBounds();

  // Check displacement
  auto strain = ValidationHelper<this->dim, this->dim_periodic>::strainRef(
      this->traction, this->E, this->nu);

  for (auto && [idx, coord, disp] :
       enumerate(make_view(coordinates, this->dim),
                 make_view(displacements, this->dim))) {
    if (not this->mesh->isLocalOrMasterNode(idx)) {
      continue;
    }

    Vector<Real> diff(this->dim);
    if (this->dim == 1) {
      diff.array() =
          (disp.array() / strain.array()).abs() - coord.array().abs();
    } else {
      diff.array() =
          disp.array() - (coord.array() - lower.array()) * strain.array();
    }
    auto error = diff.template lpNorm<Eigen::Infinity>() / strain.norm();
    EXPECT_NEAR(error, 0., displacement_threshold)
        << disp.array() << " != " << coord.array() * strain.array();
  }
}
