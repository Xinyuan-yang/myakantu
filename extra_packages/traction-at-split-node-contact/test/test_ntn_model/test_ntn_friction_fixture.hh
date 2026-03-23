/**
 * Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "mesh_utils.hh"
#include "ntn_base_contact.hh"
#include "ntn_base_friction.hh"
#include "ntn_contact.hh"
#include "ntn_contact_solvercallback.hh"
#include "ntn_initiation_function.hh"
#include "solid_mechanics_model.hh"
#include "test_ntn_fixture.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>

#ifndef AKANTU_NTN_TEST_FRICTION_FIXTURE_HH_
#define AKANTU_NTN_TEST_FRICTION_FIXTURE_HH_

namespace akantu {

class Velocity : public BC::Dirichlet::DirichletFunctor {
public:
  explicit Velocity(SolidMechanicsModel & model, Real vel, BC::Axis ax = _x)
      : DirichletFunctor(ax), model(model), vel(vel) {
    disp = vel * model.getTimeStep();
  }

  inline void operator()(Idx node, VectorProxy<bool> & flags,
                         VectorProxy<Real> & disp,
                         const VectorProxy<const Real> & /* coord*/) override {
    disp(axis) += this->disp;
    model.getVelocity()(node, axis) = vel;
    flags(axis) = true;
  }

private:
  SolidMechanicsModel & model;
  Real vel, disp;
};

class FreeBC : public BC::Dirichlet::DirichletFunctor {
public:
  explicit FreeBC(BC::Axis ax = _x) : DirichletFunctor(ax) {}

  inline void operator()(Idx /*&node*/, VectorProxy<bool> & flags,
                         VectorProxy<Real> & /*disp*/,
                         const VectorProxy<const Real> & /* coord*/) override {
    flags(axis) = false;
  }
};

class TestNTNFrictionFixture : public ::testing::Test {
public:
  static constexpr Int dim = 2;

  virtual void SetUp(const std::string & mesh_file,
                     const std::string & filename, bool periodic = false) {
    mesh = std::make_unique<Mesh>(dim);
    mesh->read(mesh_file);
    getStaticParser().parse(filename);

    if (periodic) {
      mesh->makePeriodic(_x, "slider_left", "slider_right");
      mesh->makePeriodic(_x, "base_left", "base_right");
    }

    model = std::make_unique<SolidMechanicsModel>(*mesh);

    const ParserSection & data = getUserParser();
    this->time_step_factor = data.getParameter("time_step_factor");

    Int normal_dir = 1;

    solver_ntn = std::make_unique<NTNContactSolverCallback>(
        *model, "slider_bottom", "base_top", normal_dir,
        this->time_step_factor);
  }

  virtual void TearDown() override {
    solver_ntn.reset(nullptr);
    model.reset(nullptr);
    mesh.reset(nullptr);
  }

  virtual void applyFixedBoundaryConditions() {
    this->model->applyBC(BC::Dirichlet::FlagOnly(_y), "base_bottom");
    this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "base_bottom");
    this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "base_left");
    this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "base_right");
    this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "slider_left");
    this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "slider_right");
  }

  virtual void solveStep() { model->solveStep(*solver_ntn); }

  virtual void applyVelocity(Real shear_vel) {
    Velocity velbc_top(*model, .5 * shear_vel);
    Velocity velbc_bottom(*model, -.5 * shear_vel);

    this->model->applyBC(velbc_top, "slider_top");
    this->model->applyBC(velbc_bottom, "base_bottom");
  }

  virtual void setToSteadyState(Real shear_vel) {
    const auto & slider_nodes =
        mesh->getElementGroup("slider").getNodeGroup().getNodes();
    const auto & base_nodes =
        mesh->getElementGroup("base").getNodeGroup().getNodes();

    // Specify initial nodal velocity
    auto & velo = model->getVelocity();
    auto & increment = model->getIncrement();
    auto friction = solver_ntn->getFriction();
    auto dt = model->getTimeStep();

    for (auto n : slider_nodes) {
      velo(n, _x) = 0.5 * shear_vel;
      increment(n, _x) = 0.5 * shear_vel * dt;
    }
    for (auto n : base_nodes) {
      velo(n, _x) = -0.5 * shear_vel;
      increment(n, _x) = -0.5 * shear_vel * dt;
    }

    auto contact = solver_ntn->getContact();

    auto vel_it =
        make_view(model->getVelocity(), model->getSpatialDimension()).begin();
    auto & slip_velocity = friction->getSlipVelocity();
    auto & slip_velocity_norm = friction->getSlipVelocityNorm();
    auto & is_sticking = friction->getIsSticking();
    Real phi = friction->get("friction_state");

    for (auto && [n, master, slave, slip_vel, slip_vel_n, is_sticking] :
         enumerate(contact->getMasters(), contact->getSlaves(),
                   make_view(slip_velocity, slip_velocity.getNbComponent()),
                   slip_velocity_norm, is_sticking)) {

      is_sticking = false;

      friction->updateFrictionState(n, phi);
    }
  }

protected:
  Real time_step_factor;
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<SolidMechanicsModel> model;
  std::unique_ptr<NTNContactSolverCallback> solver_ntn;
};

class ParamTestNTNFrictionFixture
    : public TestNTNFrictionFixture,
      public ::testing::WithParamInterface<std::string> {
protected:
  void SetUp(const std::string & /*mesh_file*/,
             const std::string & /*filename*/, bool periodic = false) override {
    std::string mesh_file = "ntn_test_ras.msh";
    std::string material_file = GetParam();
    TestNTNFrictionFixture::SetUp(mesh_file, material_file, periodic);
  }
};

} // namespace akantu

#endif /* AKANTU_NTN_TEST_FRICTION_FIXTURE_HH_ */
