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
    mesh.reset(nullptr);
    model.reset(nullptr);
    solver_ntn.reset(nullptr);
  }

  virtual void applyFixedBoundaryConditions() {
    this->model->applyBC(BC::Dirichlet::FlagOnly(_y), "base_bottom");
    this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "base_bottom");
    this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "base_left");
    this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "base_right");
    this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "slider_left");
    this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "slider_right");
  }

  virtual void applyFixedBoundaryConditionsInterface() {
    this->model->applyBC(BC::Dirichlet::FlagOnly(_y), "base_top");
    // this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "base_top");
    this->model->applyBC(BC::Dirichlet::FlagOnly(_y), "slider_bottom");
    // this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "slider_bottom");
  }

  virtual void freeBoundaryConditionsInterface() {
    auto && blocked_dofs = this->model->getBlockedDOFs();
    for (auto && group :
         {"slider_top", "base_bottom", "slider_bottom", "base_top"}) {
      for (auto && node :
           this->mesh->getElementGroup(group).getNodeGroup().getNodes()) {
        blocked_dofs(node, _x) = false;
        blocked_dofs(node, _y) = false;
      }
    }

    this->model->getExternalForce().set(0);
  }

  virtual void solveStep() { model->solveStep(*solver_ntn); }

  virtual void applyVelocity(Real shear_vel) {
    UInt nb_nodes = model->getFEEngine().getMesh().getNbNodes();
    Array<Real> & position = mesh->getNodes();
    Array<Real> & velo = model->getVelocity();
    const Vector<Real> & upperBounds = mesh->getUpperBounds();
    const Vector<Real> & lowerBounds = mesh->getLowerBounds();
    Real top = upperBounds(1);
    Real bottom = lowerBounds(1);
    Real stable_time_step = model->getStableTimeStep();
    Real time_step = stable_time_step * this->time_step_factor;
    Real disp_incr = shear_vel * time_step;
    Array<Real> & displacement = model->getDisplacement();
    Array<bool> & blocked = model->getBlockedDOFs();
    // directly assign velocity on the nodes
    for (UInt n = 0; n < nb_nodes; ++n) {
      if (std::abs(position(n, 1) - top) < 1e-6) {
        velo(n, _x) = 0.5 * shear_vel;
      }
      if (std::abs(position(n, 1) - bottom) < 1e-6) {
        velo(n, _x) = -0.5 * shear_vel;
      }
    }

    for (UInt n = 0; n < nb_nodes; ++n) {
      if (std::abs(position(n, 1) - top) < 1e-6) {
        displacement(n, 0) += 0.5 * disp_incr;
        blocked(n, 0) = true;
      }
      if (std::abs(position(n, 1) - bottom) < 1e-6) {
        displacement(n, 0) += -0.5 * disp_incr;
        blocked(n, 0) = true;
      }
    }
  }

  virtual void setToSteadyState(Real shear_vel) {
    const auto & slider_nodes =
        mesh->getElementGroup("slider").getNodeGroup().getNodes();
    const auto & base_nodes =
        mesh->getElementGroup("base").getNodeGroup().getNodes();

    auto slider_nb_nodes = mesh->getElementGroup("slider").getNbNodes();
    auto base_nb_nodes = mesh->getElementGroup("base").getNbNodes();

    // Specify initial nodal velocity
    auto & velo = model->getVelocity();
    auto friction = solver_ntn->getFriction();

    for (auto n = 0; n < slider_nb_nodes; ++n) {
      velo(slider_nodes(n), _x) = 0.5 * shear_vel;
    }
    for (auto n = 0; n < base_nb_nodes; ++n) {
      velo(base_nodes(n), _x) = -0.5 * shear_vel;
    }

    auto contact = solver_ntn->getContact();

    auto vel_it =
        make_view(model->getVelocity(), model->getSpatialDimension()).begin();
    auto & slip_velocity = friction->getSlipVelocity();
    auto & slip_velocity_norm = friction->getSlipVelocityNorm();
    Real fss = friction->get("friction_state");
    for (auto && [n, master, slave, slip_vel, slip_vel_n] :
         enumerate(contact->getMasters(), contact->getSlaves(),
                   make_view(slip_velocity, slip_velocity.getNbComponent()),
                   slip_velocity_norm)) {

      slip_vel = vel_it[master] - vel_it[slave];
      slip_vel_n = slip_vel.norm();

      friction->updateFrictionState(n, fss);
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
    std::string mesh_file = "mesh-H6.0-L1.0.msh";
    std::string material_file = GetParam();
    TestNTNFrictionFixture::SetUp(mesh_file, material_file, periodic);
  }
};

} // namespace akantu

#endif /* AKANTU_NTN_TEST_FRICTION_FIXTURE_HH_ */
