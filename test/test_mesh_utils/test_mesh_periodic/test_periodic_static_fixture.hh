/**
 * Copyright (©) 2018-2024 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#include "parser.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
#include "non_linear_solver.hh"
#include <gtest/gtest.h>

#ifndef AKANTU_TEST_PERIODIC_STATIC_FIXTURE_HH_
#define AKANTU_TEST_PERIODIC_STATIC_FIXTURE_HH_

using namespace akantu;


class TestPeriodicStatic : public ::testing::Test {
public:
  void SetUp() override {
    const size_t dim = 3;
    debug::setDebugLevel(dblError);
    getStaticParser().parse("material.dat");

    const auto & comm = Communicator::getStaticCommunicator();
    Int psize = comm.getNbProc();
    Int prank = comm.whoAmI();
    AKANTU_DEBUG_INFO("Running on " << psize << " processors");
    AKANTU_DEBUG_INFO("This is processor " << prank);

    mesh = std::make_unique<Mesh>(dim);
    if (prank == 0) {
      mesh->read("cube.msh");
    }
    mesh->distribute();
  }

  void TearDown() override {
    mesh.reset(nullptr);
  }

  void initAll(std::string name_suffix) {
    initModel();
    applyBC();
    initSolver();
    initDumper(name_suffix);
  }

  void solve() {
    model->dump();
    model->solveStep();
    model->dump();
  }

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<SolidMechanicsModel> model;

  void initModel() {
    model = std::make_unique<SolidMechanicsModel>(*mesh);
    model->initFull(_analysis_method = _static);
  }

  void initDumper(std::string name_suffix) {
    model->setBaseName("periodic-static-" + name_suffix);
    model->addDumpFieldVector("displacement");
    model->addDumpField("external_force");
    model->addDumpField("internal_force");
    model->addDumpField("blocked_dofs");
    model->addDumpField("grad_u");
    model->addDumpField("stress");
  }

  void applyBC() {
    model->applyBC(BC::Dirichlet::FixedValue(0.0, _x), "corner");
    model->applyBC(BC::Dirichlet::FixedValue(0.0, _y), "corner");
    model->applyBC(BC::Dirichlet::FixedValue(0.0, _y), "line");
    model->applyBC(BC::Dirichlet::FixedValue(0.0, _z), "bottom");

    Vector<Real, 3> surface_traction({1.0, 2.0, 3.0});
    model->applyBC(BC::Neumann::FromTraction(surface_traction), "top");
  }
  
  void initSolver() {
    auto & solver = model->getNonLinearSolver();
    solver.set("max_iterations", 2);
    solver.set("threshold", 2e-4);
    solver.set("convergence_type", SolveConvergenceCriteria::_solution);
  }
};

#endif /* AKANTU_TEST_PERIODIC_STATIC_FIXTURE_HH_ */
