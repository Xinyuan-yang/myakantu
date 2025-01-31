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
#include "mesh_utils.hh"
#include "non_linear_solver.hh"
#include "parser.hh"
#include "solid_mechanics_model.hh"
#include "test_mesh_periodic_fixture.hh"
#include <gtest/gtest.h>

#ifndef AKANTU_TEST_PERIODIC_STATIC_FIXTURE_HH_
#define AKANTU_TEST_PERIODIC_STATIC_FIXTURE_HH_

using namespace akantu;

template <class Params>
class TestPeriodicStatic : public TestMeshPeriodic<Params> {
  using parent = TestMeshPeriodic<Params>;

public:
  void SetUp() override {
    parent::SetUp();
    getStaticParser().parse("material.dat");
    element_type =
        *(this->mesh->elementTypes(_spatial_dimension = parent::dim).begin());
  }

  void TearDown() override {
    model.reset(nullptr);
    parent::TearDown();
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

    model->addDumpField("grad_u");
    model->addDumpField("stress");
    model->dump();
  }

protected:
  ElementType element_type;
  std::unique_ptr<SolidMechanicsModel> model;
  Vector<Real, parent::dim> traction;
  Real E;
  Real nu;

  void initModel() {
    model = std::make_unique<SolidMechanicsModel>(*this->mesh);
    model->initFull(_analysis_method = _static);
    auto & material = model->getMaterial(0);
    E = material.getParam("E");
    nu = material.getParam("nu");
  }

  void initDumper(std::string name_suffix) {
    model->setBaseName("periodic-static-" + name_suffix);
    model->addDumpFieldVector("displacement");
    model->addDumpField("external_force");
    model->addDumpField("internal_force");
    model->addDumpField("blocked_dofs");
  }

  void applyBC() {
    for (auto sp : spatial_directions<parent::dim>) {
      model->applyBC(BC::Dirichlet::FixedValue(0.0, sp), "corner");
    }

    traction.fill(0.);
    traction(parent::dim - 1) = 1e9;

    if (parent::dim == 1) {
      for (auto && node :
           this->mesh->getElementGroup("x-plus").getNodeGroup().getNodes()) {
        model->getExternalForce()(node) = traction(_x);
      }
      this->mesh->makeConsistentForPeriodicity(model->getExternalForce());
    }

    if (parent::dim == 2) {
      model->applyBC(BC::Dirichlet::FixedValue(0.0, _y), "y-minus");
      model->applyBC(BC::Neumann::FromTraction(traction), "y-plus");
    }

    if (parent::dim == 3) {
      model->applyBC(BC::Dirichlet::FixedValue(0.0, _y), "line");
      model->applyBC(BC::Dirichlet::FixedValue(0.0, _z), "z-minus");
      model->applyBC(BC::Neumann::FromTraction(traction), "z-plus");
    }
  }

  void initSolver() {
    auto & solver = model->getNonLinearSolver();
    solver.set("max_iterations", 5);
    solver.set("threshold", 2e-4);
    solver.set("convergence_type", SolveConvergenceCriteria::_solution);
  }
};

#endif /* AKANTU_TEST_PERIODIC_STATIC_FIXTURE_HH_ */
