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
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
#include "ntn_contact.hh"
#include "ntn_base_friction.hh"
#include "ntn_base_contact.hh"
#include "test_ntn_fixture.hh"
#include "ntn_contact_solvercallback.hh"
#include "ntn_initiation_function.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>

#ifndef AKANTU_NTN_TEST_FRICTION_FIXTURE_HH_
#define AKANTU_NTN_TEST_FRICTION_FIXTURE_HH_

namespace akantu{
  
  class TestNTNFrictionFixture : public ::testing::Test {
  public:
    static constexpr Int dim = 2;
    
    virtual void SetUp() override{
      mesh = std::make_unique<Mesh>(dim);
      std::string mesh_file = "ntn_test_compression.msh";
      mesh->read(mesh_file);
      std::string filename = "material_ntn_friction_slip_weakening.dat";    
      getStaticParser().parse(filename);
      model = std::make_unique<SolidMechanicsModel>(*mesh);

      Real time_step_factor = 0.1;

      Int normal_dir = 1;
      
      solver_ntn = std::make_unique<NTNContactSolverCallback>(*model,"slider_bottom","base_top",normal_dir,time_step_factor);
      
    }	
    
    virtual void TearDown() override{
      mesh.reset(nullptr);
      model.reset(nullptr);
      solver_ntn.reset(nullptr);
    }	

    virtual void applyFixedBoundaryConditions(){
    
      this->model->applyBC(BC::Dirichlet::FlagOnly(_y), "base_bottom");
      this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "base_bottom");
      this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "base_left");
      this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "base_right");
      this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "slider_left");
      this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "slider_right");

    }
    
    virtual void solveStep(){
      model->solveStep(*solver_ntn);
    }
    
  protected:

    std::unique_ptr<Mesh> mesh;
    std::unique_ptr<SolidMechanicsModel> model;
    std::unique_ptr<NTNContactSolverCallback> solver_ntn;
    
  }; // namespace akantu

}

#endif /* AKANTU_NTN_TEST_FRICTION_FIXTURE_HH_ */

