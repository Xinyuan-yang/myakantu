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
#include "tasn_contact_solvercallback.hh"
#include "ntn_initiation_function.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>

#ifndef AKANTU_NTN_TEST_FRICTION_FIXTURE_HH_
#define AKANTU_NTN_TEST_FRICTION_FIXTURE_HH_

// #define DEBUG_TEST

namespace akantu{
  
  class TestNTNFrictionFixture : public TestNTNFixture {
    using parent = TestNTNFixture;
  public:
    
    virtual void SetUp() override{
      parent::SetUp();
      solver_ntn = NULL;
    }	
    
    virtual void TearDown() override{
      friction.reset(nullptr);
      parent::TearDown();
      delete solver_ntn;
    }	
    virtual void initFrictionModel(){
      friction = initializeNTNFriction(*contact);      
    }
    
    virtual void initSolverCallBack(){
      solver_ntn = new TasnContactSolverCallback(*model,*contact,*friction);
    }

    virtual void solveStep(){
      model->solveStep(*solver_ntn);
    }
    
  protected:
    
    std::unique_ptr<NTNBaseFriction> friction;
    TasnContactSolverCallback *solver_ntn;
    
  }; // namespace akantu

}
#endif /* AKANTU_NTN_TEST_FRICTION_FIXTURE_HH_ */
