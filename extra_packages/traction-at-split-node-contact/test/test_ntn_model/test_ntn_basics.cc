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
#include "test_ntn_fixture.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>

namespace akantu{
  
  TEST_F(TestNTNFixture,TestInterfacePairing) {

    auto spatial_dimension = this->mesh->getSpatialDimension();
    
    auto & slaves = this->contact->getSlaves().getArray();
    auto & masters = this->contact->getMasters().getArray();
    const auto & position = this->mesh->getNodes();
    
    for (auto && [s,m] :
	   zip(slaves,masters)){
      
      for (Int i = 0; i < spatial_dimension;i++){
	auto distance = position(s,i)-position(m,i);	
	EXPECT_NEAR(0, distance, 1e-10);
	
      }
    }
  }; 
  
  TEST_F(TestNTNFixture,TestLumpedBoundary) {
    // Test computation of lumped boundary.
    // It computes the "surface of influence" associated
    // to an interfacial node

    const auto & l_boundary_slaves = this->ntn_contact->getLumpedBoundarySlaves();
    const auto & l_boundary_masters = this->ntn_contact->getLumpedBoundaryMasters();

    auto & slaves = this->contact->getSlaves().getArray();
    auto & masters = this->contact->getMasters().getArray();
    const auto & position = this->mesh->getNodes();
    
    const Vector<Real> lowerBounds = this->mesh->getLowerBounds();
    const Vector<Real> upperBounds = this->mesh->getUpperBounds();
    Real left = lowerBounds(0);
    Real right = upperBounds(0);
    Real system_length = right - left;    
    Int nb_contact_nodes = this->ntn_contact->getNbContactNodes();

    Real l_b_side = system_length / (nb_contact_nodes + 1);
    Real l_b_bulk = system_length / (nb_contact_nodes + 1) * 2;
    
    Real tol = 1e-11;

    for (Int n = 0; n < nb_contact_nodes; ++n){
      
      if((std::abs(position(slaves(n),_x)-left)<tol)or(std::abs(position(slaves(n),_x)-right)<tol))
	EXPECT_NEAR(l_b_side,l_boundary_slaves(n),1e-11);
      else
	EXPECT_NEAR(l_b_bulk,l_boundary_slaves(n),1e-11);
      
      if((std::abs(position(masters(n),_x)-left)<tol)or(std::abs(position(masters(n),_x)-right)<tol))
	EXPECT_NEAR(l_b_side,l_boundary_masters(n),1e-11);
      else
	EXPECT_NEAR(l_b_bulk,l_boundary_masters(n),1e-11);
      
    }
    
  };

  TEST_F(TestNTNFixture,TestImpedance) {

    Real time_step_factor = 0.1;
    this->initDynamics(time_step_factor);

    // Expected impedance is 1/ [dt/2 *(lumped_boundary[slaves]/mass[slaves] + lumped_boundary[masters]/mass[masters])]
    auto & impedance = this->ntn_contact->getImpedance();    
    Real delta_t = this->model->getTimeStep();
    
    const Vector<Real> lowerBounds = this->mesh->getLowerBounds();
    const Vector<Real> upperBounds = this->mesh->getUpperBounds();
    Real bottom = lowerBounds(1);
    Real top = upperBounds(1);
    Real system_height = top - bottom;
        
    Real density = 1;

    // Expected impedance can be reduced in this case to
    // density * system_height / (4*dt)
    Real imp_ref = density*system_height / (4*delta_t) ;
    
    for(auto && imp:
	  impedance){
      EXPECT_NEAR(imp,imp_ref,1e-11);
      
    }    
    
  };
 
}
