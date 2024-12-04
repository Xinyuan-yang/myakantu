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
#include "test_ntn_friction_fixture.hh"

/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>

namespace akantu{
  
  TEST_F(TestNTNFrictionFixture,TestCompression) {
    
    std::string mesh_file = "ntn_test_compression.msh";
    this->initMesh(mesh_file);
    
    std::string filename = "material_ntn_friction_slip_weakening.dat";    
    this->initSolidMechanicsModel(filename);
    this->initContactModel();

    Vector<Real> trac = {0,-1e6};

    Real time_factor = 0.1;

    this->initDynamics(time_factor);

    // Get time step
    Real stable_time_step = this->model->getStableTimeStep();
    Real time_step = stable_time_step * time_factor;
    
    this->model->assembleMassLumped();
    
    // Set the boundary conditions
    
    this->model->applyBC(BC::Dirichlet::FlagOnly(_y), "base_bottom");
    this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "base_bottom");
    this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "base_left");
    this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "base_right");
    this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "slider_left");
    this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "slider_right");
    
    this->initFrictionModel();
    this->initSolverCallBack();

    this->solver_ntn->assembleResidual();

    // Time of the simulation
    Real simulation_time = 5e-4;
    UInt nb_steps = UInt(simulation_time / time_step + 1);
    
    // Period of the sinusoidale compressive wave
    Real period = 5e-5;

    Vector<Real> trac_t = trac;
    Vector<Real> trac_nul = {0,0};

    const auto & pressure = this->contact->getContactPressure();
    UInt nb_contact_nodes = this->contact->getNbContactNodes();
    
    for (UInt s = 0; s < nb_steps; ++s) {
      
      if ((s*time_step) < (period/2))
	{
	  trac_t = trac * std::sin(s*time_step/period*2*M_PI);
	  this->model->applyBC(BC::Neumann::FromTraction(trac_t), "slider_top");
	}
      else
	this->model->applyBC(BC::Neumann::FromTraction(trac_nul), "slider_top");
      
      this->solveStep();
      
      if (s%5==0){
	std::cout << "Simulation at t " << s*time_step << ". Contact pressure:";
	
	for (UInt n = 0; n < nb_contact_nodes; ++n){	  
	  std::cout << " " << pressure(n,_y);
	}
	std::cout << " " << std::endl;
      }
    } 
    
  };
  
}
  
