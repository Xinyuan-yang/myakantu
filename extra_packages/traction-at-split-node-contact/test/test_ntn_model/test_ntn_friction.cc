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
#include <iomanip>
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>

namespace akantu{
  
  TEST_F(TestNTNFrictionFixture,TestCompression) {

    std::cout << std::fixed;
    std::cout << std::setprecision(12);
    
    Vector<Real> ref_pressure = {0.0, 0.0, 0.0,
      -0.000000000061,
      -0.000000072362,
      -0.000011543530,
      -0.000611366799,
      -0.015767551005,
      -0.243299242707,
      -2.545124342651,
      -19.590589983818,
      -117.424265038805,
      -570.863394812224,
      -2319.934641685013,
      -8063.335008171387,
      -24394.402146997298,
      -65124.820990435350,
      -155070.572077111341,
      -332093.403606061882,
      -643767.301224190160,
      -1135088.793070584768,
      -1826677.261436731322,
      -2689111.493139004800,
      -3625981.040552595165,
      -4480965.840503907762,
      -5077927.732134595513,
      -5286961.279143503867,
      -5090152.021338184364};

    // The reference pressure correspond to the contact node 1
    Int ref_contact_node = 1;
    // The pressure is in the direction 1, normal to the surface
    Int pressure_direction = 1;
    
    // Create the list of compression values at several time step and then EXPECT NEAR
    
    Vector<Real> trac = {0,-1e6};
    
    // Set the boundary conditions
    this->applyFixedBoundaryConditions();
          
    // Time of the simulation
    Real simulation_time = 5e-4;
    Real time_step = this->model->getTimeStep();
    UInt nb_steps = UInt(simulation_time / time_step + 1);
    
    // Period of the sinusoidale compressive wave
    Real period = 5e-5;

    Vector<Real> trac_t = trac;
    Vector<Real> trac_nul = {0,0};

    const auto & pressure = this->solver_ntn->getContactPressure();
    
    Int check_counter = 0;
    
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
	EXPECT_NEAR(pressure(ref_contact_node,pressure_direction),ref_pressure(check_counter),1e-11);	
	check_counter += 1;
      }      
    }
  }
}
