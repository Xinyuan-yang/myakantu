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
    std::string filename = "material_ntn_model.dat";

    this->initSolidMechanicsModel(filename);
    this->initContactModel();

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
}
