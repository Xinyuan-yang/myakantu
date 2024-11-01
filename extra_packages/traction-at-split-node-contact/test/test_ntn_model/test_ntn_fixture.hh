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
#include "ntn_base_contact.hh"

/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>


#ifndef AKANTU_NTN_TEST_FIXTURE_HH_
#define AKANTU_NTN_TEST_FIXTURE_HH_

// #define DEBUG_TEST

namespace akantu{
  
  class TestNTNFixture : public ::testing::Test{
  public:  
    virtual void SetUp() {
      UInt dim = 2;
      mesh = std::make_unique<Mesh>(dim);
      mesh->read("ntn_test_2d.msh");      
      model = std::make_unique<SolidMechanicsModel>(*mesh);
      ntn_contact = new NTNContact(*model);
    
    }
    
    virtual void TearDown(){
      model.reset(nullptr);
      mesh.reset(nullptr);
      delete ntn_contact;
    }
    
    virtual void initSolidMechanicsModel(const std::string & material_file){

      getStaticParser().parse(material_file);

      // Create model      
      model->initFull(SolidMechanicsModelOptions(_explicit_lumped_mass));
      model->assembleInternalForces();
      for (Int i = 0; i < model->getNbMaterials(); ++i) {
	model->getMaterial(i).setToSteadyState(_not_ghost);
	model->getMaterial(i).setToSteadyState(_ghost);
      }
      
      model->assembleInternalForces();
      model->assembleMassLumped();
      model->assembleInternalForces();
    }
    
    virtual void initContactModel(){
      
    // Create the contact and node to node interface
      
      contact = NULL;
      ntn_contact->addSurfacePair("slider_bottom","base_top",1);
      contact = ntn_contact;
      
      contact->updateNormals();
      contact->updateLumpedBoundary();
      contact->updateImpedance();
      
    }
    
    
  protected:
    std::unique_ptr<Mesh> mesh;
    std::unique_ptr<SolidMechanicsModel> model;
    NTNContact *contact;
    NTNContact *ntn_contact;
    
  }; // namespace akantu
  
}
#endif /* AKANTU_NTN_TEST_FIXTURE_HH_ */

