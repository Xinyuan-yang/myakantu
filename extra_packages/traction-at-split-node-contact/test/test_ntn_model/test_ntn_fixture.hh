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

namespace akantu{
  
  template <typename dim_>
  class TestNTNFixture : public ::testing::Test{
  public:  
    static constexpr Int dim = dim_::value;
    
    virtual void SetUp() {
      mesh = std::make_unique<Mesh>(dim);
      mesh->read("ntn_test_" + std::to_string(dim) + "d.msh");
      std::string filename = "material_ntn_model.dat";
      model = std::make_unique<SolidMechanicsModel>(*mesh);
      ntn_contact = std::make_shared<NTNContact>(*model);     
      getStaticParser().parse(filename);
      model->initFull(SolidMechanicsModelOptions(_explicit_lumped_mass));
      model->assembleMassLumped();    

      ntn_contact->addSurfacePair("slider_bottom","base_top",1);
      contact = ntn_contact;
      
      contact->updateNormals();
      contact->updateLumpedBoundary();

    }
    
    virtual void TearDown(){
      model.reset(nullptr);
      mesh.reset(nullptr);
    }
            
    virtual void initDynamics(Real time_step_factor){
      // Set the time step and compute the impedance
      Real stable_time_step = model->getStableTimeStep();
      Real time_step = stable_time_step * time_step_factor;
      model->setTimeStep(time_step);
      contact->updateImpedance();
    }
    
    virtual void checkNodesPairing(){
      auto & slaves = this->ntn_contact->getSlaves().getArray();
      auto & masters = this->ntn_contact->getMasters().getArray();
      const auto & position = this->mesh->getNodes();
      
      for (auto && [s,m] : zip(slaves,masters)){
        for (Int d = 0; d < dim; d++){
          auto distance = position(s,d)-position(m,d);
          EXPECT_NEAR(0,distance,1e-10);  
        }      
      }
    }
    
    virtual void checkLumpedBoundary(){

      const auto & l_boundary_slaves = this->ntn_contact->getLumpedBoundarySlaves();
      const auto & l_boundary_masters = this->ntn_contact->getLumpedBoundaryMasters();

      auto & slaves = this->ntn_contact->getSlaves().getArray();
      auto & masters = this->ntn_contact->getMasters().getArray();
      const auto & position = this->mesh->getNodes();

      const Vector<Real> lowerBounds = this->mesh->getLowerBounds();
      const Vector<Real> upperBounds = this->mesh->getUpperBounds();
      Int nb_contact_nodes = this->ntn_contact->getNbContactNodes();

      Real x_min = lowerBounds(0);
      Real x_max = upperBounds(0);
      Real size_x = x_max - x_min;

      Real l_b_bulk;
      Real l_b_side;
      Real tol = 1e-11;
          
      switch(dim) {
        // 2D
        case 2:	
          l_b_bulk = size_x / (nb_contact_nodes - 1);
          l_b_side = l_b_bulk/2;
          
          for (Int n = 0; n < nb_contact_nodes; ++n){
            // Slaves
            if((std::abs(position(slaves(n),_x)-x_min)<tol)or(std::abs(position(slaves(n),_x)-x_max)<tol))
              EXPECT_NEAR(l_b_side,l_boundary_slaves(n),1e-11);
            else
              EXPECT_NEAR(l_b_bulk,l_boundary_slaves(n),1e-11);
            
            // Master
            if((std::abs(position(masters(n),_x)-x_min)<tol)or(std::abs(position(masters(n),_x)-x_max)<tol))
              EXPECT_NEAR(l_b_side,l_boundary_masters(n),1e-11);
            else
              EXPECT_NEAR(l_b_bulk,l_boundary_masters(n),1e-11);	  
          }
        break;
        
        // 3D
        case 3:
          // Same number of nodes in both directions
          Int nb_contact_per_sides = std::sqrt(nb_contact_nodes);
          
          Real z_min = lowerBounds(2);
          Real z_max = upperBounds(2);
          Real size_z = z_max - z_min;
          
          l_b_bulk = size_x * size_z / ((nb_contact_per_sides - 1)*(nb_contact_per_sides - 1));
          l_b_side = l_b_bulk / 2;
          Real l_b_corner = l_b_bulk / 4;	
          
          for (Int n = 0; n < nb_contact_nodes; ++n){
            // Slaves
            if((std::abs(position(slaves(n),_x)-x_min)<tol)or(std::abs(position(slaves(n),_x)-x_max)<tol)){
              if((std::abs(position(slaves(n),_z)-z_min)<tol)or(std::abs(position(slaves(n),_z)-z_max)<tol))
                EXPECT_NEAR(l_b_corner,l_boundary_slaves(n),1e-11);
              else
                EXPECT_NEAR(l_b_side,l_boundary_slaves(n),1e-11);
            }
            else if((std::abs(position(slaves(n),_z)-z_min)<tol)or(std::abs(position(slaves(n),_z)-z_max)<tol))
              EXPECT_NEAR(l_b_side,l_boundary_slaves(n),1e-11);
            else
              EXPECT_NEAR(l_b_bulk,l_boundary_slaves(n),1e-11);
            
            // Masters
            if((std::abs(position(masters(n),_x)-x_min)<tol)or(std::abs(position(masters(n),_x)-x_max)<tol)){
              if((std::abs(position(masters(n),_z)-z_min)<tol)or(std::abs(position(masters(n),_z)-z_max)<tol))
                EXPECT_NEAR(l_b_corner,l_boundary_masters(n),1e-11);
              else
                EXPECT_NEAR(l_b_side,l_boundary_masters(n),1e-11);
            }
            else if((std::abs(position(masters(n),_z)-z_min)<tol)or(std::abs(position(masters(n),_z)-z_max)<tol))
              EXPECT_NEAR(l_b_side,l_boundary_masters(n),1e-11);
            else
              EXPECT_NEAR(l_b_bulk,l_boundary_masters(n),1e-11);
          }	
        
        break;
      }      
    }

    virtual void checkImpedance(){
      // Expected impedance is 1/ [dt/2 *(lumped_boundary[slaves]/mass[slaves] + lumped_boundary[masters]/mass[masters])]
      auto & impedance = this->ntn_contact->getImpedance();    
      Real delta_t = this->model->getTimeStep();
      
      const Vector<Real> lowerBounds = this->mesh->getLowerBounds();
      const Vector<Real> upperBounds = this->mesh->getUpperBounds();
      Real bottom = lowerBounds(1);
      Real top = upperBounds(1);
      Real system_height = top - bottom;
      
      Real density = 1;
      
      // Expected impedance can be reduced to
      // density * system_height / (4*dt)
      Real imp_ref = density*system_height / (4*delta_t) ;
      
      for(auto && imp:
	    impedance){
	      EXPECT_NEAR(imp,imp_ref,1e-11);
      }          
    }
    
  protected:
    std::unique_ptr<Mesh> mesh;
    std::unique_ptr<SolidMechanicsModel> model;
    std::shared_ptr<NTNBaseContact> contact;
    std::shared_ptr<NTNContact> ntn_contact;
    
  }; // namespace akantu
  
  using dim_types = testing::Types<dim_2_t,dim_3_t>;
  //using dim_types = gtest_list_t<TestAllDimensions>;

TYPED_TEST_SUITE(TestNTNFixture, dim_types, );

}
  
#endif /* AKANTU_NTN_TEST_FIXTURE_HH_ */
