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
    
    virtual void SetUp(const std::string &mesh_file, const std::string &filename) {
      mesh = std::make_unique<Mesh>(dim);
      mesh->read(mesh_file); 
      getStaticParser().parse(filename);
      model = std::make_unique<SolidMechanicsModel>(*mesh);

      const ParserSection & data = getUserParser();
      this->time_step_factor = data.getParameter("time_step_factor");

      Int normal_dir = 1;
      
      solver_ntn = std::make_unique<NTNContactSolverCallback>(*model,"slider_bottom","base_top",normal_dir,this->time_step_factor);
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

    virtual void applyFixedBoundaryConditionsInterface(){
      this->model->applyBC(BC::Dirichlet::FlagOnly(_y), "base_top");
      this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "base_top");
      this->model->applyBC(BC::Dirichlet::FlagOnly(_y), "slider_bottom");
      this->model->applyBC(BC::Dirichlet::FlagOnly(_x), "slider_bottom");  
    }

    virtual void freeBoundaryConditionsInterface(){
      this->model->applyBC(BC::Neumann::FreeBoundary(), "base_top");
      this->model->applyBC(BC::Neumann::FreeBoundary(), "slider_bottom");
      this->model->applyBC(BC::Neumann::FreeBoundary(), "slider_top");
      this->model->applyBC(BC::Neumann::FreeBoundary(), "base_bottom");
    }
    
    virtual void solveStep(){
      model->solveStep(*solver_ntn);
    }

    virtual void applyVelocity(Real shear_vel){
      UInt nb_nodes = model->getFEEngine().getMesh().getNbNodes();
      Array<Real> &position = mesh->getNodes();
      Array<Real> &velo = model->getVelocity();
      const Vector<Real> &upperBounds = mesh->getUpperBounds();
      const Vector<Real> &lowerBounds = mesh->getLowerBounds();
      Real top = upperBounds(1);
      Real bottom = lowerBounds(1);
      Real stable_time_step = model->getStableTimeStep();
      Real time_step = stable_time_step * this->time_step_factor;
      Real disp_incr = shear_vel * time_step;
      Array<Real> &displacement = model->getDisplacement();
      Array<bool> &blocked = model->getBlockedDOFs();
      // directly assign velocity on the nodes 
      for (UInt n = 0; n < nb_nodes; ++n) {
        if (std::abs(position(n,1)-top) < 1e-6) {
          velo(n,_x) = 0.5 * shear_vel;
        }
        if (std::abs(position(n,1)-bottom) < 1e-6) {
          velo(n,_x) = -0.5 * shear_vel;
        }
      }
      
      for (UInt n = 0; n < nb_nodes; ++n) {
        if (std::abs(position(n,1)-top) < 1e-6) {
          displacement(n,0) += 0.5 * disp_incr;
          blocked(n,0) = true;
        }
        if (std::abs(position(n,1)-bottom) < 1e-6) {
          displacement(n,0) += -0.5 * disp_incr;
          blocked(n,0) = true;
        }
      }
    }

    virtual void setToSteadyState(Real shear_vel){
      Array<Int> interface_node_pairs(0, 2);
      NTNContact::pairInterfaceNodes(mesh->getElementGroup("slider_bottom"),
                                    mesh->getElementGroup("base_top"), 1, *mesh,
                                    interface_node_pairs);
      const Array<Int> &sliderNodes =
      mesh->getElementGroup("slider").getNodeGroup().getNodes();
      UInt sliderNbNodes = mesh->getElementGroup("slider").getNbNodes();
      const Array<Int> &baseNodes =
      mesh->getElementGroup("base").getNodeGroup().getNodes();
      UInt baseNbNodes = mesh->getElementGroup("base").getNbNodes();
  
      // Specify initial nodal velocity
      Array<Real> &velo = model->getVelocity();
      auto friction = solver_ntn->getFriction();

      for (UInt n = 0; n < sliderNbNodes; ++n) {
        velo(sliderNodes(n), _x) = 0.5 * shear_vel;
      }
      for (UInt n = 0; n < baseNbNodes; ++n) {
        velo(baseNodes(n), _x) = -0.5 * shear_vel;
      }

      for (Int n = 0; n < (interface_node_pairs.size()); ++n) {
        friction->updateFrictionState(n, 0.00024999999968750003);
      }
    }
    
  protected:

    Real time_step_factor;
    std::unique_ptr<Mesh> mesh;
    std::unique_ptr<SolidMechanicsModel> model;
    std::unique_ptr<NTNContactSolverCallback> solver_ntn;
    
  }; // namespace akantu

}

#endif /* AKANTU_NTN_TEST_FRICTION_FIXTURE_HH_ */

