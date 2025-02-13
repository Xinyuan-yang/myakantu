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

namespace akantu {

TEST_F(TestNTNFrictionFixture, TestCompression) {
  SetUp("ntn_test_compression.msh", "material_ntn_friction_slip_weakening.dat");

  std::cout << std::fixed;
  std::cout << std::setprecision(12);

  Vector<Real> ref_pressure = {0.0,
                               0.0,
                               0.0,
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

  // Create the list of compression values at several time step and then EXPECT
  // NEAR
  Vector<Real> trac = {0, -1e6};

  // Set the boundary conditions
  this->applyFixedBoundaryConditions();

  // Time of the simulation
  Real simulation_time = 5e-4;
  Real time_step = this->model->getTimeStep();
  UInt nb_steps = UInt(simulation_time / time_step + 1);

  // Period of the sinusoidale compressive wave
  Real period = 5e-5;

  Vector<Real> trac_t = trac;
  Vector<Real> trac_nul = {0, 0};

  const auto & pressure = this->solver_ntn->getContactPressure();

  Int check_counter = 0;
  for (UInt s = 0; s < nb_steps; ++s) {
    if ((s * time_step) < (period / 2)) {
      trac_t = trac * std::sin(s * time_step / period * 2 * M_PI);
      this->model->applyBC(BC::Neumann::FromTraction(trac_t), "slider_top");
    } else
      this->model->applyBC(BC::Neumann::FromTraction(trac_nul), "slider_top");

    this->solveStep();

    if (s % 5 == 0) {
      EXPECT_NEAR(pressure(ref_contact_node, pressure_direction),
                  ref_pressure(check_counter), 5e-9);
      check_counter += 1;
    }
  }
}

TEST_P(ParamTestNTNFrictionFixture, TestRateAndState) {
  std::string material_file = GetParam();
  SetUp("ntn_test_ras.msh", material_file, true);

  std::cout << std::fixed;
  std::cout << std::setprecision(12);

  const ParserSection & data = getUserParser();
  const auto & mat = model->getMaterial("slider");

  Real shear_vel = data.getParameter("shear_velocity");
  Vector<Real> trac_top = data.getParameter("top_traction");
  Vector<Real> trac_bottom = data.getParameter("bot_traction");

  // Static analytical solution
  Real fss = data.getParameter("fss");
  Real E = mat.getParam("E");
  Real nu = mat.getParam("nu");  
  Real shear_modulus = E / (2. * (1. + nu));
  Real normal_strain_applied = trac_top(1) / E - nu * nu * trac_top(1) / E;
  
  Array<Real> & displacement = model->getDisplacement();
  Array<Real> & position = mesh->getNodes();
  UInt nb_nodes = model->getFEEngine().getMesh().getNbNodes();

  for (UInt n = 0; n < nb_nodes; ++n) {
    displacement(n, 0) = fss * -trac_top(1) / (shear_modulus) * position(n, 1);
    displacement(n, 1) = normal_strain_applied * position(n, 1);
  }

  // Set boundary conditions for dynamic simulation
  this->model->applyBC(BC::Neumann::FromTraction(trac_top), "slider_top");
  this->model->applyBC(BC::Neumann::FromTraction(trac_bottom), "base_bottom");

  // Set a steady state sliding
  this->setToSteadyState(shear_vel);

  // Time of the simulation
  Real stable_time_step = this->model->getStableTimeStep();
  Real time_step = stable_time_step * this->time_step_factor;
  model->setTimeStep(time_step);
  UInt nb_steps = data.getParameter("nb_steps");
  
  const Array<Int> & base_nodes =
       mesh->getElementGroup("base").getNodeGroup().getNodes();
  UInt base_nb_nodes = mesh->getElementGroup("base").getNbNodes();
  Array<Real> & velo = model->getVelocity();

  for (UInt s = 0; s < nb_steps; ++s) {
    this->applyVelocity(shear_vel);
    model->solveStep(*solver_ntn, "explicit_lumped");

    if (s % 100 == 0) {
      for (UInt n = 0; n < base_nb_nodes; ++n) {
        EXPECT_NEAR(velo(base_nodes(n), _x), -shear_vel / 2, 1e-10);
      } 
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    MaterialFileTests, ParamTestNTNFrictionFixture,
      ::testing::Values(
          "material_ntn_friction_rate_and_state_dieterich_aging_regularized.dat",
          "material_ntn_friction_rate_and_state_dieterich_aging.dat",
          "material_ntn_friction_rate_and_state_dieterich_slip.dat",
          "material_ntn_friction_rate_and_state_revised_dieterich_aging_regularized.dat",
          "material_ntn_friction_rate_and_state_revised_dieterich_aging.dat",
          "material_ntn_friction_rate_and_state_revised_dieterich_slip.dat",
          "material_ntn_friction_rate_and_state_regularized_dieterich_aging_regularized.dat",
          "material_ntn_friction_rate_and_state_regularized_dieterich_aging.dat",
          "material_ntn_friction_rate_and_state_regularized_dieterich_slip.dat",
          "material_ntn_friction_rate_and_state_n_shaped_aging_regularized.dat",
          "material_ntn_friction_rate_and_state_n_shaped_aging.dat",
          "material_ntn_friction_rate_and_state_n_shaped_slip.dat",
          "material_ntn_friction_rate_and_state_revised_n_shaped_aging_regularized.dat",
          "material_ntn_friction_rate_and_state_revised_n_shaped_aging.dat",
          "material_ntn_friction_rate_and_state_revised_n_shaped_slip.dat")); 
/*       ::testing::Values(
          "material_ntn_friction_rate_and_state_dieterich_aging.dat"));  */

/*   TEST_F(TestNTNFrictionFixture,TestRateAndStateHealing) {
    SetUp("mesh-H6.0-L1.0.msh", "material_ntn_friction_rate_and_state.dat");

    const ParserSection & data = getUserParser();
    Real trac_top = data.getParameter("top_traction"); //{0.0,-1e6};
    Real trac_bottom = data.getParameter("bot_traction"); //{-0.0, 1e6};
    this->model->applyBC(BC::Neumann::FromTraction(trac_top), "slider_top");
    this->model->applyBC(BC::Neumann::FromTraction(trac_bottom), "base_bottom");
    UInt nb_steps = 10;
    auto friction = solver_ntn->getFriction();

    for (UInt s = 0; s < nb_steps; ++s) {
      model->solveStep("explicit_lumped");
      //std::cout << "friction_state = " << friction->friction_state <<
  std::endl;
    }
  } */

} // namespace akantu
