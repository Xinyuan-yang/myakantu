/**
 * @file   phase_field_parallel.cc
 *
 * @author Mohit Pundir <mohit.pundir@ethz.ch>
 *
 * @date creation: Mon May 09 2022
 * @date last modification: Mon May 09 2022
 *
 * @brief  Example of phase field model in parallel
 *
 *
 * @section LICENSE
 *
 * Copyright (©) 2018-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
 * Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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
 *
 */

/* -------------------------------------------------------------------------- */
#include "coupler_solid_phasefield.hh"
/* -------------------------------------------------------------------------- */
#include <chrono>
#include <fstream>
#include <iostream>
/* -------------------------------------------------------------------------- */

using namespace akantu;
using clk = std::chrono::high_resolution_clock;
using second = std::chrono::duration<double>;
using millisecond = std::chrono::duration<double, std::milli>;

const UInt spatial_dimension = 2;

int main(int argc, char * argv[]) {

  initialize("material_notch.dat", argc, argv);

  // create mesh
  Mesh mesh(spatial_dimension);

  const auto & comm = Communicator::getStaticCommunicator();
  Int prank = comm.whoAmI();
  if (prank == 0) {
    // Read the mesh
    mesh.read("square_notch.msh");
  }

  mesh.distribute();

  CouplerSolidPhaseField coupler(mesh);
  auto & model = coupler.getSolidMechanicsModel();
  auto & phase = coupler.getPhaseFieldModel();

  model.initFull(_analysis_method = _static);
  auto && mat_selector =
      std::make_shared<MeshDataMaterialSelector<std::string>>("physical_names",
                                                              model);
  model.setMaterialSelector(mat_selector);

  auto && selector = std::make_shared<MeshDataPhaseFieldSelector<std::string>>(
      "physical_names", phase);
  phase.setPhaseFieldSelector(selector);

  phase.initFull(_analysis_method = _static);

  model.applyBC(BC::Dirichlet::FixedValue(0., _y), "bottom");
  model.applyBC(BC::Dirichlet::FixedValue(0., _x), "left");

  model.setBaseName("phase_notch_parallel");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpFieldVector("displacement");
  model.addDumpField("damage");
  if (mesh.isDistributed()) {
    phase.addDumpField("partitions");
  }
  phase.dump();

  UInt nbSteps = 1000;
  Real increment = 4e-6;
  UInt nb_staggered_steps = 10;

  Int nbNodes = mesh.getNbNodes();
  std::cout << "Number of nodes: " << nbNodes << std::endl;

  auto start_time = clk::now();

  for (UInt s = 1; s < nbSteps; ++s) {

    if (s % 10 == 0 && prank == 0) {
      constexpr char wheel[] = "/-\\|";
      auto elapsed = clk::now() - start_time;
      auto time_per_step = elapsed / s;
      std::cout << "\r[" << wheel[(s / 10) % 4] << "] " << std::setw(5) << s
                << "/" << nbSteps << " (" << std::setprecision(2) << std::fixed
                << std::setw(8) << millisecond(time_per_step).count()
                << "ms/step - elapsed: " << std::setw(8)
                << second(elapsed).count() << "s - ETA: " << std::setw(8)
                << second((nbSteps - s) * time_per_step).count() << "s)"
                << std::string(' ', 20) << std::flush;
    }
    model.applyBC(BC::Dirichlet::IncrementValue(increment, _y), "top");

    for (UInt i = 0; i < nb_staggered_steps; ++i) {
      coupler.solve();
    }

    if (s % 100 == 0) {
      model.dump();
    }
  }

  return 0;
}
