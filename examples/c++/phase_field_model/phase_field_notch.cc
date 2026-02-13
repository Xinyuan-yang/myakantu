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
#include "coupler_solid_phasefield.hh"
#include "phase_field_element_filter.hh"
/* -------------------------------------------------------------------------- */
#include <chrono>
#include <iostream>
/* -------------------------------------------------------------------------- */
using namespace akantu;
/* -------------------------------------------------------------------------- */

using clk = std::chrono::high_resolution_clock;
using second = std::chrono::duration<double>;
using millisecond = std::chrono::duration<double, std::milli>;

const Int spatial_dimension = 2;

/* -------------------------------------------------------------------------- */
int main(int argc, char * argv[]) {

  initialize("material_notch.dat", argc, argv);

  // create mesh
  Mesh mesh(spatial_dimension);
  mesh.read("square_notch.msh");

  // The model coupler contains the solid mechanics and the phase-field models
  CouplerSolidPhaseField coupler(mesh);
  auto & model = coupler.getSolidMechanicsModel();
  auto & phase = coupler.getPhaseFieldModel();

  // Each model can bet set separately
  model.initFull(_analysis_method = _static);
  phase.initFull(_analysis_method = _static);

  // Dirichlet BC
  model.applyBC(BC::Dirichlet::FixedValue(0., _y), "bottom");
  model.applyBC(BC::Dirichlet::FixedValue(0., _x), "left");

  // Dumper settings
  model.setBaseName("phase_notch");
  model.addDumpField("stress");
  model.addDumpField("grad_u");
  model.addDumpFieldVector("displacement");
  model.addDumpField("damage");
  model.dump();

  const Int nb_steps = 1000;
  Real increment = 4e-6;
  Int nb_staggered_steps = 10;

  auto start_time = clk::now();

  // Main loop over the loading steps
  for (Int s = 1; s < nb_steps; ++s) {

    if (s % 200 == 0) {
      constexpr std::array<char, 5> wheel{"/-\\|"};
      auto elapsed = clk::now() - start_time;
      auto time_per_step = elapsed / s;
      int idx = (s / 10) % 4;
      std::cout << "\r[" << wheel.at(idx) << "] " << std::setw(5) << s << "/"
                << nb_steps << " (" << std::setprecision(2) << std::fixed
                << std::setw(8) << millisecond(time_per_step).count()
                << "ms/step - elapsed: " << std::setw(8)
                << second(elapsed).count() << "s - ETA: " << std::setw(8)
                << second((nb_steps - s) * time_per_step).count() << "s)"
                << std::string(' ', 20) << std::flush;
    }
    model.applyBC(BC::Dirichlet::IncrementValue(increment, _y), "top");

    /* At each step, the two solvers are called alternately. Here a fixed number
     * of staggered iterations is set for simplicity but a convergence based on
     * the displacements, damage and/or energy can also be used to exit the
     * internal loop.*/
    for (Idx i = 0; i < nb_staggered_steps; ++i) {
      coupler.solve();
    }

    if (s % 100 == 0) {
      model.dump();
    }
  }

  /* Here a damage threshold is set and a mesh clustering function is called
   * using the phase-field element filter. This allows to cluster the mesh into
   * groups of elements separated by the mesh boundaries and "broken" elements
   * with damage above the threshold. */
  Real damage_limit = 0.15;
  auto global_nb_clusters = mesh.createClusters(
      spatial_dimension, "crack", PhaseFieldElementFilter(phase, damage_limit));

  model.dumpGroup("crack_0");
  model.dumpGroup("crack_1");

  std::cout << std::endl;
  std::cout << "Nb clusters: " << global_nb_clusters << std::endl;

  finalize();
  return EXIT_SUCCESS;
}
