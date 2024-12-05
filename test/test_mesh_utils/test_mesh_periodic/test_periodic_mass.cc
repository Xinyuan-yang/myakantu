/**
 * Copyright (©) 2010-2024 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_utils.hh"
#include "solid_mechanics_model.hh"
#include "dumper_iohelper_paraview.hh"
#include <gtest/gtest.h>

using namespace akantu;

TEST(TestPeriodicMass, UniformMass) {
  const size_t dim = 3;

  debug::setDebugLevel(dblError);
  getStaticParser().parse("material.dat");

  const Communicator & comm = Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();
  AKANTU_DEBUG_INFO("Running on " << psize << " processors");
  AKANTU_DEBUG_INFO("This is processor " << prank);

  // Create periodic mesh
  auto mesh = std::make_unique<Mesh>(dim);
  if (prank == 0) {
    mesh->read("cube.msh");
  }
  mesh->distribute();

  mesh->makePeriodic(_x);
  mesh->makePeriodic(_y);
  mesh->makePeriodic(_z);

  // Initialize model with lumped mass
  auto model = std::make_unique<SolidMechanicsModel>(*mesh);
  model->initFull(_analysis_method = _explicit_lumped_mass);
  model->assembleMassLumped();

  // Check that all the masses are the same
  const auto & masses = model->getMass();
  auto reference_mass = masses[0];
  for (auto & masses_i : make_view(masses, dim)) {
    for (auto & mass : masses_i) {
      EXPECT_NEAR(0., (mass - reference_mass) / reference_mass, 2e-11);
    }
  }

  // Dump mass
  auto dumper = std::make_shared<DumperParaview>("periodic-mass", "./paraview");
  mesh->registerExternalDumper(dumper, "periodic-mass", true);
  mesh->addDumpMesh(*mesh);

  // Gather info on which nodes are masters and slaves
  mesh->addDumpFieldExternalToDumper("periodic-mass", "mass", masses);
  mesh->dump();
}
