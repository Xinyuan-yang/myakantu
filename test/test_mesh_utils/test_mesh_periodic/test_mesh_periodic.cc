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
#include <cassert>

using namespace akantu;

int main(int argc, char * argv[]) {
  int dim = 3;

  initialize("material.dat", argc, argv);
  debug::setDebugLevel(akantu::dblInfo);

  const Communicator & comm = Communicator::getStaticCommunicator();
  Int psize = comm.getNbProc();
  Int prank = comm.whoAmI();
  AKANTU_DEBUG_INFO("Running on " << psize << " processors");
  AKANTU_DEBUG_INFO("This is processor " << prank);

  Mesh mesh(dim);
  if (prank == 0) {
    mesh.read("cube.msh");
  }
  mesh.distribute();
  mesh.makePeriodic(_x);
  assert(mesh.isPeriodic());

  auto dumper = std::make_shared<DumperParaview>("periodic", "./paraview");
  mesh.registerExternalDumper(dumper, "periodic", true);
  mesh.addDumpMesh(mesh);

  // Gather info on which nodes are masters and slaves
  Array<Int> periodic_id(mesh.getNbNodes(), 1, 0.);
  Array<Int> masters(mesh.getNbNodes(), 1, 0.);
  Array<Int> slaves(mesh.getNbNodes(), 1, 0.);
  int prev_node = -1;
  UInt value = 0;
  const auto & periodic_ms = mesh.getPeriodicMasterSlaves();
  for (auto & pair : periodic_ms) {
    if (prev_node != pair.first) {
      ++value;
    }

    prev_node = pair.first;
    periodic_id(pair.first) = value;
    periodic_id(pair.second) = value;

    masters(pair.first) = 1;
    slaves(pair.second) = 1;

    auto it = periodic_ms.find(pair.second);
    if (it != periodic_ms.end()) {
      AKANTU_EXCEPTION(pair.second << " is slave of " << pair.first
                                   << " and master of " << it->second);
    }
  }

  // Get processor id for each node
  Array<Int> processor(mesh.getNbNodes(), 1, 0.);
  for (auto node : arange(mesh.getNbNodes())) {
    processor(node) = prank;
  }

  mesh.addDumpFieldExternalToDumper("periodic", "periodic_id", periodic_id);
  mesh.addDumpFieldExternalToDumper("periodic", "masters", masters);
  mesh.addDumpFieldExternalToDumper("periodic", "slaves", slaves);
  mesh.addDumpFieldExternalToDumper("periodic", "processor", processor);
  mesh.dump();

  return EXIT_SUCCESS;
}
