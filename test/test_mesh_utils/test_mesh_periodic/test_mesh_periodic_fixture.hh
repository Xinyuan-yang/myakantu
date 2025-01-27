/**
 * Copyright (©) 2018-2024 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#include "aka_common.hh"
#include "dumper_iohelper_paraview.hh"
#include "mesh_utils.hh"
#include "parser.hh"
#include <gtest/gtest.h>

#ifndef AKANTU_TEST_MESH_PERIODIC_FIXTURE_HH_
#define AKANTU_TEST_MESH_PERIODIC_FIXTURE_HH_

using namespace akantu;

template <class Params> class TestMeshPeriodic : public ::testing::Test {
public:
  static constexpr Int dim{std::tuple_element_t<0, Params>{}};
  static constexpr Int dim_periodic{std::tuple_element_t<1, Params>{}};

  void SetUp() override {
    debug::setDebugLevel(dblError);

    const auto & comm = Communicator::getStaticCommunicator();
    Int psize = comm.getNbProc();
    Int prank = comm.whoAmI();
    AKANTU_DEBUG_INFO("Running on " << psize << " processors");
    AKANTU_DEBUG_INFO("This is processor " << prank);

    mesh = std::make_unique<Mesh>(dim);
    if (prank == 0) {
      if (dim == 1) {
        mesh->read("line.msh");
      } else if (dim == 2) {
        mesh->read("square.msh");
      } else if (dim == 3) {
        mesh->read("cube.msh");
      }
    }
    mesh->distribute();
  }

  void TearDown() override { mesh.reset(nullptr); }

  void checkAll() { EXPECT_TRUE(mesh->isPeriodic()); }

protected:
  void dump(std::string name_suffix) {
    auto dumper = std::make_shared<DumperParaview>("periodic",
                                                   "./paraview-" + name_suffix);
    mesh->registerExternalDumper(dumper, "periodic", true);
    mesh->addDumpMesh(*mesh);

    // Gather info on which nodes are masters and slaves
    Array<Int> periodic_id(mesh->getNbNodes(), 1, 0.);
    Array<Int> masters(mesh->getNbNodes(), 1, 0.);
    Array<Int> slaves(mesh->getNbNodes(), 1, 0.);
    int prev_node = -1;
    UInt value = 0;
    const auto & periodic_ms = mesh->getPeriodicMasterSlaves();
    for (auto && [master, slave] : periodic_ms) {
      if (prev_node != master) {
        ++value;
      }

      prev_node = master;
      periodic_id(master) = value;
      periodic_id(slave) = value;

      masters(master) = 1;
      slaves(slave) = 1;

      auto it = periodic_ms.find(slave);
      if (it != periodic_ms.end()) {
        AKANTU_EXCEPTION(slave << " is slave of " << master << " and master of "
                               << it->second);
      }
    }

    // Get processor id for each node
    Array<Int> processor(mesh->getNbNodes(), 1, 0.);
    const auto & comm = Communicator::getStaticCommunicator();
    Int prank = comm.whoAmI();
    for (auto node : arange(mesh->getNbNodes())) {
      processor(node) = prank;
    }

    mesh->addDumpFieldExternalToDumper("periodic", "periodic_id", periodic_id);
    mesh->addDumpFieldExternalToDumper("periodic", "masters", masters);
    mesh->addDumpFieldExternalToDumper("periodic", "slaves", slaves);
    mesh->addDumpFieldExternalToDumper("periodic", "processor", processor);
    mesh->dump();
  }

protected:
  std::unique_ptr<Mesh> mesh;
};

#endif /* AKANTU_TEST_MESH_PERIODIC_FIXTURE_HH_ */
