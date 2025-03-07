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
#include "communicator.hh"
#include "dumper_iohelper_paraview.hh"
#include "mesh.hh"
#include "mesh_utils.hh"
#include "parser.hh"
/* -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
/* -------------------------------------------------------------------------- */

#ifndef AKANTU_TEST_MESH_PERIODIC_FIXTURE_HH_
#define AKANTU_TEST_MESH_PERIODIC_FIXTURE_HH_

using namespace akantu;

struct dim0 : public std::integral_constant<Int, 0> {};
struct dim1 : public dim_1_t {};
struct dim2 : public dim_2_t {};
struct dim3 : public dim_3_t {};

enum class LoadType { no_load, tensile, shear };
struct no_load : public std::integral_constant<LoadType, LoadType::no_load> {};
struct tensile : public std::integral_constant<LoadType, LoadType::tensile> {};
struct shear : public std::integral_constant<LoadType, LoadType::shear> {};

template <class Params> class TestMeshPeriodic : public ::testing::Test {
public:
  static constexpr Int dim{std::tuple_element_t<0, Params>{}};
  static constexpr Int dim_periodic{std::tuple_element_t<1, Params>{}};
  static constexpr LoadType load_type{std::tuple_element_t<2, Params>{}};

  void SetUp() override {
    const auto & comm = Communicator::getStaticCommunicator();
    Int psize = comm.getNbProc();
    Int prank = comm.whoAmI();
    AKANTU_DEBUG_INFO("Running on " << psize << " processors");
    AKANTU_DEBUG_INFO("This is processor " << prank);

    mesh = std::make_unique<Mesh>(dim);
    std::vector<ID> groups;
    if (prank == 0) {
      if (dim == 1) {
        mesh->read("line.msh");
      } else if (dim == 2) {
        mesh->read("square.msh");
      } else if (dim == 3) {
        mesh->read("cube.msh");
      }

      if (dim_periodic == 1) {
        groups = {"x-minus", "x-plus"};
      } else if (dim_periodic == 2) {
        groups = {"x-minus", "x-plus", "y-minus", "y-plus"};
      } else if (dim_periodic == 3) {
        groups = {"x-minus", "x-plus",  "y-minus",
                  "y-plus",  "z-minus", "z-plus"};
      }

      auto && shared = this->mesh->getNodalData<Real>("shared");
      shared.resize(this->mesh->getNbNodes());
      shared.set(0.);

      for (auto group : groups) {
        for (auto node :
             this->mesh->getElementGroup(group).getNodeGroup().getNodes()) {
          shared(node) += 1.;
        }
      }
    }

    mesh->distribute();

    accessor = std::make_unique<MeshAccessor>(*mesh);

    periodic_string = "";
    for (auto sp : spatial_directions<dim_periodic>) {
      this->mesh->makePeriodic(sp);
      periodic_string += std::to_string(sp);
    }

    if (periodic_string.empty()) {
      periodic_string = "none";
    }
  }

  void TearDown() override { mesh.reset(nullptr); }

  void checkAll() { EXPECT_TRUE(mesh->isPeriodic()); }

protected:
  void dump(const std::string & name_suffix) {
    dumper = std::make_shared<DumperParaview>(
        "periodic", "./paraview-" + name_suffix + std::to_string(this->dim) +
                        "d-" + this->periodic_string);
    mesh->registerExternalDumper(dumper, "periodic", true);
    mesh->addDumpMesh(*mesh);

    addDumpFields("periodic");
  }

  void addDumpFields(const std::string & dumper) {
    // Gather info on which nodes are masters and slaves
    fields["periodic_id"] =
        std::make_unique<Array<Int>>(mesh->getNbNodes(), 1, 0);
    fields["masters"] = std::make_unique<Array<Int>>(mesh->getNbNodes(), 1, -1);
    fields["slaves"] = std::make_unique<Array<Int>>(mesh->getNbNodes(), 1, -1);

    auto && periodic_id = *fields["periodic_id"];
    auto && masters = *fields["masters"];
    auto && slaves = *fields["slaves"];

    Int prev_node = -1;
    Int value = 0;
    const auto & periodic_ms = mesh->getPeriodicMasterSlaves();
    for (auto && [master, slave] : periodic_ms) {
      if (prev_node != master) {
        ++value;
      }

      prev_node = master;
      periodic_id(master) = int(mesh->getNodeFlag(master)) >> 4;
      periodic_id(slave) = int(mesh->getNodeFlag(slave)) >> 4;

      masters(master) = mesh->getNodeGlobalId(master);
      slaves(slave) = mesh->getNodeGlobalId(master);

      auto it = periodic_ms.find(slave);
      if (it != periodic_ms.end()) {
        AKANTU_EXCEPTION(slave << " is slave of " << master << " and master of "
                               << it->second);
      }
    }

    // Get processor id for each node
    fields["processor"] =
        std::make_unique<Array<Int>>(mesh->getNbNodes(), 1, 0.);
    auto && processor = *fields["processor"];
    const auto & comm = Communicator::getStaticCommunicator();
    Int prank = comm.whoAmI();
    for (auto node : arange(mesh->getNbNodes())) {
      processor(node) = prank;
    }

    mesh->addDumpFieldExternalToDumper(dumper, "nodes_flags",
                                       this->accessor->getNodesFlags());

    mesh->addDumpFieldExternalToDumper(dumper, "periodic_id", periodic_id);
    mesh->addDumpFieldExternalToDumper(dumper, "masters", masters);
    mesh->addDumpFieldExternalToDumper(dumper, "slaves", slaves);
    mesh->addDumpFieldExternalToDumper(dumper, "processor", processor);
  }

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<MeshAccessor> accessor;
  std::shared_ptr<DumperParaview> dumper;

  std::map<ID, std::unique_ptr<Array<Int>>> fields;

  std::string periodic_string;
};

#endif /* AKANTU_TEST_MESH_PERIODIC_FIXTURE_HH_ */
