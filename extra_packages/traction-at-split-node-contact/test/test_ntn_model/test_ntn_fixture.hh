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
#include "aka_common.hh"
#include "dumpable.hh"
#include "mesh.hh"
#include "mesh_partition_mesh_data.hh"
#include "mesh_utils.hh"
#include "ntn_base_contact.hh"
#include "ntn_contact.hh"
#include "solid_mechanics_model.hh"

/* -------------------------------------------------------------------------- */
#include <cmath>
#include <cstdlib>
#include <gtest/gtest.h>
#include <map>
#include <memory>
#include <string>
#include <utility>

#ifndef AKANTU_NTN_TEST_FIXTURE_HH_
#define AKANTU_NTN_TEST_FIXTURE_HH_

namespace akantu {

template <typename dim_> class TestNTNFixture : public ::testing::Test {
public:
  static constexpr Int dim = dim_::value;

  void SetUp() override {
    this->mesh = std::make_unique<Mesh>(dim);

    const auto & comm = Communicator::getStaticCommunicator();
    auto psize = comm.getNbProc();
    auto prank = comm.whoAmI();

    std::shared_ptr<MeshPartition> partitioner;
    if (prank == 0) {
      this->mesh->read("ntn_test_" + std::to_string(dim) + "d.msh");

      auto length =
          (mesh->getUpperBounds()(_x) - mesh->getLowerBounds()(_x)) / psize;

      auto partition = std::make_shared<ElementTypeMapArray<Int>>("partition");
      partition->initialize(*mesh, _spatial_dimension = dim,
                            _with_nb_element = true);

      std::set<Idx> sanity_proc_check;

      for_each_element(
          *mesh,
          [&](auto && element) {
            auto barycenter = mesh->getBarycenter(element);
            auto proc = Idx(std::floor(barycenter(_x) / length));

            sanity_proc_check.insert(proc);
            (*partition)(element) = proc;
          },
          _spatial_dimension = dim);

      AKANTU_DEBUG_ASSERT(sanity_proc_check.size() == std::size_t(psize),
                          "There are too many procs for this mesh");

      partitioner =
          std::make_shared<MeshPartitionMeshData>(*mesh, partition, dim);

      partitioner->partitionate(psize);
    }

    this->mesh->distribute(partitioner);

    std::string filename = "material_ntn_model.dat";
    model = std::make_unique<SolidMechanicsModel>(*mesh);

    model->addDumpField("partitions");
    model->dump();

    ntn_contact = std::make_shared<NTNContact>(*model);
    getStaticParser().parse(filename);
    model->initFull(SolidMechanicsModelOptions(_explicit_lumped_mass));
    model->assembleMassLumped();

    ntn_contact->addSurfacePair("slider_bottom", "base_top", 1);
    contact = ntn_contact;

    contact->updateNormals();
    contact->updateLumpedBoundary();
  }

  void TearDown() override {
    model.reset(nullptr);
    mesh.reset(nullptr);
  }

  virtual void initDynamics(Real time_step_factor) {
    // Set the time step and compute the impedance
    Real stable_time_step = model->getStableTimeStep();
    Real time_step = stable_time_step * time_step_factor;
    model->setTimeStep(time_step);
    contact->updateImpedance();
  }

  virtual void checkNodesPairing() {
    auto & slaves = this->ntn_contact->getSlaves().getArray();
    auto & masters = this->ntn_contact->getMasters().getArray();
    const auto & position = this->mesh->getNodes();

    for (auto && [s, m] : zip(slaves, masters)) {
      for (Int d = 0; d < dim; d++) {
        auto distance = position(s, d) - position(m, d);
        EXPECT_NEAR(0, distance, 1e-10);
      }
    }
  }

  virtual void checkLumpedBoundary() {

    const auto & l_boundary_slaves =
        this->ntn_contact->getLumpedBoundarySlaves();
    const auto & l_boundary_masters =
        this->ntn_contact->getLumpedBoundaryMasters();

    const auto & slaves = this->ntn_contact->getSlaves().getArray();
    const auto & masters = this->ntn_contact->getMasters().getArray();
    const auto & position = this->mesh->getNodes();

    const Vector<Real> lowerBounds = this->mesh->getLowerBounds();
    const Vector<Real> upperBounds = this->mesh->getUpperBounds();
    Int nb_contact_nodes = this->ntn_contact->getNbContactNodes();

    auto x_min = lowerBounds(_x);
    auto x_max = upperBounds(_x);
    auto size_x = x_max - x_min;

    Real tol = 1e-11;

    auto args = zip(filter(slaves, make_view(position, dim)),
                    filter(masters, make_view(position, dim)),
                    l_boundary_slaves, l_boundary_masters);

    auto is_bound = [tol](auto && pos, auto && dir, auto && min, auto && max) {
      return (std::abs(pos(dir) - min) < tol) or
             (std::abs(pos(dir) - max) < tol);
    };

    switch (dim) {
    // 2D
    case 2: {
      auto l_b_bulk = size_x / (nb_contact_nodes - 1);
      auto l_b_side = l_b_bulk / 2;

      std::map<bool, Real> boun = {{true, l_b_side}, {false, l_b_bulk}};

      for (auto && [sp, mp, ls, lm] : args) {
        // Slaves
        auto l_b_s = boun[is_bound(sp, _x, x_min, x_max)];
        EXPECT_NEAR(l_b_s, ls, 1e-11);

        // Master
        auto l_b_m = boun[is_bound(mp, _x, x_min, x_max)];
        EXPECT_NEAR(l_b_m, lm, 1e-11);
      }
      break;
    }
    // 3D
    case 3: {
      // Same number of nodes in both directions
      Int nb_contact_per_sides = std::sqrt(nb_contact_nodes);

      auto z_min = lowerBounds(_z);
      auto z_max = upperBounds(_z);
      auto size_z = z_max - z_min;

      auto l_b_bulk = size_x * size_z /
                      ((nb_contact_per_sides - 1) * (nb_contact_per_sides - 1));
      auto l_b_side = l_b_bulk / 2;
      auto l_b_corner = l_b_bulk / 4;

      std::map<std::pair<bool, bool>, Real> boun = {{{true, true}, l_b_corner},
                                                    {{true, false}, l_b_side},
                                                    {{false, true}, l_b_side},
                                                    {{false, false}, l_b_bulk}};

      for (auto && [sp, mp, ls, lm] : args) {
        // Slaves
        auto l_b_s = boun[{is_bound(sp, _x, x_min, x_max),
                           is_bound(sp, _z, z_min, z_max)}];
        EXPECT_NEAR(l_b_s, ls, 1e-11);

        // Masters
        auto l_b_m = boun[{is_bound(mp, _x, x_min, x_max),
                           is_bound(mp, _z, z_min, z_max)}];
        EXPECT_NEAR(l_b_m, lm, 1e-11);
      }

      break;
    }
    }
  }

  virtual void checkImpedance() {
    // Expected impedance is 1/ [dt/2 *(lumped_boundary[slaves]/mass[slaves] +
    // lumped_boundary[masters]/mass[masters])]
    const auto & impedance = this->ntn_contact->getImpedance();
    Real delta_t = this->model->getTimeStep();

    const Vector<Real> lowerBounds = this->mesh->getLowerBounds();
    const Vector<Real> upperBounds = this->mesh->getUpperBounds();
    Real bottom = lowerBounds(1);
    Real top = upperBounds(1);
    Real system_height = top - bottom;

    Real density = 1;

    // Expected impedance can be reduced to
    // density * system_height / (4*dt)
    Real imp_ref = density * system_height / (4 * delta_t);

    for (auto && imp : impedance) {
      EXPECT_NEAR(imp, imp_ref, 1e-11);
    }
  }

protected:
  std::unique_ptr<Mesh> mesh;
  std::unique_ptr<SolidMechanicsModel> model;
  std::shared_ptr<NTNBaseContact> contact;
  std::shared_ptr<NTNContact> ntn_contact;

}; // namespace akantu

using dim_types = testing::Types<dim_2_t, dim_3_t>;
// using dim_types = gtest_list_t<TestAllDimensions>;

TYPED_TEST_SUITE(TestNTNFixture, dim_types, );

} // namespace akantu

#endif /* AKANTU_NTN_TEST_FIXTURE_HH_ */
