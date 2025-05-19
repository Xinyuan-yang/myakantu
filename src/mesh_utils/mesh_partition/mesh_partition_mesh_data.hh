/**
 * Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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

#ifndef AKANTU_MESH_PARTITION_MESH_DATA_HH_
#define AKANTU_MESH_PARTITION_MESH_DATA_HH_

/* -------------------------------------------------------------------------- */
#include "aka_common.hh"
#include "mesh.hh"
#include "mesh_partition.hh"
/* -------------------------------------------------------------------------- */
#include <functional>
#include <memory>
#include <string>
/* -------------------------------------------------------------------------- */

namespace akantu {

class MeshPartitionMeshData : public MeshPartition {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */
public:
  MeshPartitionMeshData(Mesh & mesh, Int spatial_dimension,
                        const ID & id = "MeshPartitionerMeshData");

  MeshPartitionMeshData(Mesh & mesh,
                        std::shared_ptr<ElementTypeMapArray<Idx>> & mapping,
                        Int spatial_dimension,
                        const ID & id = "MeshPartitionerMeshData");

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */
  void partitionate(
      Int nb_part,
      const std::function<Int(const Element &, const Element &)> &
          edge_load_func =
              [](auto && /*unused*/, auto && /*unused*/) { return 1; },
      const std::function<Int(const Element &)> & vertex_load_func =
          [](auto && /*unused*/) { return 1; }) override;

  void reorder() override;

  void setPartitionMapping(std::shared_ptr<ElementTypeMapArray<Idx>> & mapping);

  void setPartitionMappingFromMeshData(const std::string & data_name);

  /* ------------------------------------------------------------------------ */
  /* Class Members                                                            */
  /* ------------------------------------------------------------------------ */
private:
  std::shared_ptr<ElementTypeMapArray<Idx>> partition_mapping;
};

/* -------------------------------------------------------------------------- */
/* inline functions                                                           */
/* -------------------------------------------------------------------------- */

} // namespace akantu

#endif /* AKANTU_MESH_PARTITION_MESH_DATA_HH_ */
