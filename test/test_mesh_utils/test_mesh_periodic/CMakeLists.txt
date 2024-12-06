#===============================================================================
# Copyright (©) 2010-2024 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# This file is part of Akantu
# 
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
# 
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

add_mesh(square_mesh square.geo DIM 2 ORDER 1 OUTPUT square.msh)
add_mesh(cube_mesh cube.geo DIM 3 ORDER 1 OUTPUT cube.msh)

register_gtest_sources(
  SOURCES test_mesh_periodic.cc
  DEPENDS square_mesh cube_mesh
  FILES_TO_COPY material.dat
  PACKAGE core
  PARALLEL
  PARALLEL_LEVEL 1 2
  )

register_gtest_sources(
  SOURCES test_periodic_static.cc
  DEPENDS cube_mesh
  FILES_TO_COPY material.dat
  PACKAGE core
  PARALLEL
  PARALLEL_LEVEL 1 2
  )

register_gtest_sources(test_periodic_mass
  SOURCES test_periodic_mass.cc
  DEPENDS cube_mesh
  FILES_TO_COPY material.dat
  PACKAGE core
  PARALLEL
  PARALLEL_LEVEL 1 2
  )

register_gtest_test(test_mesh_periodic)
