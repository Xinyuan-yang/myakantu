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


#===============================================================================
# List of tests
#===============================================================================
add_akantu_test(test_mesh_io "Test mesh io object")
add_akantu_test(test_pbc_tweak "Test pbc facilities")
add_akantu_test(test_mesh_periodic "Test periodic mesh")
add_akantu_test(test_buildfacets "Tests for the generation of facets")
add_akantu_test(test_segment_nodetype "segment_nodetype")

register_test(test_purify_mesh
  SOURCES test_purify_mesh.cc
  FILES_TO_COPY purify_mesh.msh
  PACKAGE core
  )

add_mesh(test_mesh_iterators_mesh iterators_mesh.geo 3 1 OUTPUT iterators_mesh.msh)
register_test(test_mesh_iterators
  SOURCES test_mesh_iterators.cc
  PACKAGE core
  DEPENDS test_mesh_iterators_mesh
  )

add_akantu_test(test_mesh_partitionate "Test mesh partition creation")
