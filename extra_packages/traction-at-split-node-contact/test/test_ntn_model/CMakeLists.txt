#===============================================================================
# Copyright (©) 2018-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
#===============================================================================

register_gtest_sources(
  SOURCES test_ntn_basics.cc
  FILES_TO_COPY
  ntn_test_2d.msh
  ntn_test_3d.msh
  data/material_ntn_model.dat
  PACKAGE core)

register_gtest_test(test_ntn_model)

register_gtest_sources(
  SOURCES test_ntn_friction.cc
  FILES_TO_COPY
  ntn_test_compression.msh
  mesh-H6.0-L1.0.msh
  data/material_ntn_friction_slip_weakening.dat
  data/material_ntn_friction_rate_and_state.dat
  PACKAGE core)

register_gtest_test(test_ntn_friction)
