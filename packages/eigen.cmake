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


set(AKANTU_EIGEN_VERSION 3.3.7)

package_declare(Eigen3 EXTERNAL NOT_OPTIONAL DEFAULT ON
  DESCRIPTION "Add Eigen3 dependency to akantu"
  SYSTEM AUTO third-party/cmake/eigen3.cmake
  EXTRA_PACKAGE_OPTIONS ARGS ${AKANTU_EIGEN_VERSION} ARGS NO_MODULE TARGET Eigen3::Eigen
  )

mark_as_advanced(Eigen3_DIR)
package_add_third_party_script_variable(Eigen3
  EIGEN3_VERSION "5.0.1")
package_add_third_party_script_variable(Eigen3
  EIGEN3_GIT "https://gitlab.com/libeigen/eigen.git")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU"
    AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "12.0.0")
  package_set_compile_flags(Eigen3 CXX "-Wno-array-bounds -Wno-use-after-free -Wno-stringop-overread")
endif()
