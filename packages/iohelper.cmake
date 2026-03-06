#===============================================================================
# Copyright (©) 2011-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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
if(AKANTU_BYPASS_AKANTU_TARGET)
  return()
endif()

package_declare(IOHelper EXTERNAL NOT_OPTIONAL DEFAULT ON
  DESCRIPTION "Add IOHelper support in akantu"
  SYSTEM OFF third-party/cmake/iohelper.cmake
  EXTRA_PACKAGE_OPTIONS TARGET iohelper)

set(_version "1.1.1")
package_add_third_party_script_variable(IOHelper
  IOHELPER_VERSION ${_version})
package_add_third_party_script_variable(IOHelper
  IOHELPER_GIT "https://c4science.ch/source/iohelper.git")
package_add_third_party_script_variable(IOHelper
  IOHELPER_ARCHIVE "iohelper_${_version}.tar.gz")

package_declare_extra_files_to_package(dumpers
  PROJECT
  third-party/cmake/iohelper.cmake
  )
