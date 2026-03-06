#===============================================================================
# @file   iohelper.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Mar 30 2015
# @date last modification: Tue Jan 19 2016
#
# @brief  build script for iohelper
#
# @section LICENSE
#
# Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
# (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================
set(IOHELPER_TARGETS_EXPORT ${AKANTU_TARGETS_EXPORT})
add_subdirectory(${PROJECT_SOURCE_DIR}/third-party/iohelper)

set(IOHELPER_SOURCE_TYPE "internal"
  CACHE INTERNAL "internal variable for clean export")

#set(IOHELPER_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/third-party/iohelper/src" CACHE INTERNAL "IOHelper include directory")

set(IOHELPER_LIBRARIES iohelper CACHE INTERNAL "IOHelper library")
package_add_to_export_list(iohelper iohelper)
mask_package_options(IOHELPER)
