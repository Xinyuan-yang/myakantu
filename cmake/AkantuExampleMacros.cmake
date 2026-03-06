#===============================================================================
# Copyright (©) 2016-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


include(AkantuSimulationMacros)

set(AKANTU_DRIVER_SCRIPT ${AKANTU_CMAKE_DIR}/akantu_test_driver.sh)

# ==============================================================================
function(register_example example_name)
  _add_akantu_simulation(${example_name} ${ARGN} LIST_FILES _example_files)
  if(DEFINED _add_examples_pkg)
    package_add_to_variable(EXAMPLES_FILES ${_add_examples_pkg} ${_example_files})
  endif()

  if(AKANTU_TEST_EXAMPLES)
    cmake_parse_arguments(_example
      "PYTHON;PARALLEL"
      "FILES_TO_COPY;TEST_DEPENDS"
      "SCRIPT"
      ${ARGN}
      )

    if(_example_FILES_TO_COPY)
      foreach(_file ${_example_FILES_TO_COPY})
        _add_file_to_copy(${example_name} "${_file}")
      endforeach()
    endif()

    if(_example_PARALLEL)
      set(_exe ${MPIEXEC})
      if(NOT _exe)
        set(_exe ${MPIEXEC_EXECUTABLE})
      endif()
      set(_parallel_runner -p "${_exe} ${MPIEXEC_PREFLAGS} ${MPIEXEC_NUMPROC_FLAG}" -N "2")
    endif()

    set(_test_driver
      ${AKANTU_DRIVER_SCRIPT}
      ${_parallel_runner}
      -n "${example_name}-test"
      -w "${CMAKE_CURRENT_BINARY_DIR}"
      -E "${PROJECT_BINARY_DIR}/akantu_environement.sh"
    )

    if(NOT _example_SCRIPT)
      set(_executable $<TARGET_FILE:${example_name}>)
    elseif(_example_SCRIPT)
      _add_file_to_copy(${example_name} "${_example_SCRIPT}")
      if(_example_PYTHON)
        if(NOT Python_EXECUTABLE)
          find_package(Python
            ${AKANTU_PREFERRED_PYTHON_VERSION}
            COMPONENTS Interpreter
            REQUIRED)
        endif()
        list(APPEND _test_driver -P "${Python_EXECUTABLE}")
      endif()
      set(_executable "${_example_SCRIPT}")
    endif()

    add_test(
      NAME ${example_name}-test
      COMMAND ${_test_driver} -e ${_executable})

    if(_example_TEST_DEPENDS)
      set(_list_deps)
      foreach(dep ${_example_TEST_DEPENDS})
         list(APPEND _list_deps "${dep}-test")
      endforeach()
      set_tests_properties(${example_name}-test PROPERTIES DEPENDS "${_list_deps}")
    endif()
  endif()
endfunction()

# ==============================================================================
function(add_example et_name desc)
  string(TOUPPER ${et_name} upper_name)

  if(NOT _build_all_ex)
    option(AKANTU_BUILD_ALL_EXAMPLES "Activate all examples" ON)
    set( _build_all_ex TRUE)
  endif()

  option(AKANTU_BUILD_EXAMPLES_${upper_name} ${desc} OFF)

  if(AKANTU_BUILD_ALL_EXAMPLES)
    mark_as_advanced(FORCE AKANTU_BUILD_EXAMPLES_${upper_name})
  else()
    mark_as_advanced(CLEAR AKANTU_BUILD_EXAMPLES_${upper_name})
  endif()

  if(AKANTU_BUILD_EXAMPLES_${upper_name} OR AKANTU_BUILD_ALL_EXAMPLES)

    if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${et_name})
      message(FATAL_ERROR "The folder ${CMAKE_CURRENT_SOURCE_DIR}/${et_name} "
        "that you try to register as an example sub-folder, does not exists.")
    endif()

    cmake_parse_arguments(_manage_example
      ""
      ""
      "PACKAGE"
      ${ARGN}
      )

    if(_manage_example_PACKAGE)
      set(_act TRUE)
      foreach(_pkg ${_manage_example_PACKAGE})
        package_is_activated(${_pkg} _activated)
        if(NOT _activated)
          set(_act FALSE)
        endif()
      endforeach()
    else()
      message(SEND_ERROR "Examples should be associated to a package")
    endif()

    if(_act)
      if(DEFINED _add_examples_pkg)
        set(_save_add_examples_pkg ${_add_examples_pkg})
      endif()
      list(GET _manage_example_PACKAGE 0 _pkg)
      set(_add_examples_pkg ${_pkg})

      add_subdirectory(${et_name})

      unset(_add_examples_pkg)
      if(DEFINED _save_add_examples_pkg)
        set(_add_examples_pkg ${_save_add_examples_pkg})
        unset(_save_add_examples_pkg)
      endif()
    endif()
  endif()
endfunction()
