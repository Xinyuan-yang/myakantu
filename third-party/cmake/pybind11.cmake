set(CMAKE_CXX_STANDARD ${AKANTU_CXX_STANDARD}) # Otherwhy pybind11 default to cxx14

include(FetchContent)
FetchContent_Declare(
  pybind11
  SYSTEM
  GIT_REPOSITORY    ${PYBIND11_GIT}
  GIT_TAG           ${PYBIND11_VERSION}
)

set(PYBIND11_PYTHON_VERSION ${AKANTU_PREFERRED_PYTHON_VERSION})
set(PYBIND11_FINDPYTHON NEW)
FetchContent_MakeAvailable(pybind11)

mask_package_options(PYBIND11)
mask_package_options(Python)
