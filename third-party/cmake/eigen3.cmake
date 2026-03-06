include(FetchContent)
FetchContent_Declare(
  eigen3
  SYSTEM
  GIT_REPOSITORY ${EIGEN3_GIT}
  GIT_TAG        ${EIGEN3_VERSION}
)

set(CMAKE_BUILD_TYPE Release)
set(BUILD_TESTING OFF CACHE BOOL "Eigen Tests" FORCE)

FetchContent_MakeAvailable(eigen3)

set(EIGEN3_LIBRARIES Eigen3::Eigen)
package_set_libraries(Eigen3 ${EIGEN3_LIBRARIES})
package_add_to_export_list(Eigen3 eigen)

install(TARGETS eigen EXPORT ${AKANTU_TARGETS_EXPORT})


mask_package_options(EIGEN3)
mask_package_options(EIGEN)
mask_package_options(PASTIX)
mark_as_advanced(
  BUILD_TESTING
  INCLUDE_INSTALL_DIR
  PKGCONFIG_INSTALL_DIR
  PATCH_COMMAND
  QT_QMAKE_EXECUTABLE
  RT_LIBRARY
  CMAKEPACKAGE_INSTALL_DIR
  BLA_VENDOR
  )
