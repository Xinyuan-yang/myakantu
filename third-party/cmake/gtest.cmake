include(FetchContent)
FetchContent_Declare(
  googletest
  SYSTEM
  GIT_REPOSITORY ${GTEST_GIT}
  GIT_TAG        ${GTEST_VERSION}
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(gtest_build_samples OFF CACHE INTERNAL "" FORCE)
set(BUILD_GTEST ON CACHE INTERNAL "" FORCE)
set(BUILD_GMOCK OFF CACHE INTERNAL "" FORCE)

FetchContent_MakeAvailable(googletest)

set(GTEST_LIBRARIES GTest::gtest_main CACHE INTERNAL "" FORCE)
