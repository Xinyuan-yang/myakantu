include(FetchContent)
FetchContent_Declare(
  gbenchmark
  SYSTEM
  GIT_REPOSITORY https://github.com/google/benchmark.git
  GIT_TAG        "v1.9.5"
)

set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE INTERNAL "" FORCE)
set(BENCHMARK_ENABLE_TESTING OFF CACHE INTERNAL "" FORCE)

FetchContent_MakeAvailable(gbenchmark)

set_property(TARGET benchmark PROPERTY CXX_STANDARD 14)
add_library(benchmark::benchmark ALIAS benchmark)

set(GBENCHMARK_LIBRARIES benchmark::benchmark CACHE INTERNAL "" FORCE)

mask_package_options(BENCHMARK)
mark_as_advanced(gbenchmark_DIR LIBRT LLVM_FILECHECK_EXE)
