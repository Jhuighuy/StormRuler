# Copyright (C) 2020-2023 Oleg Butakov
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.

# ------------------------------------------------------------------------------

# Benchmarks executable.
file(GLOB CXX_SOURCES "./*.cpp")
add_executable(StormRuler_Benchmarks ${CXX_SOURCES})

# ------------------------------------------------------------------------------

# Link with the core library.
target_link_libraries(StormRuler_Benchmarks StormRuler)

# ------------------------------------------------------------------------------

# Find doctest.
find_package(doctest CONFIG REQUIRED)
target_link_libraries(StormRuler_Benchmarks doctest::doctest)

# Find nanobench.
# TODO: Something is broken with nanobench + vcpkg + CMake.
#find_package(nanobench CONFIG REQUIRED)
#target_link_libraries(StormRuler_Benchmarks nanobench)

# ------------------------------------------------------------------------------

# Find Armadillo.
find_package(Armadillo CONFIG QUIET)
if(Armadillo_FOUND)
  set(STORM_BENCH_ARMADILLO_ENABLED 1)
  message(STATUS "| Benchmarks: Armadillo enabled.")
  target_link_libraries(StormRuler_Benchmarks ${ARMADILLO_LIBRARIES})
endif()

# Find Blaze.
find_package(blaze CONFIG QUIET)
if(blaze_FOUND)
  set(STORM_BENCH_BLAZE_ENABLED 1)
  message(STATUS "| Benchmarks: Blaze enabled.")
  target_link_libraries(StormRuler_Benchmarks blaze::blaze)
endif()

# Find Eigen.
find_package(Eigen3 3.4 NO_MODULE QUIET)
if(Eigen3_FOUND)
  set(STORM_BENCH_EIGEN_ENABLED 1)
  message(STATUS "| Benchmarks: Eigen enabled.")
  target_link_libraries(StormRuler_Benchmarks Eigen3::Eigen)
endif()

# Find XTensor.
find_package(xtensor QUIET)
if(xtensor_FOUND)
  set(STORM_BENCH_XTENSOR_ENABLED 1)
  message(STATUS "| Benchmarks: XTensor enabled.")
  target_link_libraries(StormRuler_Benchmarks xtensor)
endif()

# ------------------------------------------------------------------------------

# Configuration file.
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/_BenchmarksConfig.hpp.in"
    "${CMAKE_CURRENT_SOURCE_DIR}/_BenchmarksConfig.hpp")

# ------------------------------------------------------------------------------
