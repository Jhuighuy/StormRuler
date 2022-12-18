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
# FITNESS FOR Allocator PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
# SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# TODO:
# • clang-format integration:
#   Format the code on build.
# • clang-tidy integration:
#   Impossible now, since clang-tidy doesn't support C++23.
#   (But now we support C++20.. again)
# • User friendly compilation script:
#   It is inconvenient to write `CXX=.. cmake . -DCMAKE_BUILD_TYPE=..`
#   every time.
# • Add sanitizers for the Coverage build.

# ------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.20..3.26)

message(STATUS "\n" [=[
-- ╔═══════════════════════════════════════════════════════════════════╗
-- ║ ╔═══════════════════════════════════════════════════════════════╗ ║
-- ║ ║      _____ __                       ____        __            ║ ║
-- ║ ║     / ___// /_____   ____ ___ ___  / __ \__  __/ ___  _____   ║ ║
-- ║ ║     \__ \/ __/`__ \/ ___/`__ `__ \/ /_/ / / / / / _ \/ ___/   ║ ║
-- ║ ║    ___/ / /_/ /_/ / /  / / / / / / _, _/ /_/ / /  __/ /       ║ ║
-- ║ ║   /____/\__/\____/_/  /_/ /_/ /_/_/ |_|\__,_/_/\___/_/        ║ ║
-- ║ ║                                                               ║ ║
-- ║ ╚═══════════════════════════════════════════════════════════════╝ ║
-- ╚═══════════════════════════════════════════════════════════════════╝
--]=])

# ------------------------------------------------------------------------------

# Setup the project.
project(StormRuler
        VERSION 0.1
        DESCRIPTION "A very high order FVM framework"
        LANGUAGES CXX)

# A container library that stores all the properties.
add_library(StormRuler_BASE INTERFACE)

# ------------------------------------------------------------------------------

# Set output directories.
# (Set globally, since these properties are not inherited by targets.)
set(CMAKE_PDB_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

# ------------------------------------------------------------------------------

# Check C/C++ compiler.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(Storm_CXX_COMPILER_GNU_LIKE YES)
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "12.1.0")
    message(FATAL_ERROR "GCC version 12.1 or more recent is required!")
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(Storm_CXX_COMPILER_CLANG_LIKE YES)
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "16.0.0")
    message(FATAL_ERROR "Clang version 16.0 or more recent is required!")
  endif()
else()
  message(FATAL_ERROR "Unsupported C++ compiler!")
endif()

# Require C++23 and disable extensions.
# (Set globally, since these properties are not inherited by targets.)
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 23)
endif()
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)

# Lists of compile options.
set(StormRuler_CXX_COMPILE_OPTIONS "")
set(StormRuler_CXX_DEBUG_COMPILE_OPTIONS "")
set(StormRuler_CXX_RELEASE_COMPILE_OPTIONS "")
set(StormRuler_CXX_COVERAGE_COMPILE_OPTIONS "")

# Lists of link options.
set(StormRuler_CXX_LINK_OPTIONS "")
set(StormRuler_CXX_DEBUG_LINK_OPTIONS "")
set(StormRuler_CXX_RELEASE_LINK_OPTIONS "")
set(StormRuler_CXX_COVERAGE_LINK_OPTIONS "")

# Set the general options.
if(Storm_CXX_COMPILER_GNU_LIKE)
  # Extended concepts diagnostics.
  list(APPEND StormRuler_CXX_COMPILE_OPTIONS -fconcepts-diagnostics-depth=3)
elseif(Storm_CXX_COMPILER_CLANG_LIKE)
endif()

# Set the warnings and diagnostics options.
if(Storm_CXX_COMPILER_GNU_LIKE OR Storm_CXX_COMPILER_CLANG_LIKE)
  # Enable maaany warnings, treated as errors.
  list(APPEND StormRuler_CXX_COMPILE_OPTIONS 
       -Werror -Wall -Wextra -Wpedantic -Wno-unused-function -Wno-psabi)
  # TODO: to be removed!
  list(APPEND StormRuler_CXX_COMPILE_OPTIONS -Wno-unused-parameter)
endif()

# Set the optimization options.
if(Storm_CXX_COMPILER_GNU_LIKE)
  # Debug: debug information, no optimizations.
  list(APPEND StormRuler_CXX_DEBUG_COMPILE_OPTIONS -g -O0)
  # Release: debug information, full optimizations.
  list(APPEND StormRuler_CXX_RELEASE_COMPILE_OPTIONS -g -Ofast -march=native)
  # Coverage: coverage, debug information, no optimizations, no inlining.
  list(APPEND StormRuler_CXX_COVERAGE_COMPILE_OPTIONS 
       --coverage -g -O0 -fno-inline -fno-default-inline
       -fno-inline-small-functions -fno-elide-constructors)
  list(APPEND StormRuler_CXX_COVERAGE_LINK_OPTIONS --coverage) 
elseif(Storm_CXX_COMPILER_CLANG_LIKE)
endif()

# Apply the options (for all build types).
target_compile_options(
    StormRuler_BASE INTERFACE 
    ${StormRuler_CXX_COMPILE_OPTIONS})
target_link_options(
    StormRuler_BASE INTERFACE 
    ${StormRuler_CXX_LINK_OPTIONS})

# Apply the options (for Debug build type).
target_compile_options(
    StormRuler_BASE INTERFACE
    $<$<CONFIG:Debug>:${StormRuler_CXX_DEBUG_COMPILE_OPTIONS}>)
target_link_options(
    StormRuler_BASE INTERFACE
    $<$<CONFIG:Debug>:${StormRuler_CXX_DEBUG_LINK_OPTIONS}>)

# Apply the options (for Release build type).
target_compile_options(
    StormRuler_BASE INTERFACE
    $<$<CONFIG:Release>:${StormRuler_CXX_RELEASE_COMPILE_OPTIONS}>)
target_link_options(
    StormRuler_BASE INTERFACE
    $<$<CONFIG:Release>:${StormRuler_CXX_RELEASE_LINK_OPTIONS}>)

# Apply the options (for Coverage build type).
target_compile_options(
    StormRuler_BASE INTERFACE
    $<$<CONFIG:Coverage>:${StormRuler_CXX_COVERAGE_COMPILE_OPTIONS}>)
target_link_options(
    StormRuler_BASE INTERFACE
    $<$<CONFIG:Coverage>:${StormRuler_CXX_COVERAGE_LINK_OPTIONS}>)

# ------------------------------------------------------------------------------

message(STATUS "Looking for libraries..")

# Add the CMake modules.
list(APPEND CMAKE_MODULE_PATH "./cmake")

# Find Python interpreter.
find_package(Python3 REQUIRED COMPONENTS Interpreter)

# ------------------------------------------------------------------------------

# Find {fmt}.
find_package(fmt CONFIG REQUIRED)
target_link_libraries(StormRuler_BASE INTERFACE fmt::fmt-header-only)

# Find spdlog.
find_package(spdlog CONFIG REQUIRED)
target_link_libraries(StormRuler_BASE INTERFACE spdlog::spdlog_header_only)

# ------------------------------------------------------------------------------

# Find TBB.
if(StormRuler_NO_TBB)
  message(STATUS "| TBB was disabled by user.")
else()
  find_package(TBB QUIET)
  if(TBB_FOUND)
    set(STORM_TBB_ENABLED 1)
    message(STATUS "| TBB enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE TBB::tbb)
  endif()
endif()

# Find OpenMP.
if(StormRuler_NO_OpenMP)
  message(STATUS "| OpenMP was disabled by user.")
else()
  find_package(OpenMP QUIET)
  if(OpenMP_FOUND)
    set(STORM_OpenMP_ENABLED 1)
    message(STATUS "| OpenMP enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE OpenMP::OpenMP_CXX)
  endif()
endif()

# ------------------------------------------------------------------------------

# Find MKL.
if(StormRuler_NO_MKL)
  message(STATUS "| MKL was disabled by user.")
else()
  # Configure OpenMP library to the Fortran compiler native one.
  #if(STORM_OpenMP_ENABLED)
  #  if(CMAKE_CXX_COMPILER_ID STREQUAL "INTEL")
  #    set(MKL_THREADING "intel_thread")
  #  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  #    set(MKL_THREADING "gnu_thread")
  #  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
  #    set(MKL_THREADING "pgi_thread")
  #  endif()
  #endif()
  # Configure MKL to use 32-bit integers.
  set(MKL_INTERFACE lp64)
  find_package(MKL CONFIG QUIET)
  if(MKL_FOUND)
    set(STORM_MKL_ENABLED 1)
    message(STATUS "| MKL enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE MKL::MKL)
  endif()
endif() 

# Find BLAS.
if(StormRuler_NO_BLAS)
  message(STATUS "| BLAS was disabled by user.")
elseif(STORM_MKL_ENABLED)
  # Separate BLAS is not needed in case MKL is used.
  set(STORM_BLAS_ENABLED 1)
  message(STATUS "| BLAS enabled (MKL).")
else()
  find_package(BLAS QUIET)
  if(BLAS_FOUND)
    set(STORM_BLAS_ENABLED 1)
    message(STATUS "| BLAS enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE BLAS::BLAS)
  endif()
endif()

# Find LAPACK.
if(StormRuler_NO_LAPACK)
  message(STATUS "| LAPACK was disabled by user.")
elseif(STORM_MKL_ENABLED)
  # Separate LAPACK is not needed in case MKL is used.
  set(STORM_LAPACK_ENABLED 1)
  message(STATUS "| LAPACK enabled (MKL).")
else()
  find_package(LAPACK QUIET)
  if(LAPACK_FOUND)
    set(STORM_LAPACK_ENABLED 1)
    message(STATUS "| LAPACK enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE LAPACK::LAPACK)
  endif()
endif()

# Find METIS.
if(StormRuler_NO_METIS)
  message(STATUS "| METIS was disabled by user.")
else()
  find_package(METIS QUIET)
  if(METIS_FOUND)
    message(STATUS "| METIS enabled.")
    set(STORM_METIS_ENABLED 1)
    target_link_libraries(StormRuler_BASE INTERFACE METIS)
  endif()
endif()

# ------------------------------------------------------------------------------

# Find ZLIB.
if(StormRuler_NO_ZLIB)
  message(STATUS "| ZLIB was disabled by user.")
else()
  find_package(ZLIB QUIET)
  if(ZLIB_FOUND)
    message(STATUS "| ZLIB enabled.")
    set(STORM_ZLIB_ENABLED 1)
    target_link_libraries(StormRuler_BASE INTERFACE ZLIB::ZLIB)
  endif()
endif()

# Find OpenGL, GLEW, GLFW and GLM.
if(StormRuler_NO_OpenGL)
  message(STATUS "| OpenGL was disabled by user.")
else()
  find_package(OpenGL QUIET)
  if(OpenGL_FOUND)
    message(STATUS "| OpenGL enabled.")
    set(STORM_OpenGL_ENABLED 1)
    target_link_libraries(StormRuler_BASE INTERFACE OpenGL::GL)
  endif()
  find_package(GLEW QUIET)
  if(GLEW_FOUND)
    set(STORM_GLEW_ENABLED 1)
    message(STATUS "| GLEW enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE GLEW::GLEW)
  endif()
  find_package(glfw3 3.3 QUIET)
  if(glfw3_FOUND)
    set(STORM_GLFW_ENABLED 1)
    message(STATUS "| GLFW enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE glfw)
  endif()
  find_package(glm QUIET)
  if(glm_FOUND)
    set(STORM_GLM_ENABLED 1)
    message(STATUS "| GLM enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE glm::glm)
  endif()
endif()

# ------------------------------------------------------------------------------

# Add the core library.
message(STATUS "Adding StormRuler core library..")
add_subdirectory("./source")

# ------------------------------------------------------------------------------

# Add the applications.
message(STATUS "Adding StormRuler applications..")
add_subdirectory("./source_apps")

# ------------------------------------------------------------------------------

# Add the documentation.
message(STATUS "Adding StormRuler documentation..")
add_subdirectory("./source_docs")

# ------------------------------------------------------------------------------

# Add the tests.
message(STATUS "Adding StormRuler tests..")
add_subdirectory("./tests")

# ------------------------------------------------------------------------------
