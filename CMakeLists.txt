# Copyright (C) 2022 Oleg Butakov
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
# • User friendly compilation script:
#   It is inconvenient to write `CXX=.. cmake . -DCMAKE_BUILD_TYPE=..`
#   every time.
# • Problems with GTest linking on ARM-based Mac machines.

cmake_minimum_required(VERSION 3.20..3.26)

# -----------------------------------------------------------------------------

# Setup the project.
project(StormRuler
        VERSION 0.1
        DESCRIPTION "A very high order CFD solver"
        LANGUAGES CXX)

# A container library that stores all the properties.
add_library(StormRuler_BASE INTERFACE)

# Print the banner.
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

# -----------------------------------------------------------------------------

# Set output directories.
# (Set globally, since these properties are not inherited by targets.)
set(CMAKE_PDB_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")

# Require C++23 and disable extensions.
# (Set globally, since these properties are not inherited by targets.)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)

# Check build type is valid.
if(NOT (CMAKE_BUILD_TYPE STREQUAL "Debug" OR
        CMAKE_BUILD_TYPE STREQUAL "Release"))
  message(FATAL_ERROR "Invalid or unspecified build type!")
endif()

# Check C/C++ compiler.
if(NOT (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
        CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
  message(WARNING "Unknown C++ compiler, correct options may not be set!")
endif()

# Compiler warning and diagnostics.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  list(APPEND StormRuler_CXX_OPTIONS -Werror -Wall -Wextra -Wpedantic 
                                     -Wno-unused-function 
                                     -Wno-unused-parameter 
                                     -Wno-comment 
                                     -Wno-unknown-pragmas 
                                     -fconcepts-diagnostics-depth=3)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  list(APPEND StormRuler_CXX_OPTIONS -Werror -Wall -Wextra -Wpedantic 
                                     -Wno-unused-function 
                                     -Wno-unused-parameter 
                                     -Wno-mismatched-tags
                                     -Wno-gnu-zero-variadic-macro-arguments)
endif()

# Compiler optimization options.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  list(APPEND StormRuler_CXX_DEBUG_OPTIONS -g -O0)
  list(APPEND StormRuler_CXX_RELEASE_OPTIONS -g -Ofast -march=native) 
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  list(APPEND StormRuler_CXX_DEBUG_OPTIONS -g -O0)
  list(APPEND StormRuler_CXX_RELEASE_OPTIONS -g -Ofast -march=native)
endif()

# Set the C/C++ options.
target_compile_options(StormRuler_BASE INTERFACE 
                       ${StormRuler_CXX_OPTIONS})
target_compile_options(StormRuler_BASE INTERFACE
                       $<$<CONFIG:Debug>:${StormRuler_CXX_DEBUG_OPTIONS}>)
target_compile_options(StormRuler_BASE INTERFACE
                       $<$<CONFIG:Release>:${StormRuler_CXX_RELEASE_OPTIONS}>)

# -----------------------------------------------------------------------------

# Add the CMake modules.
list(APPEND CMAKE_MODULE_PATH "./cmake")

# Find Git and update submodules.
find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    RESULT_VARIABLE GIT_SUBMODULE_RESULT)
  if(NOT GIT_SUBMODULE_RESULT EQUAL "0")
    message(FATAL_ERROR "`git submodule update --init --recursive`"
                        " failed with ${GIT_SUBMODULE_RESULT}.")
  endif()
endif()

# Find Python interpreter.
find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)

# -----------------------------------------------------------------------------

# Find {fmt}.
find_package(fmt CONFIG REQUIRED)
target_link_libraries(StormRuler_BASE INTERFACE fmt::fmt-header-only)

# Find spdlog.
find_package(spdlog CONFIG REQUIRED)
target_link_libraries(StormRuler_BASE INTERFACE spdlog::spdlog_header_only)

# -----------------------------------------------------------------------------

# Find TBB.
if(StormRuler_NO_TBB)
  message(STATUS "StormRuler TBB was disabled by user.")
else()
  find_package(TBB)
  if(TBB_FOUND)
    set(STORM_TBB_ENABLED 1)
    message(STATUS "TBB enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE TBB::tbb)
  endif()
endif()

# Find OpenMP.
if(StormRuler_NO_OpenMP)
  message(STATUS "StormRuler OpenMP was disabled by user.")
else()
  find_package(OpenMP)
  if(OpenMP_FOUND)
    set(STORM_OpenMP_ENABLED 1)
    message(STATUS "OpenMP enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE OpenMP::OpenMP_CXX)
  endif()
endif()

# -----------------------------------------------------------------------------

# Find MKL.
if(StormRuler_NO_MKL)
  message(STATUS "StormRuler MKL was disabled by user.")
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
  find_package(MKL CONFIG)
  if(MKL_FOUND)
    set(STORM_MKL_ENABLED 1)
    message(STATUS "MKL enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE MKL::MKL)
  endif()
endif() 

# Find BLAS.
if(StormRuler_NO_BLAS)
  message(STATUS "StormRuler BLAS was disabled by user.")
elseif(STORM_MKL_ENABLED)
  # Separate BLAS is not needed in case MKL is used.
  set(STORM_BLAS_ENABLED 1)
  message(STATUS "BLAS enabled (MKL).")
else()
  find_package(BLAS)
  if(BLAS_FOUND)
    set(STORM_BLAS_ENABLED 1)
    message(STATUS "BLAS enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE BLAS::BLAS)
  endif()
endif()

# Find LAPACK.
if(StormRuler_NO_LAPACK)
  message(STATUS "StormRuler LAPACK was disabled by user.")
elseif(STORM_MKL_ENABLED)
  # Separate LAPACK is not needed in case MKL is used.
  set(STORM_LAPACK_ENABLED 1)
  message(STATUS "LAPACK enabled (MKL).")
else()
  find_package(LAPACK)
  if(LAPACK_FOUND)
    set(STORM_LAPACK_ENABLED 1)
    message(STATUS "LAPACK enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE LAPACK::LAPACK)
  endif()
endif()

# Find METIS.
if(StormRuler_NO_METIS)
  message(STATUS "StormRuler METIS was disabled by user.")
else()
  find_package(METIS)
  if(METIS_FOUND)
    message(STATUS "METIS enabled.")
    set(STORM_METIS_ENABLED 1)
    target_link_libraries(StormRuler_BASE INTERFACE METIS)
  endif()
endif()

# -----------------------------------------------------------------------------

# Find ZLIB.
if(StormRuler_NO_ZLIB)
  message(STATUS "StormRuler ZLIB was disabled by user.")
else()
  find_package(ZLIB)
  if(ZLIB_FOUND)
    message(STATUS "ZLIB enabled.")
    set(STORM_ZLIB_ENABLED 1)
    target_link_libraries(StormRuler_BASE INTERFACE ZLIB::ZLIB)
  endif()
endif()

# Find OpenGL, GLEW, GLFW and GLM.
if(StormRuler_NO_OpenGL)
  message(STATUS "StormRuler OpenGL was disabled by user.")
else()
  find_package(OpenGL)
  if(OpenGL_FOUND)
    message(STATUS "OpenGL enabled.")
    set(STORM_OpenGL_ENABLED 1)
    target_link_libraries(StormRuler_BASE INTERFACE OpenGL::GL)
  endif()
  find_package(GLEW)
  if(GLEW_FOUND)
    set(STORM_GLEW_ENABLED 1)
    message(STATUS "GLEW enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE GLEW::GLEW)
  endif()
  find_package(glfw3 3.3)
  if(glfw3_FOUND)
    set(STORM_GLFW_ENABLED 1)
    message(STATUS "GLFW enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE glfw)
  endif()
  find_package(glm)
  if(glm_FOUND)
    set(STORM_GLM_ENABLED 1)
    message(STATUS "GLM enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE glm::glm)
  endif()
endif()

# -----------------------------------------------------------------------------

# Add the core library.
message(STATUS "Adding StormRuler core library..")
add_subdirectory("./source")

# Add the applications.
message(STATUS "Adding StormRuler applications..")
add_subdirectory("./source_apps")

# Add the documentation.
message(STATUS "Adding StormRuler documentation..")
add_subdirectory("./source_docs")

# Add the tests.
message(STATUS "Adding StormRuler tests..")
add_subdirectory("./tests")
