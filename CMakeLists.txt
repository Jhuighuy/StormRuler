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

# TODO:
# • more native mold integration.
# • clang-format integration:
#   Format the code on build.
# • clang-tidy integration:
#   Impossible now, since clang-tidy doesn't support C++23.
#   (But now we support C++20.. again)
# • Add sanitizers for the Coverage build.
# • Share project description with `source_docs`.

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

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# ------------------------------------------------------------------------------

# Setup the project.
project(StormRuler
        VERSION 0.1
        DESCRIPTION "A very high order FVM framework"
        LANGUAGES CXX)

# A container library that stores all the properties.
add_library(StormRuler_Base INTERFACE)
add_library(StormRuler::Base ALIAS StormRuler_Base)

# ------------------------------------------------------------------------------

# Include directories.
target_include_directories(StormRuler_Base INTERFACE "./source")

# ------------------------------------------------------------------------------

# Set output directories.
# (Set globally, since these properties are not inherited by targets.)
set(CMAKE_PDB_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# ------------------------------------------------------------------------------

# Check C/C++ compiler.
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(CXX_COMPILER_GNU_LIKE YES)
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "12.1.0")
    message(FATAL_ERROR "GCC version 12.1 or more recent is required!")
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  set(CXX_COMPILER_CLANG_LIKE YES)
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "16.0.0")
    message(FATAL_ERROR "Clang version 16.0 or more recent is required!")
  endif()
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  set(CXX_COMPILER_MSVC_LIKE YES)
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "19.34.0")
    message(FATAL_ERROR "MSVC version 19.34.0 or more recent is required!")
  endif()
else()
  message(FATAL_ERROR "Unsupported C++ compiler!")
endif()

# Require C++20 and disable extensions.
# (Set globally, since these properties are not inherited by targets.)
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
elseif(CMAKE_CXX_STANDARD LESS 20)
  message(FATAL_ERROR "C++20 or more recent is required!")
endif()
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_EXTENSIONS NO)

# Lists of defines.
set(CXX_DEFINITIONS "")
set(CXX_DEBUG_DEFINITIONS "")
set(CXX_RELEASE_DEFINITIONS "")
set(CXX_COVERAGE_DEFINITIONS "")

# Lists of compile options.
set(CXX_COMPILE_OPTIONS "")
set(CXX_DEBUG_COMPILE_OPTIONS "")
set(CXX_RELEASE_COMPILE_OPTIONS "")
set(CXX_COVERAGE_COMPILE_OPTIONS "")

# Lists of link options.
set(CXX_LINK_OPTIONS "")
set(CXX_DEBUG_LINK_OPTIONS "")
set(CXX_RELEASE_LINK_OPTIONS "")
set(CXX_COVERAGE_LINK_OPTIONS "")

# Set the general options.
if(CXX_COMPILER_GNU_LIKE)
  # Enable checks in STL.
  list(APPEND CXX_DEBUG_DEFINITIONS _GLIBCXX_DEBUG)
  list(APPEND CXX_COVERAGE_DEFINITIONS _GLIBCXX_DEBUG)
  # Extended concepts diagnostics.
  list(APPEND CXX_COMPILE_OPTIONS -fconcepts-diagnostics-depth=3)
elseif(CXX_COMPILER_CLANG_LIKE)
elseif(CXX_COMPILER_MSVC_LIKE)
  # Enable new preprocessor.
  list(APPEND CXX_COMPILE_OPTIONS /Zc:preprocessor)
endif()
list(APPEND CXX_COVERAGE_DEFINITIONS STORM_COVERAGE=1)

# Set the warnings and diagnostics options.
if(CXX_COMPILER_GNU_LIKE OR CXX_COMPILER_CLANG_LIKE)
  # Enable maaany warnings, treated as errors.
  list(APPEND CXX_COMPILE_OPTIONS
       -Werror -Wall -Wextra -Wpedantic -Wno-unused-function -Wno-psabi)
  # TODO: to be removed!
  list(APPEND CXX_COMPILE_OPTIONS -Wno-unused-parameter)
elseif(CXX_COMPILER_MSVC_LIKE)
  # Enable maaany warnings, treated as errors.
  list(APPEND CXX_COMPILE_OPTIONS
       /WX /Wall
       /wd4061  # Enumerator '...' in switch of enum '...' is not explicitly
                # handled by a case label.
       /wd4365  # Signed/unsigned mismatch.
       /wd4458  # Declaration of '...' hides class member.
       /wd4459  # Declaration of '...' hides global declaration.
       /wd4464  # Relative include path contains '...'.
       /wd4514  # Unreferenced inline function has been removed.
       /wd4557  # '__assume' contains effect '...'.
       /wd4582  # Constructor is not implicitly called.
       /wd4623  # Default constructor was implicitly defined as deleted.
       /wd4625  # Copy constructor was implicitly defined as deleted.
       /wd4626  # Assignment operator was implicitly defined as deleted.
       /wd4710  # '...': function not inlined.
       /wd4711  # Function '...' selected for automatic inline expansion.
       /wd4820  # '...' bytes padding added after data member '...'.
       /wd4866  # Compiler may not enforce left-to-right evaluation order
                # for call to '...'.
       /wd4868  # Compiler may not enforce left-to-right evaluation order in
                # braced initializer list.
       /wd5026  # Move constructor was implicitly defined as deleted.
       /wd5027  # Move assignment operator was implicitly defined as deleted.
       /wd5039  # Pointer or reference to potentially throwing function passed
                # to extern C function under -EHc. Undefined behavior may occur
                # if this function throws an exception.
       /wd5045  # Compiler will insert Spectre mitigation for memory load if
                # /Qspectre switch specified.
       /wd5105  # Macro expansion producing 'defined' has undefined behavior
       /wd5214  # Applying '...' to an operand with a volatile qualified type
                # is deprecated in C++20.
       /wd5246  # The initialization of a subobject should be wrapped in braces
       /wd5262  # Implicit fall-through occurs here.
       /wd5264) # 'const' variable is not used.
  # TODO: to be removed!
  list(APPEND CXX_COMPILE_OPTIONS
       /wd4100  # Unused parameter.
       /wd5219) # Implicit conversion from '..' to '..', possible loss of data.
endif()

# Set the optimization options.
if(CXX_COMPILER_GNU_LIKE)
  # Debug: debug information, no optimizations.
  list(APPEND CXX_DEBUG_COMPILE_OPTIONS -g -O0)
  # Release: debug information, full optimizations.
  list(APPEND CXX_RELEASE_COMPILE_OPTIONS -g -Ofast
       -march=native -finline-limit=10000000 -ffp-contract=fast)
  # Coverage: coverage, debug information, no optimizations, no inlining.
  list(APPEND CXX_COVERAGE_COMPILE_OPTIONS
       --coverage -g -O0 -fno-inline -fno-default-inline
       -fno-inline-small-functions -fno-elide-constructors)
  list(APPEND CXX_COVERAGE_LINK_OPTIONS --coverage)
elseif(CXX_COMPILER_CLANG_LIKE)
elseif(CXX_COMPILER_MSVC_LIKE)
  # TODO: https://gitlab.kitware.com/cmake/cmake/-/issues/20812
  # Debug: no optimizations.
  # TODO: debug information?
  list(APPEND CXX_DEBUG_COMPILE_OPTIONS /Od)
  # Release: full optimizations.
  # TODO: AVX512?
  list(APPEND CXX_RELEASE_COMPILE_OPTIONS /O2 /arch:AVX2)
  # TODO: Coverage options.
endif()

# Apply the options.
target_compile_definitions(
    StormRuler_Base INTERFACE
    ${CXX_DEFINITIONS}
    $<$<CONFIG:Debug>:${CXX_DEBUG_DEFINITIONS}>
    $<$<CONFIG:Release>:${CXX_RELEASE_DEFINITIONS}>
    $<$<CONFIG:Coverage>:${CXX_COVERAGE_DEFINITIONS}>)
target_compile_options(
    StormRuler_Base INTERFACE
    ${CXX_COMPILE_OPTIONS}
    $<$<CONFIG:Debug>:${CXX_DEBUG_COMPILE_OPTIONS}>
    $<$<CONFIG:Release>:${CXX_RELEASE_COMPILE_OPTIONS}>
    $<$<CONFIG:Coverage>:${CXX_COVERAGE_COMPILE_OPTIONS}>)
target_link_options(
    StormRuler_Base INTERFACE
    ${CXX_LINK_OPTIONS}
    $<$<CONFIG:Debug>:${CXX_DEBUG_LINK_OPTIONS}>
    $<$<CONFIG:Release>:${CXX_RELEASE_LINK_OPTIONS}>
    $<$<CONFIG:Coverage>:${CXX_COVERAGE_LINK_OPTIONS}>)

# ------------------------------------------------------------------------------

# Find Python interpreter.
# (Here we wrap `find_package` into a separate function
#  to avoid any problems when pybind11 would look for Python itself later.)
function(_find_python_interpreter)
  find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter QUIET)
  set(Python3_EXECUTABLE ${Python3_EXECUTABLE} PARENT_SCOPE)
endfunction()
_find_python_interpreter()
unset(_find_python_interpreter)

# Custom configuration function.
function(storm_configure_file INPUT_FILE OUTPUT_FILE DEFAULTS_OUTPUT_FILE)
  configure_file(${INPUT_FILE} ${OUTPUT_FILE})
  add_custom_command(
      OUTPUT ${DEFAULTS_OUTPUT_FILE}
      COMMAND ${Python3_EXECUTABLE}
              ${CMAKE_SOURCE_DIR}/scripts/make_config_defaults.py
              ${INPUT_FILE} ${DEFAULTS_OUTPUT_FILE}
      MAIN_DEPENDENCY ${INPUT_FILE}
      VERBATIM)
endfunction()

# ------------------------------------------------------------------------------

# Find the mold linker.
if(STORM_NO_MOLD)
  message(STATUS "| mold was disabled by user.")
else()
  # TODO: once `mold` has a CMake integration.
  # find_package(mold CONFIG QUIET)
  # if(mold_FOUND)
  #   message(STATUS "| mold enabled.")
  #   target_use_mold(StormRuler_Base)
  # endif()
  # For now we have to use this:
  if(NOT APPLE AND (CXX_COMPILER_GNU_LIKE OR CXX_COMPILER_CLANG_LIKE))
    find_program(MOLD_EXECUTABLE NAMES mold)
    if(MOLD_EXECUTABLE)
      message(STATUS "| mold enabled.")
      target_link_options(StormRuler_Base INTERFACE -fuse-ld=mold)
    endif()
  endif()
endif()

# ------------------------------------------------------------------------------

# Find {fmt}.
find_package(fmt CONFIG REQUIRED QUIET)
target_link_libraries(StormRuler_Base INTERFACE fmt::fmt-header-only)

# Find spdlog.
find_package(spdlog CONFIG REQUIRED QUIET)
target_link_libraries(StormRuler_Base INTERFACE spdlog::spdlog_header_only)

# ------------------------------------------------------------------------------

# Find TBB.
if(STORM_NO_TBB)
  message(STATUS "| TBB was disabled by user.")
else()
  find_package(TBB QUIET)
  if(TBB_FOUND)
    set(STORM_TBB_ENABLED 1)
    message(STATUS "| TBB enabled.")
    target_link_libraries(StormRuler_Base INTERFACE TBB::tbb)
  endif()
endif()

# Find OpenMP.
if(STORM_NO_OpenMP)
  message(STATUS "| OpenMP was disabled by user.")
else()
  find_package(OpenMP QUIET)
  if(OpenMP_FOUND)
    set(STORM_OpenMP_ENABLED 1)
    message(STATUS "| OpenMP enabled.")
    target_link_libraries(StormRuler_Base INTERFACE OpenMP::OpenMP_CXX)
  endif()
endif()

# ------------------------------------------------------------------------------

# Find MKL.
if(STORM_NO_MKL)
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
    target_link_libraries(StormRuler_Base INTERFACE MKL::MKL)
  endif()
endif()

# Find BLAS.
if(STORM_NO_BLAS)
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
    target_link_libraries(StormRuler_Base INTERFACE BLAS::BLAS)
  endif()
endif()

# Find LAPACK.
if(STORM_NO_LAPACK)
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
    target_link_libraries(StormRuler_Base INTERFACE LAPACK::LAPACK)
  endif()
endif()

# Find METIS.
if(STORM_NO_METIS)
  message(STATUS "| METIS was disabled by user.")
else()
  # `FindMETIS.cmake` looks broken.
  #find_package(METIS QUIET)
  #if(METIS_FOUND)
  #  message(STATUS "| METIS enabled.")
  #  set(STORM_METIS_ENABLED 1)
  #  target_link_libraries(StormRuler_Base INTERFACE METIS)
  #endif()
endif()

# ------------------------------------------------------------------------------

# Find ZLIB.
if(STORM_NO_ZLIB)
  message(STATUS "| ZLIB was disabled by user.")
else()
  find_package(ZLIB QUIET)
  if(ZLIB_FOUND)
    message(STATUS "| ZLIB enabled.")
    set(STORM_ZLIB_ENABLED 1)
    target_link_libraries(StormRuler_Base INTERFACE ZLIB::ZLIB)
  endif()
endif()

# ------------------------------------------------------------------------------

# Find pybind11.
if(STORM_NO_PYTHON)
  message(STATUS "| Python was disabled by user.")
else()
  find_package(pybind11 QUIET)
  if(pybind11_FOUND)
    set(STORM_PYTHON_ENABLED 1)
    message(STATUS "| Python enabled.")
    target_link_libraries(StormRuler_Base INTERFACE pybind11::pybind11)
  endif()
endif()

# ------------------------------------------------------------------------------

# Find OpenGL, GLEW, GLFW and GLM.
if(STORM_NO_OpenGL)
  message(STATUS "| OpenGL was disabled by user.")
else()
  find_package(OpenGL QUIET)
  if(OpenGL_FOUND)
    message(STATUS "| OpenGL enabled.")
    set(STORM_OpenGL_ENABLED 1)
    target_link_libraries(StormRuler_Base INTERFACE OpenGL::GL)
  endif()
  find_package(GLEW QUIET)
  if(GLEW_FOUND)
    set(STORM_GLEW_ENABLED 1)
    message(STATUS "| GLEW enabled.")
    target_link_libraries(StormRuler_Base INTERFACE GLEW::GLEW)
  endif()
  find_package(glfw3 3.3 QUIET)
  if(glfw3_FOUND)
    set(STORM_GLFW_ENABLED 1)
    message(STATUS "| GLFW enabled.")
    target_link_libraries(StormRuler_Base INTERFACE glfw)
  endif()
  find_package(glm QUIET)
  if(glm_FOUND)
    set(STORM_GLM_ENABLED 1)
    message(STATUS "| GLM enabled.")
    target_link_libraries(StormRuler_Base INTERFACE glm::glm)
  endif()
endif()

# ------------------------------------------------------------------------------

# Add the core library.
message(STATUS "Adding StormRuler libraries..")
add_subdirectory(./source)

# ------------------------------------------------------------------------------

# Add the applications.
message(STATUS "Adding StormRuler applications..")
add_subdirectory(./source_apps)

# ------------------------------------------------------------------------------

# Add the tests.
message(STATUS "Adding StormRuler tests..")
add_subdirectory(./tests)

# ------------------------------------------------------------------------------

# Add the documentation.
message(STATUS "Adding StormRuler documentation..")
add_subdirectory(./source_docs)

# ------------------------------------------------------------------------------
