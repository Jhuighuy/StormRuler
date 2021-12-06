## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ##
## Copyright (C) 2021 Oleg Butakov
## 
## Permission is hereby granted, free of charge, to any person 
## obtaining a copy of this software and associated documentation 
## files (the "Software"), to deal in the Software without 
## restriction, including without limitation the rights  to use, 
## copy, modify, merge, publish, distribute, sublicense, and/or
## sell copies of the Software, and to permit persons to whom the  
## Software is furnished to do so, subject to the following 
## conditions:
## 
## The above copyright notice and this permission notice shall be 
## included in all copies or substantial portions of the Software.
## 
## THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
## EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
## OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
## HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
## WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
## FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
## OTHER DEALINGS IN THE SOFTWARE.
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ##

cmake_minimum_required(VERSION 3.12)

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Setup StormRuler project, using Fortran and C++.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

project(StormRuler
  VERSION 1.0
  DESCRIPTION "A very high order CFD solver and beyond"
  LANGUAGES Fortran CXX)

# A container library that stores all the properties.
add_library(StormRuler_BASE INTERFACE)

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Compilation options.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

# Detect Apple Silicon CPU.
# Fortran compilers on Apple Silicon are not production-ready yet,
# so we need to disable some features (for now).
if(APPLE)
  include(CheckCxxCompilerFlag)
  check_cxx_compiler_flag("-arch arm64" APPLE_SILICON)
  message(STATUS "Apple Silicon detected, this is bad..")
endif()

# Warning and diagnostics options (Fortran).
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(StormRuler_BASE INTERFACE
    $<$<COMPILE_LANGUAGE:Fortran>:
      -Wall -Wextra -Wpedantic -Wno-unused-dummy-argument>)
  target_compile_options(StormRuler_BASE INTERFACE
    $<$<COMPILE_LANGUAGE:Fortran>:$<$<CONFIG:Debug>:-fcheck=all>>)
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  target_compile_options(StormRuler_BASE INTERFACE
    $<$<COMPILE_LANGUAGE:Fortran>: -warn all,nounused>)
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI")
    #target_compile_options(StormRuler_BASE INTERFACE
    #  $<$<COMPILE_LANGUAGE:Fortran>: -warn all,nounused>)
else()
  message(WARNING 
    "Unknown Fortran compiler ${CMAKE_Fortran_COMPILER_ID}, warning options not set.")
endif()
# Warning and diagnostics options (C/C++).
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(StormRuler_BASE INTERFACE
    $<$<COMPILE_LANGUAGE:C,CXX>: -Werror -Wall -Wextra 
       -Wpedantic -Wno-unused-function -Wno-unused-parameter>)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR
       "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  target_compile_options(StormRuler_BASE INTERFACE
    $<$<COMPILE_LANGUAGE:C,CXX>:-Werror -Wall -Wextra -Wpedantic 
      -Wno-unused-function -Wno-unused-parameter 
      -Wno-gnu-zero-variadic-macro-arguments>)
else()
  message(WARNING 
    "Unknown C++ compiler ${CMAKE_CXX_COMPILER_ID}, warning options not set.")
endif()

# Optimization options (Fortran).
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(StormRuler_BASE INTERFACE
    $<$<COMPILE_LANGUAGE:Fortran>:$<$<CONFIG:Debug>:-g -O0>>)
  if (NOT APPLE_SILICON)
    target_compile_options(StormRuler_BASE INTERFACE
      $<$<COMPILE_LANGUAGE:Fortran>:$<$<CONFIG:Release>:-g -Ofast -march=native>>)
  endif()
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
  target_compile_options(StormRuler_BASE INTERFACE
    $<$<COMPILE_LANGUAGE:Fortran>:$<$<CONFIG:Debug>:-g -O0>>)
  target_compile_options(StormRuler_BASE INTERFACE
    $<$<COMPILE_LANGUAGE:Fortran>:$<$<CONFIG:Release>:-g -Ofast -march=native>>)
elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI")
    target_compile_options(StormRuler_BASE INTERFACE
      $<$<COMPILE_LANGUAGE:Fortran>:$<$<CONFIG:Debug>:-g -C -O0>>)
    target_compile_options(StormRuler_BASE INTERFACE
      $<$<COMPILE_LANGUAGE:Fortran>:$<$<CONFIG:Release>:-g -O4 -march=native>>)
else()
  message(WARNING 
    "Unknown Fortran compiler ${CMAKE_Fortran_COMPILER_ID}, optimization options not set.")
endif()
# Optimization options (C/C++).
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(StormRuler_BASE INTERFACE
    $<$<COMPILE_LANGUAGE:C,CXX>:$<$<CONFIG:Debug>:-g -O0>>)
  target_compile_options(StormRuler_BASE INTERFACE
    $<$<COMPILE_LANGUAGE:C,CXX>:$<$<CONFIG:Release>:-g -Ofast -march=native>>)
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR
       "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  target_compile_options(StormRuler_BASE INTERFACE
    $<$<COMPILE_LANGUAGE:C,CXX>:$<$<CONFIG:Debug>:-g -O0>>)
  target_compile_options(StormRuler_BASE INTERFACE
    $<$<COMPILE_LANGUAGE:C,CXX>:$<$<CONFIG:Release>:-g -Ofast -march=native>>)
else()
  message(WARNING 
    "Unknown C++ compiler ${CMAKE_CXX_COMPILER_ID}, warning options not set.")
endif()

# Other options.
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
  # NAG does not use byte-based kind system by default.
  target_compile_options(StormRuler_BASE INTERFACE 
    $<$<COMPILE_LANGUAGE:Fortran>:-kind=byte>)
endif()

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Libraries.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

# CMake modules.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMakeModules")

# Library include directories.
target_include_directories(
  StormRuler_BASE INTERFACE "${CMAKE_SOURCE_DIR}/src")
target_include_directories(
  StormRuler_BASE INTERFACE "${CMAKE_SOURCE_DIR}/ext")
  
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

# Find OpenMP.
if(StormRuler_NO_OpenMP)
  message(STATUS "StormRuler OpenMP was disabled by user.")
else()
  find_package(OpenMP)
  if(OpenMP_Fortran_FOUND)
    # TODO: NAG has linking problems for OpenMP builds.
    if(NOT "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
      set(OpenMP_Fortran_ENABLED 1)
    endif()
    if(OpenMP_Fortran_ENABLED)
      message(STATUS "Fortran OpenMP enabled.")
      target_link_libraries(StormRuler_BASE INTERFACE OpenMP::OpenMP_Fortran)
    endif()
  endif()
  if(OpenMP_CXX_FOUND)
    if(OpenMP_Fortran_ENABLED)
      # Enable C++ OpenMP only in case Fortran and C++ compiler vendors match.
      # Mixing different OpenMP runtimes is generally not a good idea.
      if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "${CMAKE_CXX_COMPILER_ID}")
        set(OpenMP_CXX_ENABLED 1)
      endif()
      # NAG uses system C compiler as a backend,
      # which can should theoretically match used C compiler.
      if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
        set(OpenMP_CXX_ENABLED 1)
      endif()
    endif()
    if(OpenMP_CXX_ENABLED)
      message(STATUS "C++ OpenMP enabled.")
      target_link_libraries(StormRuler_BASE INTERFACE OpenMP::OpenMP_CXX)
    else()
      message(STATUS "C++ OpenMP disabled, different Fortran and C++ compiler vendors.")
    endif()
  endif()
endif()

# Find MKL.
if(StormRuler_NO_MKL)
  message(STATUS "StormRuler MKL was disabled by user.")
else()
  # Configure OpenMP library to the Fortran compiler native one.
  if(OpenMP_Fortran_ENABLED)
    if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "INTEL")
      set(MKL_THREADING "intel_thread")
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
      set(MKL_THREADING "gnu_thread")
    elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI")
      set(MKL_THREADING "pgi_thread")
    endif()
  endif()
  # Configure MKL to use 32-bit integers.
  set(MKL_INTERFACE lp64)
  find_package(MKL CONFIG)
  if(MKL_FOUND)
    set(MKL_ENABLED 1)
    message(STATUS "MKL enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE MKL::MKL)
  endif()
endif()

# TODO:
# Find BLAS.
# Find Lapack.

# Find METIS.
if(StormRuler_NO_METIS)
  message(STATUS "StormRuler METIS was disabled by user.")
else()
  find_package(METIS)
  if(METIS_FOUND)
    message(STATUS "METIS enabled.")
    set(METIS_ENABLED 1)
    target_link_libraries(StormRuler_BASE INTERFACE METIS)
  endif()
endif()

# Find ZLIB.
if(StormRuler_NO_ZLIB)
  message(STATUS "StormRuler ZLIB was disabled by user.")
else()
  find_package(ZLIB)
  if(ZLIB_FOUND)
    message(STATUS "ZLIB enabled.")
    set(ZLIB_ENABLED 1)
    target_link_libraries(StormRuler_BASE INTERFACE ZLIB::ZLIB)
  endif()
endif()

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Main shared library and the executable.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

message(STATUS "Adding StormRuler shared library.")
add_subdirectory(
  "${CMAKE_SOURCE_DIR}/src")

message(STATUS "Adding StormRuler main executable.")
add_subdirectory(
  "${CMAKE_SOURCE_DIR}/src_aux/StormRuler_Main")

message(STATUS "Adding StormRuler LBM executable.")
add_subdirectory(
  "${CMAKE_SOURCE_DIR}/src_aux/StormRuler_LbMain")
