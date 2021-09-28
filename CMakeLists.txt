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
  LANGUAGES Fortran C CXX)

# A container library that stores all the properties.
add_library(StormRuler_BASE INTERFACE)

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Compilation options.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")

if(APPLE)
  # Detect Apple Silicon CPU.
  # Fortran compilers on Apple Silicon are not production-ready yet,
  # so we need to disable some features (for now).
  include(CheckCxxCompilerFlag)
  check_cxx_compiler_flag("-arch arm64" APPLE_SILICON)
  message(STATUS "Apple Silicon detected, this is bad..")
endif()

# NAG (unexpectedly) does not use byte-based kind system by default.
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
  target_compile_options(
    StormRuler_BASE INTERFACE 
    $<$<COMPILE_LANGUAGE:Fortran>:-kind=byte>)
endif()

#target_compile_options(your_target
#  PRIVATE
#  $<$<COMPILE_LANGUAGE:CXX>:$<$<CONFIG:Debug>:${CMAKE_CXX_FLAGS_RELEASE}>>
#  $<$<COMPILE_LANGUAGE:C>:$<$<CONFIG:Debug>:${CMAKE_C_FLAGS_RELEASE}>>

# TODO: add compile options based on Debug/Release configurations.
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
  target_compile_options(
    StormRuler_BASE INTERFACE $<$<COMPILE_LANGUAGE:Fortran>:-g -fcheck=all>)
  # Enabled optimization fail on GFORTRAN/M1.
  if(NOT APPLE_SILICON)
    target_compile_options(
      StormRuler_BASE INTERFACE $<$<COMPILE_LANGUAGE:Fortran>:-Ofast -march=native>)
  endif()
endif()

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Libraries.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

# Find Git and update submodules.
find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    RESULT_VARIABLE GIT_SUBMODULE_RESULT)
  if(NOT GIT_SUBMODULE_RESULT EQUAL "0")
    message(FATAL_ERROR 
      "`git submodule update --init --recursive` failed with ${GIT_SUBMODULE_RESULT}.")
  endif()
endif()

# TODO: we should check for quality of the OpenMP implementation,
# e.g. support of block constructs inside of the parallel loops.

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
      # which can should theoretically match used C++ compiler.
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
    # Intel suggest the other more explicit way, see: https://intel.ly/3ABZNsK
    #target_compile_options(
    #  StormRuler_BASE INTERFACE 
    #    $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    #target_include_directories(
    #  StormRuler_BASE INTERFACE 
    #    $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    #target_link_libraries(
    #  StormRuler_BASE INTERFACE $<LINK_ONLY:MKL::MKL>)
  endif()
endif()

# Find LAPACK.
if(StormRuler_NO_LAPACK)
  message(STATUS "StormRuler LAPACK was disabled by user.")
elseif(MKL_ENABLED)
  # Separate LAPACK is not needed in case MKL is used.
  set(LAPACK_ENABLED 1)
  message(STATUS "LAPACK enabled (MKL).")
else()
  find_package(LAPACK)
  if(LAPACK_FOUND)
    set(LAPACK_ENABLED 1)
    message(STATUS "LAPACK enabled.")
    target_link_libraries(StormRuler_BASE INTERFACE LAPACK::LAPACK)
  endif()
endif()

# Find Lua.
if(StormRuler_NO_Lua)
  message(STATUS "StormRuler Lua was disabled by user.")
else()
  find_package(Lua)
  if(LUA_FOUND)
    message(STATUS "Lua enabled.")
    set(Lua_ENABLED 1)
    # Lua target is missing.
    #target_link_libraries(StormRuler_BASE INTERFACE Lua::Lua)
    target_link_libraries(
      StormRuler_BASE INTERFACE "${LUA_LIBRARIES}")
    target_include_directories(
      StormRuler_BASE INTERFACE "${LUA_INCLUDE_DIR}")
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

if(Lua_ENABLED)
  message(STATUS "Adding StormRuler Lua executable.")
  add_subdirectory(
    "${CMAKE_SOURCE_DIR}/src_aux/StormRuler_Main_Lua")
endif()
