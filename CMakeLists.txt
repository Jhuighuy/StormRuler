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
## NONINFRINGEMENT. IN NO EVENmT SHALL THE AUTHORS OR COPYRIGHT 
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
        DESCRIPTION "A very high order CFD solver"
        LANGUAGES Fortran CXX)

# A container library that stores all the properties.
add_library(StormRuler_BASE INTERFACE)

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Compilation options.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

# NAG (unexpectedly) uses some wierd notation for kinds by default.
# Switch to the normal system.
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
  target_compile_options(
    StormRuler_BASE INTERFACE $<$<COMPILE_LANGUAGE:Fortran>:-kind=byte>)
endif()

target_compile_options(
  StormRuler_BASE INTERFACE $<$<COMPILE_LANGUAGE:Fortran>:-Ofast -march=native>)

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Libraries.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

set(Fortiel_DEFINES "")

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
      set(Fortiel_DEFINES "${Fortiel_DEFINES} -DHAS_OPENMP=True")
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
# Configure OpenMP library to the Fortran compiler native one.
if(StormRuler_NO_MKL)
  message(STATUS "StormRuler MKL was disabled by user.")
else()
  if(OpenMP_Fortran_FOUND)
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
    # Intel suggest the other, more explicit, way, see: https://intel.ly/3ABZNsK
    #target_compile_options(
    #  StormRuler_BASE INTERFACE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    #target_include_directories(
    #  StormRuler_BASE INTERFACE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    #target_link_libraries(
    #  StormRuler_BASE INTERFACE $<LINK_ONLY:MKL::MKL>)
  endif()
endif()

# Build MUDPACK.
if(StormRuler_NO_MUDPACK)
  message(STATUS "StormRuler MUDPACK was disabled by user.")
else()
  set(MUDPACK_ENABLED 1)
  message(STATUS "MUDPACK enabled.")
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/ext/mudpack")
  target_link_libraries(StormRuler_BASE INTERFACE MUDPACK)
endif()

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Sources and preprocessing.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

# Locate the Fortran and C++ sources.
file(GLOB_RECURSE Fortran_SOURCES "src/*.f90") 
file(GLOB_RECURSE CXX_SOURCES "src/*.cpp")

# Find Python interpreter.
find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)

# Find Fortiel preprocessor.
set(Fortiel_PATH "${CMAKE_CURRENT_SOURCE_DIR}/ext/fortiel/src")
if(NOT EXISTS "${Fortiel_PATH}/fortiel.py")
  message(FATAL_ERROR "Broken Fortiel submodule.")
endif()

# Preprocess Fortran sources.
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src_gen")
foreach(Fortran_SOURCE IN LISTS Fortran_SOURCES)
  file(RELATIVE_PATH Fortran_SOURCE_PREPROCESSED
       "${CMAKE_CURRENT_SOURCE_DIR}/src" "${Fortran_SOURCE}")
  get_filename_component(
    Fortran_SOURCE_PREPROCESSED "${Fortran_SOURCE_PREPROCESSED}" NAME_WE)
  set(Fortran_SOURCE_PREPROCESSED 
      "${CMAKE_CURRENT_SOURCE_DIR}/src_gen/${Fortran_SOURCE_PREPROCESSED}.out.f90")
  add_custom_command(
    OUTPUT "${Fortran_SOURCE_PREPROCESSED}"
    COMMAND "${Python3_EXECUTABLE}" "${Fortiel_PATH}/fortiel.py" 
            "${Fortran_SOURCE}" -o "${Fortran_SOURCE_PREPROCESSED}"
    MAIN_DEPENDENCY "${Fortran_SOURCE}"
    VERBATIM)
  list(APPEND Fortran_SOURCES_PREPROCESSED "${Fortran_SOURCE_PREPROCESSED}")
endforeach()

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Main executable.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

add_executable(StormRuler 
  ${Fortran_SOURCES} ${Fortran_SOURCES_PREPROCESSED} ${CXX_SOURCES})

# Link with the base target.
target_link_libraries(StormRuler PUBLIC StormRuler_BASE)

# Specify compilation output directories.
set_target_properties(
  StormRuler PROPERTIES 
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")
# Bug in CMake or NAG?
if(NOT "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
set_target_properties(
  StormRuler PROPERTIES 
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
else()
  set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")
endif()

# Require C++17.
set_target_properties(
  StormRuler PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# Main is located in C++, so let C++ linker work.
set_target_properties(
  StormRuler PROPERTIES LINKER_LANGUAGE Fortran)
