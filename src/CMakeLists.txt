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

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## StormRuler shared library.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

project(StormRuler_Lib
  VERSION 1.0
  DESCRIPTION "StormRuler shared library"
  LANGUAGES Fortran C)

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Sources and preprocessing.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

# Locate the Fortran and C++ sources.
file(GLOB Fortran_SOURCES "./*.f90") 
file(GLOB C_SOURCES "./*.c")

# Build a list of Fortiel options.
set(Fortiel_DEFINES "")
if(OpenMP_Fortran_ENABLED)
  list(APPEND Fortiel_DEFINES "-DHAS_OpenMP=True")
endif()
if(MKL_ENABLED)
  list(APPEND Fortiel_DEFINES "-DHAS_MKL=True")
endif()
if(LAPACK_ENABLED)
  list(APPEND Fortiel_DEFINES "-DHAS_LAPACK=True")
endif()

# Find Fortiel preprocessor.
if(NOT "$ENV{StormRuler_Fortiel_PATH}" STREQUAL "")
  set(Fortiel_PATH "$ENV{StormRuler_Fortiel_PATH}")
  message(STATUS "Overriding submodule Fortiel path to ${Fortiel_PATH}")
else()
  set(Fortiel_PATH "${CMAKE_SOURCE_DIR}/ext/fortiel/src")
endif()
if(NOT EXISTS "${Fortiel_PATH}/fortiel.py")
  message(FATAL_ERROR "Broken Fortiel submodule.")
endif()

# Find Python interpreter.
find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)

# Preprocess Fortran sources.
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/gen")
foreach(Fortran_SOURCE IN LISTS Fortran_SOURCES)
  get_filename_component(
    Fortran_SOURCE_PREPROCESSED "${Fortran_SOURCE}" NAME_WE)
  set(Fortran_SOURCE_PREPROCESSED 
      "${CMAKE_CURRENT_SOURCE_DIR}/gen/${Fortran_SOURCE_PREPROCESSED}.out.f90")
  add_custom_command(
    OUTPUT "${Fortran_SOURCE_PREPROCESSED}"
    COMMAND "${Python3_EXECUTABLE}" "${Fortiel_PATH}/fortiel.py"
            ${Fortiel_DEFINES} -o "${Fortran_SOURCE_PREPROCESSED}" "${Fortran_SOURCE}"
    MAIN_DEPENDENCY "${Fortran_SOURCE}"
    VERBATIM)
  list(APPEND Fortran_SOURCES_PREPROCESSED "${Fortran_SOURCE_PREPROCESSED}")
endforeach()

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Main shared library.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

add_library(StormRuler SHARED ${Fortran_SOURCES}
  ${Fortran_SOURCES_PREPROCESSED} ${C_SOURCES})

# Link with the base target.
target_link_libraries(StormRuler PUBLIC StormRuler_BASE)

# Require C11.
set_target_properties(
  StormRuler PROPERTIES 
  C_STANDARD 11 C_STANDARD_REQUIRED YES)

# Let Fortran linker work.
set_target_properties(
  StormRuler PROPERTIES LINKER_LANGUAGE Fortran)

# Specify compilation output directories:
# binary output directory,
set_target_properties(
  StormRuler PROPERTIES 
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set_target_properties(
  StormRuler PROPERTIES 
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
# module output directiry,
set_target_properties(
  StormRuler PROPERTIES 
  Fortran_MODULE_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "NAG")
  # Bug in CMake or NAG?
  set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
endif()
