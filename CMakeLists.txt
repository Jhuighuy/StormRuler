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

file(GLOB SOURCES_Fortran "src/*.f90") 
file(GLOB SOURCES_CXX "src/*.cpp")

add_executable(StormRuler ${SOURCES_Fortran} ${SOURCES_CXX})

# Require C++17.
set_target_properties(
  StormRuler PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED YES)

# Main is located in C++, so let C++ linker work.
set_target_properties(
  StormRuler PROPERTIES LINKER_LANGUAGE CXX)

# Specify compilation output directories.
set_target_properties(
  StormRuler PROPERTIES 
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib"
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/lib")

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Libraries and preprocessing.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

set(Fortiel_Defines "")

# Find OpenMP.
find_package(OpenMP)
if(OpenMP_Fortran_FOUND)
  message(STATUS "Fortran OpenMP enabled.")
  set(Fortiel_Defines "${Fortiel_Defines} -DHAS_OPENMP=True")
  target_link_libraries(StormRuler PUBLIC OpenMP::OpenMP_Fortran)
endif()
if(OpenMP_CXX_FOUND)
  # Enable C++ OpenMP only in case Fortran and C++ compiler vendors match.
  # Mixing different OpenMP runtimes is generally not a good idea.
  # TODO: here may be some exceptions, e.g. NAG uses system C compiler as a backend,
  # which can theoretically match used C++ compiler.
  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "${CMAKE_CXX_COMPILER_ID}")
    message(STATUS "C++ OpenMP enabled.")
    target_link_libraries(StormRuler PUBLIC OpenMP::OpenMP_CXX)
  else()
    message(STATUS "C++ OpenMP disabled, different Fortran and C++ compiler vendors.")
  endif()
endif()

# Find MKL.
# Configure OpenMP library to the Fortran compiler native one.
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
  target_link_libraries(StormRuler PRIVATE MKL::MKL)
  # Intel suggest the other, more explicit, way, see: https://intel.ly/3ABZNsK
  #target_compile_options(
  #  StormRuler PRIVATE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
  #target_include_directories(
  #  StormRuler PRIVATE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
  #target_link_libraries(
  #  StormRuler PRIVATE $<LINK_ONLY:MKL::MKL>)
endif()

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

# Find Python interpreter.
find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)

# Find Fortiel preprocessor.
set(Fortiel_PATH "${CMAKE_CURRENT_SOURCE_DIR}/ext/fortiel/src")
if(NOT EXISTS "${Fortiel_PATH}/fortiel.py")
  message(FATAL_ERROR "Broken Fortiel submodule.")
endif()

# Preprocess Fortran sources.
file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src_gen")
foreach(SOURCE_Fortran IN LISTS SOURCES_Fortran)
  file(RELATIVE_PATH SOURCE_Fortran_PREPROCESSED
       "${CMAKE_CURRENT_SOURCE_DIR}/src" "${SOURCE_Fortran}")
  get_filename_component(
    SOURCE_Fortran_PREPROCESSED "${SOURCE_Fortran_PREPROCESSED}" NAME_WE)
  set(SOURCE_Fortran_PREPROCESSED 
      "${CMAKE_CURRENT_SOURCE_DIR}/src_gen/${SOURCE_Fortran_PREPROCESSED}.out.f90")
  add_custom_command(
    OUTPUT "${SOURCE_Fortran_PREPROCESSED}"
    COMMAND "${Python3_EXECUTABLE}" "${Fortiel_PATH}/fortiel.py" 
            "${SOURCE_Fortran}" -o "${SOURCE_Fortran_PREPROCESSED}"
    MAIN_DEPENDENCY "${SOURCE_Fortran}"
    VERBATIM)
  list(APPEND SOURCES_Fortran_PREPROCESSED "${SOURCE_Fortran_PREPROCESSED}")
endforeach()
target_sources(StormRuler PRIVATE ${SOURCES_Fortran_PREPROCESSED})

## -----------------------------------------------------------------  

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
