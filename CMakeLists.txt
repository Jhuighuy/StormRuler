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
set(CMAKE_CXX_STANDARD 17)

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Setup StormRuler project, using Fortran and C++.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

project(StormRuler
        VERSION 1.0
        DESCRIPTION "A very high order CFD solver"
        LANGUAGES Fortran CXX)

## Gather Fortran and C++ sources.
file(GLOB fortran_source_list "src/*.f90") 
file(GLOB cxx_source_list "src/*.cpp")

add_executable(StormRuler ${fortran_sources} ${cxx_source_list})

set_target_properties(StormRuler PROPERTIES LINKER_LANGUAGE Fortran)

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Set binary output directory.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY 
  "${CMAKE_CURRENT_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY 
  "${CMAKE_CURRENT_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY 
  "${CMAKE_CURRENT_BINARY_DIR}/lib")
set(CMAKE_Fortran_MODULE_DIRECTORY 
  "${CMAKE_CURRENT_BINARY_DIR}/lib")

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Libraries and preprocessing.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

set(fortiel_defines "")

# Find OpenMP.
find_package(OpenMP)
if(OpenMP_Fortran_FOUND)
  #set(fortiel_defines "${fortiel_defines} -DHAS_OPENMP=True")
  target_link_libraries(StormRuler PUBLIC OpenMP::OpenMP_Fortran)
endif()
if(OpenMP_CXX_FOUND)
  target_link_libraries(StormRuler PUBLIC OpenMP::OpenMP_CXX)
endif()

# Find MKL.
set(MKL_INTERFACE "lp64")
if(OpenMP_Fortran_FOUND)
  if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "INTEL")
    set(MKL_THREADING "intel_thread")
  elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set(MKL_THREADING "gnu_thread")
  elseif("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "PGI")
    set(MKL_THREADING "pgi_thread")
  endif()
endif()
find_package(MKL)
if(MKL_FOUND)
  #set(fortiel_defines "${fortiel_defines} -DHAS_MKL=True")
  target_link_libraries(StormRuler PUBLIC MKL::MKL)
endif()

# Find Python interpreter.
find_package(Python3 3.9 REQUIRED COMPONENTS Interpreter)

# Find Fortiel preprocessor.
if(NOT DEFINED ENV{FORTIEL_PATH})
  message(FATAL_ERROR 
          "Environmental variable $FORTIEL_PATH is not set.")
endif()
set(fortiel_path "$ENV{FORTIEL_PATH}/fortiel.py")
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/src_gen")

# Preprocess Fortran sources.
foreach(fortran_source IN LISTS fortran_source_list)
  file(RELATIVE_PATH fortran_preprocessed_source
       "${CMAKE_CURRENT_BINARY_DIR}/src" "${fortran_source}")
  set(fortran_preprocessed_source 
      "${CMAKE_CURRENT_BINARY_DIR}/src_gen/${fortran_preprocessed_source}")
  add_custom_command(
    OUTPUT "${fortran_preprocessed_source}"
    COMMAND "${Python3_EXECUTABLE}" "${fortiel_path}" 
            "${fortran_source}" -o "${fortran_preprocessed_source}"
    MAIN_DEPENDENCY "${fortran_source}"
    VERBATIM)
  list(APPEND fortran_preprocessed_source_list "${fortran_preprocessed_source}")
endforeach()
target_sources(StormRuler PRIVATE ${fortran_preprocessed_source_list})

## -----------------------------------------------------------------  
