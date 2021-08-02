## <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
## >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## -----------------------------------------------------------------  
## StormRuler CMake configuration.
cmake_minimum_required(VERSION 3.5.1)
## -----------------------------------------------------------------  

## -----------------------------------------------------------------  
## StormRuler project definitions.
project(StormRuler)
enable_language(Fortran CXX)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/bin)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/lib)
## -----------------------------------------------------------------  

## -----------------------------------------------------------------  
## StormRuler executable.
#set(CMAKE_Fortran_FLAGS "-g -O0 -fbounds-check")
set(CMAKE_Fortran_FLAGS "-g -Ofast -fopenmp -march=native -Wall")
file(GLOB SOURCES src/*.f90 src/*.cpp) 
add_executable("StormRuler.out" ${SOURCES})
set_target_properties("StormRuler.out" PROPERTIES LINKER_LANGUAGE Fortran)
## -----------------------------------------------------------------  
