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
## StormRuler C executable.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

project(StormRuler_Main
  VERSION 1.0
  DESCRIPTION "StormRuler main executable"
  LANGUAGES CXX)

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Sources.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

# Locate the C++ sources.
file(GLOB_RECURSE CXX_SOURCES "./*.cpp")

## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##
## Main executable library.
## -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- ##

add_executable(StormRuler_Main ${CXX_SOURCES})

# Require C++20.
set_target_properties(
  StormRuler_Main PROPERTIES 
  CXX_STANDARD 20 CXX_STANDARD_REQUIRED YES)

# Link with the main library.
target_link_libraries(StormRuler_Main PUBLIC StormRuler)

# Specify binary output directory.
set_target_properties(
  StormRuler_Main PROPERTIES 
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
set_target_properties(
  StormRuler_Main PROPERTIES
  VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
