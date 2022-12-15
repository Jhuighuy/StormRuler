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
# FITNESS FOR Allocator PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
# SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# ------------------------------------------------------------------------------

# Core header-only library.
file(GLOB_RECURSE StormRuler_CXX_SOURCES "Storm/*.hpp")
set(StormRuler_CXX_SOURCES ${StormRuler_CXX_SOURCES} PARENT_SCOPE)
add_library(StormRuler INTERFACE ${StormRuler_CXX_SOURCES})

# Include directories.
target_include_directories(StormRuler_BASE INTERFACE ".")

# ------------------------------------------------------------------------------

# Configuration file.
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/Storm/Config.hpp.in" 
    "${CMAKE_CURRENT_SOURCE_DIR}/Storm/Config.hpp")

# "Default" configuration file. 
# (It should be in the repository, used for non-CMake builds.)
add_custom_command(
    OUTPUT "./Storm/ConfigDefault.hpp"
    COMMAND ${Python3_EXECUTABLE} 
        "../scripts/MakeConfigDefault.py" 
        "${CMAKE_CURRENT_SOURCE_DIR}/Storm/Config.hpp.in" 
        "${CMAKE_CURRENT_SOURCE_DIR}/Storm/ConfigDefault.hpp"
    MAIN_DEPENDENCY "./Storm/Config.hpp.in"
    VERBATIM)

# ------------------------------------------------------------------------------

# Link with the base target.
target_link_libraries(StormRuler INTERFACE StormRuler_BASE)

# ------------------------------------------------------------------------------
