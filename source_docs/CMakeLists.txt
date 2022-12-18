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

# Configure so that `cmake .` can be ran directly from this folder.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  cmake_minimum_required(VERSION 3.20..3.26)
  project(StormRuler
          VERSION 0.1
          DESCRIPTION "A very high order FVM framework" 
          LANGUAGES NONE)
endif()

# ------------------------------------------------------------------------------

# Find Doxygen.
find_package(Doxygen QUIET)
if(NOT DOXYGEN_FOUND)
  return()
endif()

# Check if Doxygen is recent enough.
if(DOXYGEN_VERSION VERSION_LESS "1.9.4")
  message(STATUS "| Doxygen 1.9.4 or more recent is required!")
  return()
endif()

# ------------------------------------------------------------------------------

# Configure the documentation extraction.
set(DOXYGEN_DISTRIBUTE_GROUP_DOC      YES)
set(DOXYGEN_PREDEFINED "STORM_DOXYGEN_=1")
set(DOXYGEN_SHOW_USED_FILES           NO )

# ------------------------------------------------------------------------------

# Configure the required documenation components.
# (No alphabetical index, and a source browser with comments.)
set(DOXYGEN_ALPHABETICAL_INDEX  NO )
set(DOXYGEN_SOURCE_BROWSER      YES)
set(DOXYGEN_STRIP_CODE_COMMENTS NO )

# ------------------------------------------------------------------------------

# Configure graphs.
# (Enable it, and specify the parameters and style.)
set(DOXYGEN_CALLER_GRAPH         NO )
set(DOXYGEN_CALL_GRAPH           NO )
set(DOXYGEN_COLLABORATION_GRAPH  NO )
set(DOXYGEN_DIRECTORY_GRAPH      NO )
set(DOXYGEN_GENERATE_LEGEND      NO )
set(DOXYGEN_GRAPHICAL_HIERARCHY  NO )
set(DOXYGEN_INCLUDED_BY_GRAPH    NO )
set(DOXYGEN_INCLUDE_GRAPH        NO )
set(DOXYGEN_TEMPLATE_RELATIONS   NO )

# ------------------------------------------------------------------------------

# Generate the HTML documentation.
set(DOXYGEN_GENERATE_HTML YES)
set(DOXYGEN_HTML_OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/../docs)

# Configure the HTML pages layout.
# (Disable the navigation panel, enable the 
#  side bar with a tree-view that is folded by default.)
set(DOXYGEN_DISABLE_INDEX          YES)
set(DOXYGEN_FULL_SIDEBAR           YES)
set(DOXYGEN_GENERATE_TREEVIEW      YES)
set(DOXYGEN_HTML_DYNAMIC_MENUS     YES)
set(DOXYGEN_HTML_DYNAMIC_SECTIONS  YES)
set(DOXYGEN_HTML_INDEX_NUM_ENTRIES 1  )
set(DOXYGEN_HTML_TIMESTAMP         NO )

# Other HTML options.
set(DOXYGEN_ENUM_VALUES_PER_LINE 1   )
set(DOXYGEN_HTML_COLORSTYLE      DARK)
set(DOXYGEN_OBFUSCATE_EMAILS     NO  )
set(DOXYGEN_USE_MATHJAX          YES )

# ------------------------------------------------------------------------------

# Add the C++ sources and Markdown pages.
file(GLOB_RECURSE MD_SOURCES ./*.md)
list(APPEND StormRuler_Doxygen_INPUTS ${MD_SOURCES})
file(GLOB_RECURSE CXX_SOURCES ../source/*.hpp)
list(APPEND StormRuler_Doxygen_INPUTS ${CXX_SOURCES})

# Add the Doxygen target.
message(STATUS "| Adding the Doxygen documentation..")
doxygen_add_docs(StormRuler_Doxygen ALL
                 USE_STAMP_FILE ${StormRuler_Doxygen_INPUTS}
                 COMMENT "Generate Doxygen documentation")

# ------------------------------------------------------------------------------
