# Copyright (C) 2020 - 2023 Oleg Butakov
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

# -----------------------------------------------------------------------------

# Find Doxygen.
find_package(Doxygen REQUIRED dot)
if(NOT DOXYGEN_FOUND)
  message(STATUS "| Doxygen was not found! Skipping documentation.")
  return()
endif()

# Check if Doxygen is recent enough.
set(DOXYGEN_MIN_REQUIRED_VERSION 1.9.5)
if(DOXYGEN_VERSION VERSION_LESS DOXYGEN_MIN_REQUIRED_VERSION)
  message(WARNING "| Outdated Doxygen ${DOXYGEN_VERSION} found! Minimum "
                    "required version is ${DOXYGEN_MIN_REQUIRED_VERSION}.")
  message(STATUS "| Skipping documentation.")
  return()
endif()

# -----------------------------------------------------------------------------

# Configure the required documenation components.
# (No alphabetical index, and a source browser with comments.)
set(DOXYGEN_ALPHABETICAL_INDEX  NO )
set(DOXYGEN_SOURCE_BROWSER      YES)
set(DOXYGEN_STRIP_CODE_COMMENTS NO )
set(DOXYGEN_TAB_SIZE            2  )

# -----------------------------------------------------------------------------

# Configure the DOT utility.
# (Enable it, and specify the parameters and style.)
set(DOXYGEN_CLASS_GRAPH          YES)
set(DOXYGEN_HIDE_UNDOC_RELATIONS YES)
set(DOXYGEN_COLLABORATION_GRAPH  NO )
set(DOXYGEN_GROUP_GRAPHS         YES)
set(DOXYGEN_TEMPLATE_RELATIONS   NO )
set(DOXYGEN_INCLUDE_GRAPH        NO )
set(DOXYGEN_INCLUDED_BY_GRAPH    NO )
set(DOXYGEN_CALL_GRAPH           NO )
set(DOXYGEN_CALLER_GRAPH         NO )
set(DOXYGEN_GRAPHICAL_HIERARCHY  NO )
set(DOXYGEN_DIRECTORY_GRAPH      NO )
set(DOXYGEN_GENERATE_LEGEND      NO )

# -----------------------------------------------------------------------------

# Configure the HTML pages layout.
# (Disable the navigation panel, enable the 
#  side bar with a tree-view that is folded by default.)
set(DOXYGEN_HTML_TIMESTAMP         NO )
set(DOXYGEN_HTML_DYNAMIC_MENUS     YES)
set(DOXYGEN_HTML_DYNAMIC_SECTIONS  YES)
set(DOXYGEN_HTML_INDEX_NUM_ENTRIES 1  )
set(DOXYGEN_DISABLE_INDEX          YES)
set(DOXYGEN_GENERATE_TREEVIEW      YES)
set(DOXYGEN_FULL_SIDEBAR           YES)

# Configure the HTML style.
# TODO: we definitely need better documentation styling :)
set(DOXYGEN_HTML_COLORSTYLE DARK)

# Configure the HTML formulas rendering.
# (No LaTeX, render formulas with MathJax.)
set(DOXYGEN_USE_MATHJAX YES)

# Other HTML options.
set(DOXYGEN_ENUM_VALUES_PER_LINE 1 )
set(DOXYGEN_OBFUSCATE_EMAILS     NO)

# -----------------------------------------------------------------------------

# Generate the HTML documentation.
set(DOXYGEN_GENERATE_HTML YES)
set(DOXYGEN_HTML_OUTPUT   "${CMAKE_SOURCE_DIR}/bin/docs")

# -----------------------------------------------------------------------------

# Generate the Doxygen inputs list.
set(StormRuler_DOXYGEN_INPUTS)

# Add the core library headers.
message(STATUS "| Adding the core library documentation..")
list(APPEND StormRuler_DOXYGEN_INPUTS ${StormRuler_CXX_SOURCES})
list(APPEND StormRuler_DOXYGEN_INPUTS ${StormRuler_MD_SOURCES})

# Add the Doxygen target.
doxygen_add_docs(StormRuler_DOXYGEN ALL
                 USE_STAMP_FILE ${StormRuler_CXX_SOURCES}
                 COMMENT "Generate Doxygen documentation")
