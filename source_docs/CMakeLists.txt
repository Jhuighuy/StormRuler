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
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# ------------------------------------------------------------------------------

# This may be a standalone source directory. Because of this, consider using
# paths relative to `CMAKE_CURRENT_SOURCE_DIR` instead of `CMAKE_SOURCE_DIR`.
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

# Doxygen in too verbose by default.
set(DOXYGEN_QUIET YES)

# Configure the project logo.
set(DOXYGEN_PROJECT_LOGO ${CMAKE_CURRENT_SOURCE_DIR}/Images/LogoIcon.svg)

# Configure the documentation extraction.
set(DOXYGEN_DISTRIBUTE_GROUP_DOC   YES )
set(DOXYGEN_PREDEFINED "STORM_DOXYGEN_=1")
set(DOXYGEN_SHOW_USED_FILES        NO  )
set(DOXYGEN_STRIP_CODE_COMMENTS    NO  )
set(DOXYGEN_STRIP_FROM_INC_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../source)

# Configure the required documenation components.
set(DOXYGEN_ALPHABETICAL_INDEX     NO  )
set(DOXYGEN_SOURCE_BROWSER         YES )

# Configure the HTML output.
set(DOXYGEN_HTML_OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/../docs)

# Configure the HTML pages layout.
set(DOXYGEN_DISABLE_INDEX          NO  )
set(DOXYGEN_FULL_SIDEBAR           NO  )
set(DOXYGEN_GENERATE_TREEVIEW      YES )

# Configure the HTML pages components.
set(DOXYGEN_ENUM_VALUES_PER_LINE   1   )
set(DOXYGEN_HTML_COLORSTYLE        DARK)
set(DOXYGEN_HTML_DYNAMIC_MENUS     YES )
set(DOXYGEN_HTML_DYNAMIC_SECTIONS  YES )
set(DOXYGEN_HTML_INDEX_NUM_ENTRIES 1   )
set(DOXYGEN_HTML_TIMESTAMP         NO  )
set(DOXYGEN_OBFUSCATE_EMAILS       NO  )
set(DOXYGEN_USE_MATHJAX            YES )

# Configure graphs.
set(DOXYGEN_DOT_IMAGE_FORMAT       svg )
set(DOXYGEN_DOT_TRANSPARENT        YES )
set(DOXYGEN_CALLER_GRAPH           NO  )
set(DOXYGEN_CALL_GRAPH             NO  )
set(DOXYGEN_COLLABORATION_GRAPH    NO  )
set(DOXYGEN_DIRECTORY_GRAPH        NO  )
set(DOXYGEN_GENERATE_LEGEND        NO  )
set(DOXYGEN_GRAPHICAL_HIERARCHY    NO  )
set(DOXYGEN_INCLUDED_BY_GRAPH      NO  )
set(DOXYGEN_INCLUDE_GRAPH          NO  )
set(DOXYGEN_TEMPLATE_RELATIONS     NO  )

# ------------------------------------------------------------------------------

# Fetch Doxygen Awesome CSS.
include(FetchContent)
FetchContent_Declare(
    doxygen_awesome_css
    GIT_REPOSITORY https://github.com/jothepro/doxygen-awesome-css
    GIT_TAG        a3c119b4797be2039761ec1fa0731f038e3026f6) # release-2.1.0
FetchContent_MakeAvailable(doxygen_awesome_css)

# Enable Doxygen Awesome theme.
set(DOXYGEN_HTML_HEADER PageHeader.html.in)
set(DOXYGEN_HTML_FOOTER PageFooter.html.in)
set(DOXYGEN_HTML_EXTRA_FILES
    ${doxygen_awesome_css_SOURCE_DIR}/doxygen-awesome-darkmode-toggle.js
    ${doxygen_awesome_css_SOURCE_DIR}/doxygen-awesome-fragment-copy-button.js
    ${doxygen_awesome_css_SOURCE_DIR}/doxygen-awesome-paragraph-link.js
    ${doxygen_awesome_css_SOURCE_DIR}/doxygen-awesome-interactive-toc.js)
set(DOXYGEN_HTML_EXTRA_STYLESHEET
    ${doxygen_awesome_css_SOURCE_DIR}/doxygen-awesome.css
    ${doxygen_awesome_css_SOURCE_DIR}/doxygen-awesome-sidebar-only.css
    ${doxygen_awesome_css_SOURCE_DIR}/doxygen-awesome-sidebar-only-darkmode-toggle.css)

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
