#!/usr/bin/env bash

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

# Setup C++ compiler and standard.
CXX=${CXX:-g++-12}
CXX_STD=${CXX_STD:-20}
echo "C++$CXX_STD compiler: $CXX"

# Setup the configuration.
CONFIG=${1:-Release}
echo "Configuration: $CONFIG"

# Setup vcpck root directory.
VCPKG_ROOT=${VCPKG_ROOT:-${VCPKG_INSTALLATION_ROOT:-$HOME/vcpkg}}
echo "vcpkg root directory: $VCPKG_ROOT"
VCPKG_TOOLCHAIN_FILE=$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake
if [ ! -f "$VCPKG_TOOLCHAIN_FILE" ]; then
    echo "Broken vcpkg installation!"
    exit 1
fi

# ------------------------------------------------------------------------------

# Remove old build directory.
rm -rf ./bin
rm -rf ./docs
rm -rf ./build

# ------------------------------------------------------------------------------

# Configure CMake.
export CXX
cmake -S . \
      -B build \
      -DCMAKE_BUILD_TYPE="$CONFIG" \
      -DCMAKE_CXX_STANDARD="$CXX_STD" \
      -DCMAKE_TOOLCHAIN_FILE="$VCPKG_TOOLCHAIN_FILE"

# ------------------------------------------------------------------------------
