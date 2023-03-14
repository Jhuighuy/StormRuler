#!/usr/bin/env python3

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

import argparse
import os
import shutil
import subprocess

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Configure StormRuler with CMake.",
    )
    parser.add_argument(
        "-vcpkg",
        "--vcpkg_root",
        metavar="PATH",
        action="store",
        default=None,
        help="vcpkg installation root path",
    )
    parser.add_argument(
        "-cxx",
        "--cxx_compiler",
        metavar="EXE",
        action="store",
        default=None,
        help="C++ compiler executable",
    )
    parser.add_argument(
        "-std",
        "--cxx_standard",
        metavar="STD",
        action="store",
        default=None,
        choices=["20", "23", "c++20", "c++23"],
        help="C++ standard",
    )
    parser.add_argument(
        "-cfg",
        "--configuration",
        metavar="CONFIG",
        action="store",
        default=None,
        choices=["Debug", "Release", "Coverage"],
        help="build configuration",
    )
    parser.add_argument(
        "-args",
        "--arguments",
        metavar="ARGS",
        action="store",
        nargs=argparse.REMAINDER,
        default=None,
        help="extra CMake arguments",
    )
    args = parser.parse_args()

    # Remove the old build directories.
    for output_dir in ["bin", "docs", "build"]:
        shutil.rmtree(output_dir, ignore_errors=True)

    # Prepare the CMake arguments.
    cmake_args = ["cmake"]

    # Setup the source and build directories.
    cmake_args += ["-S", ".", "-B", "build"]

    # Setup the vcpkg root.
    vcpkg_root = (
        args.vcpkg_root
        or os.environ.get("VCPKG_ROOT")
        or os.environ.get("VCPKG_INSTALLATION_ROOT")
    )
    if vcpkg_root is None:
        vcpkg_root_candidate = os.path.join(os.path.expanduser("~"), "vcpkg")
        if os.path.exists(vcpkg_root_candidate) and os.path.isdir(vcpkg_root_candidate):
            vcpkg_root = vcpkg_root_candidate
    if vcpkg_root is not None:
        vcpkg_toolchain_file = os.path.join(
            vcpkg_root, "scripts", "buildsystems", "vcpkg.cmake"
        )
        cmake_args.append(f"-DCMAKE_TOOLCHAIN_FILE={vcpkg_toolchain_file}")

    # Setup the C++ configuration.
    if args.cxx_compiler is not None:
        cxx_compiler = args.cxx_compiler
        cmake_args.append(f"-DCMAKE_CXX_COMPILER={cxx_compiler}")
    if args.cxx_standard is not None:
        cxx_standard = args.cxx_standard.removeprefix("c++")
        cmake_args.append(f"-DCMAKE_CXX_STANDARD={cxx_standard}")

    # Setup the build configuration.
    if args.configuration is not None:
        cmake_args.append(f"-DCMAKE_BUILD_TYPE={args.configuration}")

    # Append the extra arguments.
    if args.arguments is not None:
        cmake_args += args.arguments

    # Run CMake!
    subprocess.check_call(cmake_args)

# ------------------------------------------------------------------------------
