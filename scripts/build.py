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
# FITNESS FOR Allocator PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
# SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# ------------------------------------------------------------------------------

import argparse
import multiprocessing
import os
import re
import subprocess

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build StormRuler with CMake.",
    )
    parser.add_argument(
        "-j",
        "--jobs",
        metavar="N",
        action="store",
        nargs="?",
        type=int,
        const=0,
        default=None,
        help="number of threads to parallelize the build",
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
    args = parser.parse_args()

    # Prepare the CMake arguments.
    cmake_args = ["cmake", "--build", "build"]

    # Setup the build configuration.
    if args.configuration is not None:
        cmake_args += ["--config", args.configuration]

    # Setup the number of threads.
    if args.jobs is not None:
        # (`-j` is a valid option for Makefile or Ninja generators only.)
        def jobs_supported():
            try:
                cmake_cache_path = os.path.join("build", "CMakeCache.txt")
                with open(cmake_cache_path, "r") as cmake_cache_file:
                    cmake_cache = cmake_cache_file.read()
                cmake_generator = re.search(
                    r"^CMAKE_GENERATOR:INTERNAL=(.+)$", cmake_cache, re.MULTILINE
                )[1]
                return "Ninja" in cmake_generator or "Makefiles" in cmake_generator
            except (IOError, TypeError):
                return False

        if jobs_supported():
            num_threads = args.jobs or multiprocessing.cpu_count()
            cmake_args += ["--", "-j", str(num_threads)]

    # Run CMake!
    subprocess.check_call(cmake_args)

# ------------------------------------------------------------------------------
