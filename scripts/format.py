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
import os
import subprocess

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Count lines of code of the repository files."
    )
    parser.add_argument(
        "directory",
        metavar="DIR",
        action="store",
        nargs="?",
        default=None,
        help="Directory to format files in",
    )
    args = parser.parse_args()

    # Get the indexed files.
    git_ls_files_args = ["git", "ls-files"]
    if args.directory is not None:
        git_ls_files_args.append(args.directory)
    indexed_files = subprocess.check_output(git_ls_files_args).decode().splitlines()

    # Format the files!
    for indexed_file in indexed_files:
        _, extension = os.path.splitext(indexed_file)
        if extension == ".py":
            print(f"Format Python file {indexed_file} with Black..")
            black_args = ["black", "-q", indexed_file]
            results = subprocess.check_call(black_args)
        elif extension in [".h", ".c", ".hpp", ".cpp"]:
            print(f"Format C++ file {indexed_file} wih clang-format..")
            clang_format_args = ["clang-format", "-i", indexed_file]
            results = subprocess.check_call(clang_format_args)

# ------------------------------------------------------------------------------
