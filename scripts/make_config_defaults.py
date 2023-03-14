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
import re

# ------------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate the default values for the CMake configuration file."
    )
    parser.add_argument("input_file_path", help="Input CMake configuration file path.")
    parser.add_argument("output_file_path", help="Generated output file path")
    args = parser.parse_args()

    FIND_REGEXP = r"#cmakedefine01\s*(\w+)$"
    REPLACE_REGEXP = r"#ifndef \1\n#  define \1 0\n#endif"

    with open(args.input_file_path, "r") as input_file:
        contents = input_file.read()
        generated = re.sub(FIND_REGEXP, REPLACE_REGEXP, contents, flags=re.M)

    with open(args.output_file_path, "w") as output_file:
        output_file.write(generated)


# ------------------------------------------------------------------------------
