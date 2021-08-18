<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->
# MUDPACK: Multigrid Software for Elliptic Partial Differential Equations
<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->

**MUDPACK** was first released in 1990, and has remained at version 5.0.1 
the past five years. It is a collection of portable, mostly Fortran 77 
subprograms (the code also employs a few Fortran 90 extensions). 
Its purpose is the efficient solving of linear elliptic Partial 
Differential Equations (PDEs) — both separable and nonseparable —
using multigrid iteration. 

**MUDPACK** solvers can achieve parallel speedup 
via OpenMP directives in the 5.0.1 subroutines, but the user will need 
to activate the directives by providing appropriate OpenMP compiler 
options when the library and application are built. Speedup is dependent 
on problem size and characteristics of the shared multi-processor 
platform where the code is compiled and run.

<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->
# Custom modifications
<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->
**MUDPACK** was modified in order to utilize the double precision 
real arithmetics.

<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->
# Licence
<!--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-->

Copyright © 2004 the University Corporation for 
Atmospheric Research ("UCAR"). All rights reserved. 
Developed by NCAR's Computational and Information Systems 
Laboratory, UCAR.

Redistribution and use of the Software in source and binary forms, 
with or without modification, is permitted provided that 
the following conditions are met:

* Neither the names of NCAR's Computational and Information Systems 
  Laboratory, the University Corporation for Atmospheric Research, 
  nor the names of its sponsors or contributors may be used to endorse 
  or promote products derived from this Software without specific 
  prior written permission.
    
* Redistributions of source code must retain the above copyright 
  notices, this list of conditions, and the disclaimer below.
    
* Redistributions in binary form must reproduce the above 
  copyright notice, this list of conditions, and the disclaimer below 
  in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO THE WARRANTIES OF 
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY 
CLAIM, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR 
THE USE OR OTHER DEALINGS WITH THE SOFTWARE.
