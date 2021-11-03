// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// Copyright (C) 2021 Oleg Butakov
// 
// Permission is hereby granted, free of charge, to any person 
// obtaining a copy of this software and associated documentation 
// files (the "Software"), to deal in the Software without 
// restriction, including without limitation the rights  to use, 
// copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the  
// Software is furnished to do so, subject to the following 
// conditions:
// 
// The above copyright notice and this permission notice shall be 
// included in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

#define SR_MATLAB 0
#include "StormRuler_API.h"

#include <stdio.h>
#include <stdlib.h>

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

void SR_PrintPointer(const void* p) {
  fprintf(stdout, "PRINT_PTR=%p\n", p);
  fflush(stdout);
} // SR_PrintPointer

// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< //
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> //

#if 0
struct SR_tRequestMRCIStruct {
  SR_STRING String;
  SR_tFieldR Ay, y;
}; // struct SR_tRequestMRCIStruct

SR_tRequestMRCI SR_MRCI_LinSolveR(SR_tMesh mesh,
    SR_STRING method, SR_STRING precondMethod,
    SR_tFieldR x, SR_tFieldR b) {
  static SR_tRequestMRCIStruct request;
  request.String = SR_RCI_LinSolveR(mesh, 
    method, precondMethod, x, b, &request.Ay, &request.y);
  return &request;
} // SR_MRCI_LinSolveR

SR_STRING SR_MRCI_ReadRequest(SR_tRequestMRCI request) {
  return request->String != NULL ? request->String : "Done";
} // SR_MRCI_ReadRequest

SR_tFieldR SR_MRCI_ReadFields(SR_tRequestMRCI request, SR_INTEGER index) {
  if (index == 1) return request->Ay;
  if (index == 2) return request->y;
  return NULL;
} // SR_MRCI_ReadFields
#endif
