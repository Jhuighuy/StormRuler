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

#include "StormRuler_API.h"

#include <stdio.h>
#include <stdlib.h>

void SR_PrintPointer(const void* p) {
  fprintf(stdout, "PRINT_PTR=%p\n", p);
  fflush(stdout);
} // SR_PrintPointer

SR_MRCI_REQUEST_R SR_MRCI_LinSolveR(SR_tMesh mesh,
    SR_STRING method, SR_STRING precondMethod,
    SR_tFieldR x, SR_tFieldR b) {

  SR_sMRCI_RequestR* pRequest = 
    (SR_sMRCI_RequestR*)malloc(sizeof(*pRequest)); 

  pRequest->request = SR_RCI_LinSolveR(mesh, 
    method, precondMethod, x, b, &pRequest->Ay, &pRequest->y);

  return pRequest;
} // SR_MRCI_LinSolveR

SR_STRING SR_MRCI_ReadRequest(SR_MRCI_REQUEST_R pRequest) {
  return pRequest->request != NULL ? pRequest->request : "Done";
} // SR_MRCI_ReadRequest

SR_tFieldR SR_MRCI_ReadFields(SR_MRCI_REQUEST_R pRequest, SR_INTEGER index) {
  return (index == 1) ? pRequest->Ay : pRequest->y;
} // SR_MRCI_ReadRequest

void SR_MRCI_Free(SR_MRCI_REQUEST_R pRequest) {
  free(pRequest);
} // SR_MRCI_Free
