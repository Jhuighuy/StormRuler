%@ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< @%
%@ Copyright (C) 2021 Oleg Butakov
%@ 
%@ Permission is hereby granted, free of charge, to any person 
%@ obtaining a copy of this software and associated documentation 
%@ files (the "Software"), to deal in the Software without 
%@ restriction, including without limitation the rights  to use, 
%@ copy, modify, merge, publish, distribute, sublicense, and/or
%@ sell copies of the Software, and to permit persons to whom the  
%@ Software is furnished to do so, subject to the following 
%@ conditions:
%@ 
%@ The above copyright notice and this permission notice shall be 
%@ included in all copies or substantial portions of the Software.
%@ 
%@ THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
%@ EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
%@ OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
%@ NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
%@ HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
%@ WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
%@ FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
%@ OTHER DEALINGS IN THE SOFTWARE.
%@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> @%

function SR = StormRuler_Matlab()

  if ~libisloaded('libStormRuler')
    addpath('bin', 'src');
    loadlibrary('libStormRuler.so', 'StormRuler_API.h');
  end

  libfunctions('libStormRuler','-full')
  
  clear SR;

  SR.InitMesh = @SR_InitMesh;
  SR.Mesh_NumCells = @SR_Mesh_NumCells;

  SR.AllocR = @SR_AllocR;
  SR.Alloc_Mold = @SR_Alloc_Mold;
  SR.Free = @SR_Free;
  SR.At = @SR_At;
  SR.SetAt = @SR_SetAt;

  SR.IO_Begin = @SR_IO_Begin;
  SR.IO_Add = @SR_IO_Add;
  SR.IO_Flush = @SR_IO_Flush;

  SR.Fill = @SR_Fill;
  SR.Fill_Random = @SR_Fill_Random;
  SR.Set = @SR_Set;
  SR.Scale = @SR_Scale;
  SR.Add = @SR_Add;
  SR.Sub = @SR_Sub;
  SR.Mul = @SR_Mul;

  SR.RCI_LinSolve = @SR_RCI_LinSolve;

  SR.ApplyBCs = @SR_ApplyBCs;
  SR.ApplyBCs_Dirichlet = @SR_ApplyBCs_Dirichlet;
  SR.ApplyBCs_Neumann = @SR_ApplyBCs_Neumann;
  SR.ApplyBCs_PureDirichlet = @SR_ApplyBCs_PureDirichlet;
  SR.ApplyBCs_PureNeumann = @SR_ApplyBCs_PureNeumann;
  SR.Div = @SR_Div;
  SR.Grad = @SR_Grad;
  SR.Conv = @SR_Conv;
  SR.DivGrad = @SR_DivGrad;
  SR.DivKGrad = @SR_DivKGrad;

  % Reset the RCI state.
  nullMesh = libpointer('SR_tMeshStructPtr');
  nullField = libpointer('SR_tFieldRStructPtr');
  calllib('libStormRuler', 'SR_MRCI_LinSolveR', ...
    nullMesh, '', '', nullField, nullField);

end

%@ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< @%
%@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> @%

function [mesh] = SR_InitMesh()
  mesh = calllib('libStormRuler', 'SR_InitMesh');
end

function [numCells] = SR_Mesh_NumCells(mesh)
  numCells = calllib('libStormRuler', 'SR_Mesh_NumCells', mesh);
end

%@ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< @%
%@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> @%

function [field] = SR_AllocR(mesh, numVars, rank)
  field = calllib('libStormRuler', 'SR_AllocR', mesh, numVars, rank);
end

function [field] = SR_Alloc_Mold(mold)
  field = calllib('libStormRuler', 'SR_Alloc_MoldR', mold);
end

function SR_Free(field)
  calllib('libStormRuler', 'SR_FreeR', field);
end

function value = SR_At(field, index)
  pointer = calllib('libStormRuler', 'SR_AtR', field, index);
  setdatatype(pointer,'doublePtr',1);
  value = pointer.Value;
end

function SR_SetAt(field, index, value)
  pointer = calllib('libStormRuler', 'SR_AtR', field, index);
  setdatatype(pointer,'doublePtr',1);
  pointer.Value = value;
end

%@ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< @%
%@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> @%

function [ioList] = SR_IO_Begin()
  ioList = calllib('libStormRuler', 'SR_IO_Begin');
end

function SR_IO_Add(ioList, x, name)
  calllib('libStormRuler', 'SR_IO_Add', ioList, x, name);
end

function SR_IO_Flush(ioList, mesh, filename)
  calllib('libStormRuler', 'SR_IO_Flush', ioList, mesh, filename);
end

%@ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< @%
%@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> @%

function SR_Fill(mesh, x, alpha, beta)
  calllib('libStormRuler', 'SR_FillR', mesh, x, alpha, beta);
end

function SR_Fill_Random(mesh, x, a, b)
  calllib('libStormRuler', 'SR_Fill_RandomR', mesh, x, a, b);
end

function SR_Set(mesh, y, x)
  calllib('libStormRuler', 'SR_SetR', mesh, y, x);
end

function SR_Scale(mesh, y, x, alpha)
  calllib('libStormRuler', 'SR_ScaleR', mesh, y, x, alpha);
end

function SR_Add(mesh, z, y, x, alpha, beta)
  calllib('libStormRuler', 'SR_AddR', mesh, z, y, x, alpha, beta);
end

function SR_Sub(mesh, z, y, x, alpha, beta)
  SR_Add(mesh, z, y, x, -alpha, beta)
end

function SR_Mul(mesh, z, y, x)
  calllib('libStormRuler', 'SR_MulR', mesh, z, y, x);
end

%@ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< @%
%@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> @%

function [request, Ay, y] = SR_RCI_LinSolve(mesh, method, precondMethod, x, b)
  pRequest = calllib( ...
    'libStormRuler', 'SR_MRCI_LinSolveR', mesh, method, precondMethod, x, b);
  
  request = calllib( ...
    'libStormRuler', 'SR_MRCI_ReadRequest', pRequest);
  Ay = calllib('libStormRuler', 'SR_MRCI_ReadFields', pRequest, 1);
  y = calllib('libStormRuler', 'SR_MRCI_ReadFields', pRequest, 2);
end

%@ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< @%
%@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> @%

function SR_ApplyBCs(mesh, u, iBC, alpha, beta, gamma)
  calllib('libStormRuler', 'SR_ApplyBCsR', ...
    mesh, u, iBC, alpha, beta, gamma);
end

function SR_ApplyBCs_Dirichlet(mesh, u, iBC, gamma)
  SR_ApplyBCs(mesh, u, iBC, 1.0, 0.0, gamma);
end

function SR_ApplyBCs_Neumann(mesh, u, iBC, gamma)
  SR_ApplyBCs(mesh, u, iBC, 0.0, 1.0, gamma);
end

function SR_ApplyBCs_PureDirichlet(mesh, u, iBC)
  SR_ApplyBCs_Dirichlet(mesh, u, iBC, 0.0);
end

function SR_ApplyBCs_PureNeumann(mesh, u, iBC)
  SR_ApplyBCs_Neumann(mesh, u, iBC, 0.0);
end

function SR_Grad(mesh, vVec, lambda, u)
  calllib('libStormRuler', 'SR_GradR', mesh, vVec, lambda, u);
end

function SR_Div(mesh, v, lambda, uVec)
  calllib('libStormRuler', 'SR_DivR', mesh, v, lambda, uVec);
end

function SR_Conv(mesh, v, lambda, u, a)
  calllib('libStormRuler', 'SR_ConvR', mesh, v, lambda, u, a);
end

function SR_DivGrad(mesh, v, lambda, u)
  calllib('libStormRuler', 'SR_DivGradR', mesh, v, lambda, u);
end

function SR_DivKGrad(mesh, v, lambda, k, u)
  calllib('libStormRuler', 'SR_DivKGradR', mesh, v, lambda, k, u);
end
