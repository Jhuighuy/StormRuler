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
    loadlibrary('libStormRuler', 'src/StormRuler_API.h');
  end

  SR.InitMesh = @SR_InitMesh;

  SR.AllocR = @SR_AllocR;
  SR.AllocC = @SR_AllocC;
  SR.AllocS = @SR_AllocS;

end

%@ <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< @%
%@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> @%

function [mesh] = SR_InitMesh()
  mesh = calllib('libStormRuler', 'SR_InitMesh');
end

%@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> @%
%@ >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> @%

function [field] = SR_AllocR(mesh, numVars, rank)
  field = calllib('libStormRuler', 'SR_AllocR', mesh, numVars, rank);
end

function [field] = SR_AllocC(mesh, numVars, rank)
  field = calllib('libStormRuler', 'SR_AllocC', mesh, numVars, rank);
end

function [field] = SR_AllocS(mesh, numVars, rank)
  field = calllib('libStormRuler', 'SR_AllocS', mesh, numVars, rank);
end
