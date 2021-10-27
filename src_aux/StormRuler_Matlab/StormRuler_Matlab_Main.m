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

tau = (pi/50)*(pi/50); 

sigma = 1.0;
gamma = 0.01;

rho = 1.0;
mu = 0.1;

SR = StormRuler_Matlab();

system('rm out/*');

mesh = SR.InitMesh();

c = SR.AllocR(mesh, 1, 0);
p = SR.AllocR(mesh, 1, 0);
v = SR.AllocR(mesh, 1, 1);
c_hat = SR.Alloc_Mold(c);
w_hat = SR.Alloc_Mold(c);
p_hat = SR.Alloc_Mold(p);
v_hat = SR.Alloc_Mold(v);

%SR.Fill_Random(mesh, c, -1.0, +1.0);
SR.Fill(mesh, c, 1.0, 0.0);
SR.Fill(mesh, v, 0.0, 0.0);
SR.Fill(mesh, p, 0.0, 0.0);

for time = 0:200
  if time ~= 0
    for frac = 1:10
      CahnHilliard_Step(SR, mesh, ...
        c, c, c_hat, w_hat, tau, gamma, sigma)
      %NavierStokes_Step(SR, mesh, ...
      % p, v, c, w_hat, p_hat, v_hat, tau, rho, mu);
      NavierStokes_VaRho_Step(SR, mesh, ...
        p, v, c, w_hat, p_hat, v_hat, tau, rho, 2.0*rho, mu, 2.0*mu);

      [c, c_hat] = deal(c_hat, c);
      [p, p_hat] = deal(p_hat, p);
      [v, v_hat] = deal(v_hat, v);
    end
  end

  filename = convertStringsToChars(sprintf("out/fld-%d.vtk", time));
  io = SR.IO_Begin();
  SR.IO_Add(io, v, 'velocity');
  SR.IO_Add(io, c, 'phase');
  SR.IO_Add(io, p, 'pressure');
  SR.IO_Add(io, c, 'density');
  SR.IO_Flush(io, mesh, filename);
end

function SetBCs_c(SR, mesh, c)
  SR.ApplyBCs_PureNeumann(mesh, c, 0);
  SR.ApplyBCs_Dirichlet(mesh, c, 3, 1.0);
end

function SetBCs_w(SR, mesh, w)
  SR.ApplyBCs_PureNeumann(mesh, w, 0);
end

function SetBCs_p(SR, mesh, p)
  SR.ApplyBCs_PureNeumann(mesh, p, 0);
  SR.ApplyBCs_Dirichlet(mesh, p, 2, 1.0);
  SR.ApplyBCs_Dirichlet(mesh, p, 4, 1.0);
end

function SetBCs_v(SR, mesh, v)
  SR.ApplyBCs_PureDirichlet(mesh, v, 0);
  SR.ApplyBCs_PureNeumann(mesh, v, 2);
  SR.ApplyBCs_PureNeumann(mesh, v, 3);
  SR.ApplyBCs_PureNeumann(mesh, v, 4);
end

function CahnHilliard_Step(SR, mesh, c, v, c_hat, w_hat, tau, gamma, sigma)

  % 
  % Compute a single time step of the
  % Cahn-Hilliard equation with convection term:
  %
  % ∂𝑐/∂𝑡 + ∇⋅𝑐𝒗 = Δ𝑤, 
  % 𝑤 = 𝑊'(𝑐) - 𝛾Δ𝑐, 𝑊(𝑐) = ¼(1 - 𝑐²)²,
  %
  % with the linear-implicit scheme:
  % 
  % 𝑐̃ + 𝜏∇⋅𝑐̃𝒗 = 𝑐,
  % 𝓠𝑐̂ = (1 - 2𝜏𝜎)𝑐̂ + 𝜏𝛾Δ²𝑐̂ = (1 - 2𝜏𝜎)𝑐̃ + 𝜏∇⋅𝑐̃𝒗 + 𝜏Δ𝑊'(𝑐̃),
  % 𝑤̂ = 𝑊'(𝑐̂) - 𝛾Δ𝑐̂.
  %

  SetBCs_c(SR, mesh, c);
  SetBCs_v(SR, mesh, v);

  %SR.FuncProd(mesh, w_hat, c, dWdC, NULL);
  for iCell = 1:SR.Mesh_NumCells(mesh)
    ic = SR.At(c, iCell);
    SR.SetAt(w_hat, iCell, -ic*(1.0 - ic*ic))
  end
  SetBCs_w(SR, mesh, w_hat);

  rhs = SR.Alloc_Mold(c);
  SR.Scale(mesh, rhs, c, 1.0 - 2.0*tau*sigma);
  SR.Conv(mesh, rhs, tau, c, v);
  SR.DivGrad(mesh, rhs, tau, w_hat);

  SR.Set(mesh, c_hat, c);

  while true
    [request, Qc, c] = SR.RCI_LinSolve(mesh, 'CG', '', c_hat, rhs);
    if strcmp(request, 'Done') ~= 0, break, end
    
    SetBCs_c(SR, mesh, c)

    tmp = SR.Alloc_Mold(c);

    SR.Fill(mesh, tmp, 0.0, 0.0);
    SR.DivGrad(mesh, tmp, 1.0, c);
    
    SR.Scale(mesh, Qc, c, 1.0 - 2.0*tau*sigma);
    SetBCs_w(SR, mesh, tmp);
    SR.DivGrad(mesh, Qc, tau*gamma, tmp);

    SR.Free(tmp);
  end

  SR.Free(rhs);

  SetBCs_c(SR, mesh, c_hat);
  %SR.FuncProd(mesh, w_hat, c_hat, dWdC, NULL);
  for iCell = 1:SR.Mesh_NumCells(mesh)
    ic = SR.At(c_hat, iCell);
    SR.SetAt(w_hat, iCell, -ic*(1.0 - ic*ic))
  end
  SR.DivGrad(mesh, w_hat, -gamma, c_hat);

end

function NavierStokes_Step(SR, mesh, p, v, c, w, p_hat, v_hat, tau, rho, mu)

  % 
  % Compute a single time step of the incompressible
  % Navier-Stokes equation with convection term:
  %
  % 𝜌(∂𝒗/∂𝑡 + 𝒗(∇⋅𝒗)) + ∇𝑝 = 𝜇Δ𝒗 + 𝙛,
  % ∇⋅𝒗 = 0, 𝙛 = -𝑐∇𝑤,
  %
  % with the semi-implicit scheme:
  % 
  % 𝒗̂ ← 𝒗 - 𝜏𝒗(∇⋅𝒗) + (𝜏𝜇/𝜌)Δ𝒗 + (𝜏/𝜌)𝙛,
  % Δ𝑝̂ = 𝜌/𝜏∇⋅𝒗,
  % 𝒗̂ ← 𝒗̂ - (𝜏/𝜌)∇𝑝̂.
  % 

  %
  % Compute 𝒗̂ prediction.
  %
  SetBCs_v(SR, mesh, v);
  SR.Set(mesh, v_hat, v);
  SR.Conv(mesh, v_hat, tau, v, v);
  SR.DivGrad(mesh, v_hat, tau*mu/rho, v);

  %
  % Compute (𝜏/𝜌)𝙛 = -(𝜏/𝜌)𝑐∇𝑤, 𝒗̂ ← 𝒗̂ + (𝜏/𝜌)𝙛.
  %
  f = SR.Alloc_Mold(v);

  SetBCs_w(SR, mesh, w);
  SR.Fill(mesh, f, 0.0, 0.0);
  SR.Grad(mesh, f, 1.0, w);
  SR.Mul(mesh, f, c, f);

  SR.Add(mesh, v_hat, v_hat, f, tau/rho, 1.0);

  SR.Free(f);

  %
  % Solve pressure equation and correct 𝒗̂.
  % 
  SetBCs_v(SR, mesh, v_hat);
  rhs = SR.Alloc_Mold(p);
  SR.Fill(mesh, rhs, 0.0, 0.0);
  SR.Div(mesh, rhs, -rho/tau, v_hat);

  SR.Set(mesh, p_hat, p);
  while true
    [request, Lp, p] = SR.RCI_LinSolve(mesh, 'CG', '', p_hat, rhs);
    if strcmp(request, 'Done') ~= 0, break, end
    
    SetBCs_p(SR, mesh, p);
    SR.Fill(mesh, Lp, 0.0, 0.0);
    SR.DivGrad(mesh, Lp, 1.0, p);
  end

  SR.Free(rhs);

  SetBCs_p(SR, mesh, p_hat);
  SR.Grad(mesh, v_hat, tau/rho, p_hat);

end

function NavierStokes_VaRho_Step(SR, mesh, p, v, c, w, ...
  p_hat, v_hat, tau, rho_1, rho_2, mu_1, mu_2)

  % 
  % Compute a single time step of the incompressible
  % Navier-Stokes equation with convection term:
  %
  % 𝜌(∂𝒗/∂𝑡 + 𝒗(∇⋅𝒗)) + ∇𝑝 = 𝙛,
  % ∇⋅𝒗 = 0, 𝙛 = 𝜇Δ𝒗 - 𝑐∇𝑤,
  % 𝜌 = ½𝜌₁(1 - 𝑐) + ½𝜌₂(1 + 𝑐),
  % 𝜇 = ½𝜇₁(1 - 𝑐) + ½𝜇₂(1 + 𝑐).
  %
  % with the semi-implicit scheme:
  % 
  % 𝒗̂ ← 𝒗 - 𝜏𝒗(∇⋅𝒗) + 𝜏𝙛/𝜌,
  % ∇⋅(∇𝑝̂/𝜌) = ∇⋅𝒗/𝜏,
  % 𝒗̂ ← 𝒗̂ - (𝜏/𝜌)∇𝑝̂.
  % 

  %
  % Compute 𝜌, 𝜇:
  % 𝜌 = ½(𝜌₁ + 𝜌₂) + ½(𝜌₂ - 𝜌₁)𝑐.
  % 𝜇 = ½(𝜇₁ + 𝜇₂) + ½(𝜇₂ - 𝜇₁)𝑐.
  %
  rho = SR.Alloc_Mold(c);
  SR.Fill(mesh, rho, 0.5*(rho_1 + rho_2), 0.0);
  SR.Add(mesh, rho, rho, c, 0.5*(rho_2 - rho_1), 1.0);

  rho_inv = SR.Alloc_Mold(rho);
  %SR.Inv(mesh, rho_inv, rho);
  for iCell = 1:SR.Mesh_NumCells(mesh)
    ir = SR.At(rho, iCell);
    SR.SetAt(rho_inv, iCell, 1.0/ir)
  end

  mu = SR.Alloc_Mold(c);
  SR.Fill(mesh, mu, 0.5*(mu_1 + mu_2), 0.0);
  SR.Add(mesh, mu, mu, c, 0.5*(mu_2 - mu_1), 1.0);

  %
  % Compute 𝒗̂ prediction.
  %
  SetBCs_v(SR, mesh, v);
  SR.Set(mesh, v_hat, v);
  SR.Conv(mesh, v_hat, tau, v, v);

  %
  % Compute 𝙛 = 𝜇Δ𝒗 - 𝑐∇𝑤, 𝒗̂ ← 𝒗̂ + (𝜏/𝜌)𝙛.
  %
  f = SR.Alloc_Mold(v);

  SetBCs_w(SR, mesh, w);
  SR.Fill(mesh, f, 0.0, 0.0);
  SR.Grad(mesh, f, 1.0, w);
  SR.Mul(mesh, f, c, f);

  g = SR.Alloc_Mold(v);
  SR.Fill(mesh, g, 0.0, 0.0);
  SR.DivGrad(mesh, g, 1.0, v);
  SR.Mul(mesh, g, mu, g);
  SR.Add(mesh, f, f, g, 1.0, 1.0);
  SR.Free(g);

  SR.Mul(mesh, f, rho_inv, f);

  SR.Add(mesh, v_hat, v_hat, f, tau, 1.0);

  SR.Free(f);
  
  %
  % Solve pressure equation and correct 𝒗̂.
  % 
  SetBCs_v(SR, mesh, v_hat);
  rhs = SR.Alloc_Mold(p);
  SR.Fill(mesh, rhs, 0.0, 0.0);
  SR.Div(mesh, rhs, -1.0/tau, v_hat);

  SR.Set(mesh, p_hat, p);
  while true
    [request, Lp, p] = SR.RCI_LinSolve(mesh, 'CG', '', p_hat, rhs);
    if strcmp(request, 'Done') ~= 0, break, end
    
    SetBCs_p(SR, mesh, p);
    SR.Fill(mesh, Lp, 0.0, 0.0);
    SR.DivKGrad(mesh, Lp, 1.0, rho_inv, p);
  end

  SR.Free(rhs);

  SetBCs_p(SR, mesh, p_hat);
  g = SR.Alloc_Mold(v);
  SR.Fill(mesh, g, 0.0, 0.0);
  SR.Grad(mesh, g, 1.0, p_hat);
  SR.Mul(mesh, g, rho_inv, g);
  SR.Add(mesh, v_hat, v_hat, g, tau, 1.0);
  SR.Free(g);

  SR.Free(rho);
  SR.Free(rho_inv);
  SR.Free(mu);

end
