clear all;

global tau; global mu; global rho;

rho = 1.0;
mu = 0.1;

tau = (pi/50)*(pi/50); 
Gamma = 0.01;
sigma = 1.0;

global SR
SR = StormRuler_Matlab();

mesh = SR.InitMesh();

c = SR.AllocR(mesh, 1, 0);
p = SR.AllocR(mesh, 1, 0);
v = SR.AllocR(mesh, 1, 1);
c_hat = SR.Alloc_Mold(c);
w_hat = SR.Alloc_Mold(c);
p_hat = SR.Alloc_Mold(p);
v_hat = SR.Alloc_Mold(v);

%SR.Fill_Random(mesh, c, -1.0, +1.0);
SR.Fill(mesh, c, 0.0, 0.0);
SR.Fill(mesh, w_hat, 0.0, 0.0);
SR.Fill(mesh, v, 0.0, 0.0);
SR.Fill(mesh, p, 0.0, 0.0);

for time = 0:200
  if time ~= 0
    for frac = 1:10
      time
      frac
      NavierStokes_Step(mesh, p, v, c, w_hat, p_hat, v_hat);

      %swap(c, c_hat);
      swap(p, p_hat);
      swap(v, v_hat);
    end
  end
end

function NavierStokes_Step(mesh, p, v, c, w, p_hat, v_hat)

  % 
  % Compute a single time step of the incompressible
  % Navier-Stokes equation with convection term:
  %
  % 𝜌(∂𝒗/∂𝑡 + 𝒗(∇⋅𝒗)) + ∇𝑝 = 𝜇Δ𝒗 + 𝙛,
  % ∇⋅𝒗 = 0, 𝙛 = 𝑐∇𝑤,
  %
  % with the semi-implicit scheme:
  % 
  % 𝒗̂ ← 𝒗 - 𝜏𝒗(∇⋅𝒗) + (𝜏𝜇/𝜌)Δ𝒗 + (𝜏/𝜌)𝙛,
  % Δ𝑝̂ = 𝜌/𝜏∇⋅𝒗,
  % 𝒗̂ ← 𝒗̂ - (𝜏/𝜌)∇𝑝̂.
  % 

  global SR;
  global tau; global mu; global rho;

  %
  % Compute 𝒗̂ prediction.
  %
  SR.ApplyBCs_PureNeumann(mesh, w, 0);
  SR.ApplyBCs_PureDirichlet(mesh, v, 0);
  SR.ApplyBCs_PureNeumann(mesh, v, 2);
  SR.ApplyBCs_PureNeumann(mesh, v, 3);
  SR.ApplyBCs_PureNeumann(mesh, v, 4);

  SR.Set(mesh, v_hat, v);
  SR.Conv(mesh, v_hat, tau, v, v);
  SR.DivGrad(mesh, v_hat, tau*mu/rho, v);

  %
  % Compute 𝙛 = -𝑐∇𝑤, 𝒗̂ ← 𝒗̂ + 𝙛.
  %
  f = SR.Alloc_Mold(v);

  SR.Fill(mesh, f, 0.0, 0.0);
  SR.Grad(mesh, f, 1.0, w);
  SR.Mul(mesh, f, c, f);

  SR.Add(mesh, v_hat, v_hat, f, tau/rho, 1.0);

  SR.Free(f);

  SR.ApplyBCs_PureDirichlet(mesh, v_hat, 0);
  SR.ApplyBCs_PureNeumann(mesh, v_hat, 2);
  SR.ApplyBCs_PureNeumann(mesh, v_hat, 3);
  SR.ApplyBCs_PureNeumann(mesh, v_hat, 4);

  %
  % Solve pressure equation and correct 𝒗̂.
  % 
  RHS = SR.Alloc_Mold(p);
  SR.Fill(mesh, RHS, 0.0, 0.0);
  SR.Div(mesh, RHS, -rho/tau, v_hat);

  SR.Set(mesh, p_hat, p);
  while true
    [request, Lp, p] = SR.RCI_LinSolve(mesh, 'CG', '', p_hat, RHS);
    if strcmp(request, 'Done') ~= 0
      break
    end
    
    SR.ApplyBCs_PureNeumann(mesh, p, 0);
    SR.ApplyBCs_Dirichlet(mesh, p, 2, 1.0);
    SR.ApplyBCs_Dirichlet(mesh, p, 4, 10.0);

    SR.Fill(mesh, Lp, 0.0, 0.0);
    SR.DivGrad(mesh, Lp, 1.0, p);
  end

  SR.ApplyBCs_PureNeumann(mesh, p_hat, 0);
  SR.ApplyBCs_Dirichlet(mesh, p_hat, 2, 1.0);
  SR.ApplyBCs_Dirichlet(mesh, p_hat, 4, 10.0);

  SR.Grad(mesh, v_hat, tau/rho, p_hat);

  SR.Free(RHS);
end

function [b, a] = swap(a, b)
end
