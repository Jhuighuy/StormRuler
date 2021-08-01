!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!! Copyright (C) 2021 Oleg Butakov
!! 
!! Permission is hereby granted, free of charge, to any person 
!! obtaining a copy of this software and associated documentation 
!! files (the "Software"), to deal in the Software without 
!! restriction, including without limitation the rights  to use, 
!! copy, modify, merge, publish, distribute, sublicense, and/or
!! sell copies of the Software, and to permit persons to whom the  
!! Software is furnished to do so, subject to the following 
!! conditions:
!! 
!! The above copyright notice and this permission notice shall be 
!! included in all copies or substantial portions of the Software.
!! 
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
!! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
!! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
!! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, 
!! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
!! OTHER DEALINGS IN THE SOFTWARE.
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
module StormRuler_NavierStokes

use StormRuler_Helpers
use StormRuler_BLAS
use StormRuler_FDM_Operators
use StormRuler_FDM_BCs
use StormRuler_ConvParams
use StormRuler_KrylovSolvers

implicit none
    
contains
   
subroutine SolvePoisson(mesh,u,f)
  class(tMesh), intent(in) :: mesh
  real(8), dimension(:), intent(inout) :: u(:),f(:)
  class(tConvParams), allocatable :: Params
  ! ----------------------
  ! Initialize iterations.
  allocate(Params)
  call Params%Init(mesh%Dl(1)*mesh%Dl(2)*1.0D-4, mesh%Dl(1)*mesh%Dl(1)*1.0D-4, 100000)
  !call Params%Init(1.0D-8, 1.0D-8, 100000)
  ! ----------------------
  !call Fill(mesh,u,0.0_dp)
  call Solve_CG(mesh,u,f,PoissonOperator,Params,Params)
  !call Solve_BiCGStab(mesh,u,f,PoissonOperator,Params,Params)
contains
  subroutine PoissonOperator(mesh,v,u,opParams)
    class(tMesh), intent(in) :: mesh
    real(8), dimension(:), intent(inout) :: v,u
    class(*), intent(in) :: opParams
    !call FDM_ApplyBCs(mesh,2,u,0.00_dp,1.0_dp,0.0_dp)
    !call FDM_ApplyBCs(mesh,1,u,1.00_dp,0.0_dp,0.0_dp)
    call Fill(mesh,v,0.0_dp)
    call FDM_Laplacian_Central(mesh,v,1.0_dp,u)
  end subroutine PoissonOperator
end subroutine SolvePoisson

subroutine NavierStokes_PredictVelocity(mesh,v,p,c,s,g,nu,rho,beta)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(inout) :: v(:,:),g(:,:),p(:),c(:),s(:)
  real(dp), intent(in) :: nu,rho,beta
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), allocatable :: z(:,:),w(:,:),d(:),h(:)
  allocate(w,z,mold=v)
  allocate(d,h,mold=p)
  associate(dt=>mesh%dt)
    ! z ← 0,
    ! z ← z - dt/ρ⋅∇s,
    ! z ← c⋅z
    call Fill(mesh,z,0.0_dp)
    call FDM_Gradient_Central(mesh,z,dt/rho,s)
    call Mul(mesh,z,c,z)
    ! w ← v + z,
    ! w ← w - dt⋅(v⋅∇)v,
    ! w ← w - β⋅dt/ρ⋅∇p,
    ! w ← w + dt⋅ν⋅Δv.
    call Add(mesh,w,v,z)
    call FDM_Convection_Central(mesh,w,dt,v,v)
    call FDM_Gradient_Forward(mesh,w,beta*dt/rho,p)
    call FDM_Laplacian_Central(mesh,w,dt*nu,v)
    !call Add(mesh,w,w,g,dt)
    ! h ← 0,
    ! h ← h + (-ρ/dt)⋅(∇⋅w),
    ! solve Δd = h.
    call Fill(mesh,h,0.0_dp)
    call FDM_Divergence_Backward(mesh,h,-rho/dt,w)
    call SolvePoisson(mesh,d,h)
    ! p ← d + β⋅p,
    ! w ← w - dt/ρ⋅∇d,
    ! v ← w.
    call Add(mesh,p,d,p,beta)
    call FDM_Gradient_Forward(mesh,v,dt/rho,d)
    call Set(mesh,v,w)
  end associate
end subroutine NavierStokes_PredictVelocity

end module StormRuler_NavierStokes
