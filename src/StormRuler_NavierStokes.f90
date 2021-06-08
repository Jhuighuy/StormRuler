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
use StormRuler_FiniteDifferences
use StormRuler_KrylovSolvers

implicit none
    
contains
   
subroutine SolvePoisson(mesh,u,f)
  class(Mesh2D), intent(in) :: mesh
  real(8), dimension(:), intent(inout) :: u(:),f(:)
  class(ConvParameters), allocatable :: Params
  ! ----------------------
  ! Initialize iterations.
  allocate(Params)
  call Params%Init(mesh%Dx(1)*mesh%Dx(1)*1.0D-4, mesh%Dx(1)*mesh%Dx(1)*1.0D-4, 100000)
  ! ----------------------
  call Fill(mesh,u)
  call Solve_BiCGStab(mesh,u,f&
    ,PoissonOperator,Params,Params)
contains
  subroutine PoissonOperator(mesh,f,u,opParams)
    class(Mesh2D), intent(in) :: mesh
    real(8), dimension(:), intent(inout) :: f,u
    class(*), intent(in) :: opParams
    call Fill(mesh,f)
    call FDM_Laplacian_Central(mesh,f,1.0_dp,u)
  end subroutine PoissonOperator
end subroutine SolvePoisson

subroutine NavierStokes_PredictVelocity(mesh,v,p,c,s,g,nu,rho,beta)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(inout) :: v(:,:),g(:,:),p(:),c(:),s(:)
  real(dp), intent(in) :: nu,rho,beta
  ! >>>>>>>>>>>>>>>>>>>>>>
  real(dp), allocatable :: z(:,:),w(:,:),d(:),q(:),h(:)
  allocate(w,z, mold=v)
  allocate(d,h,q, mold=p)
  associate(dt=>mesh%dt)
    ! z ← 0,
    ! z ← z - dt/ρ⋅∇s,
    ! z ← c⋅z
    call Fill(mesh,z)
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
    call Fill(mesh,h)
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
    