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
module StormRuler_CahnHilliard

use StormRuler_Helpers
use StormRuler_Mesh
use StormRuler_BLAS
use StormRuler_FDM_Operators
use StormRuler_ConvParams
use StormRuler_Solvers_Krylov
  
implicit none
   
type CahnHilliardParams
  real(dp) :: EpsSqr
contains
  procedure, nopass :: F => CahnHilliardParams_DoubleWell
  procedure, nopass :: dFdC => CahnHilliardParams_ddCDoubleWell
end type CahnHilliardParams

contains  

pure function CahnHilliardParams_DoubleWell(C) result(F)
  real(dp), intent(in) :: C
  real(dp) :: F
  F = 0.25*(1 - C**2)**2
end function CahnHilliardParams_DoubleWell

pure function CahnHilliardParams_ddCDoubleWell(C) result(F)
  real(dp), intent(in) :: C
  real(dp) :: F
  F = -C*(1 - C**2)
end function CahnHilliardParams_ddCDoubleWell

!! -----------------------------------------------------------------  
!! Estimate the RHS for the implicit scheme SLAE 
!! b ← φ + dt⋅(ΔF'(φ)-∇⋅φv).   
subroutine CahnHilliard_ImplicitSchemeRHS(mesh, C,v,B,physParams)
  class(tMesh), intent(in) :: mesh
  real(dp), dimension(:), intent(inout), target :: C
  real(dp), intent(inout), target :: v(:,:)
  real(dp), dimension(:), intent(out), target :: B
  class(CahnHilliardParams), intent(in) :: physParams
  associate(dt=>mesh%dt)
    ! ----------------------
    ! b ← φ,
    ! b ← b - dt⋅∇⋅φv,
    ! b ← b + dt⋅ΔF'(φ).
    call Set(mesh, B,C)
    call FDM_Convection_Central(mesh, B,dt,c,v)
    call FDM_LaplacianF_Central(mesh, B,dt,CahnHilliardParams_ddCDoubleWell,C)
  end associate
end subroutine CahnHilliard_ImplicitSchemeRHS

!! -----------------------------------------------------------------  
!! Apply a operator estimate for the implicit scheme SLAE.
subroutine CahnHilliard_ImplicitSchemeOperator(mesh,U,C,CHPhysParams)
  class(tMesh), intent(in) :: mesh
  real(dp), dimension(:), intent(in), target :: c
  real(dp), dimension(:), intent(inout), target :: u
  class(CahnHilliardParams), intent(in) :: CHPhysParams
  associate(dt=>mesh%dt,eps=>CHPhysParams%EpsSqr)
    ! ----------------------
    ! u ← c
    ! u ← u + ε⋅dt⋅Δ²c.
    call Set(mesh,u,c)
    call FDM_Bilaplacian_Central(mesh,u,eps*dt,c)
  end associate
end subroutine CahnHilliard_ImplicitSchemeOperator
subroutine CahnHilliard_ImplicitSchemeOperatorHelper(mesh,u,c,aCHPhysParams)
  class(tMesh), intent(in) :: mesh
  real(dp), dimension(:), intent(in), target :: c
  real(dp), dimension(:), intent(inout), target :: u
  class(*), intent(in) :: aCHPhysParams
  select type(aCHPhysParams)
    class is (CahnHilliardParams)
      call CahnHilliard_ImplicitSchemeOperator(mesh,u,c,aCHPhysParams)
    class default
      error stop 1
  end select
end subroutine CahnHilliard_ImplicitSchemeOperatorHelper

!! -----------------------------------------------------------------  
!! Solve the SLAE of the implicit scheme (using Conjugate Gradients method).
subroutine CahnHilliard_ImplicitSchemeSolve(mesh, C,v, CHPhysParams)
  class(tMesh), intent(in) :: mesh
  real(dp), dimension(:), intent(inout), target :: C
  real(dp), intent(inout) :: v(:,:)
  class(CahnHilliardParams), intent(in) :: CHPhysParams
  class(tConvParams), allocatable :: Params
  real(dp), allocatable, target :: b(:)
  ! ----------------------
  ! Initialize iterations.
  allocate(Params)
  call Params%Init(mesh%dl(1)*mesh%dl(1)*1.0D-4, mesh%dl(1)*mesh%dl(1)*1.0D-4, 100000)
  ! ----------------------
  allocate(b, mold=c)
  call CahnHilliard_ImplicitSchemeRHS(mesh,c,v,b,CHPhysParams)
  call Solve_BiCGStab(mesh,c,b&
    ,CahnHilliard_ImplicitSchemeOperatorHelper,CHPhysParams &
    ,Params)
end subroutine CahnHilliard_ImplicitSchemeSolve

!! -----------------------------------------------------------------  
!! Compute Cahn-Hilliard time step with an implicit scheme.
subroutine CahnHilliard_ImplicitSchemeStep(mesh, C,S,v, CHPhysParams)
  class(tMesh), intent(in) :: mesh
  real(dp), dimension(:), intent(inout), target :: C,S
  real(dp), intent(inout) :: v(:,:)
  class(CahnHilliardParams), intent(in) :: CHPhysParams
  call CahnHilliard_ImplicitSchemeSolve(mesh, C,v,CHPhysParams)
  ! s ← F'(c)
  ! s ← s + (-ε)⋅Δc
  call FuncProd(mesh,s,c,CahnHilliardParams_ddCDoubleWell)
  call FDM_Laplacian_Central(mesh,s,-CHPhysParams%EpsSqr,c)
end subroutine CahnHilliard_ImplicitSchemeStep
  
end module StormRuler_CahnHilliard
  