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
use StormRuler_FiniteDifferences
use StormRuler_KrylovSolvers
  
implicit none
   
type CahnHilliardParams
  real(8) :: EpsSqr
contains
  procedure, nopass :: F => CahnHilliardParams_DoubleWell
  procedure, nopass :: dFdC => CahnHilliardParams_ddCDoubleWell
end type CahnHilliardParams

contains  
pure function CahnHilliardParams_DoubleWell(C) result(F)
  real(8), intent(in) :: C
  real(8) :: F
  F = 0.25*(1 - C**2)**2
end function CahnHilliardParams_DoubleWell
pure function CahnHilliardParams_ddCDoubleWell(C) result(F)
  real(8), intent(in) :: C
  real(8) :: F
  F = -C*(1 - C**2)
end function CahnHilliardParams_ddCDoubleWell

!! -----------------------------------------------------------------  
!! Estimate the Ginzburg-Landau free energy functional value.
function CahnHilliard_FreeEnergy(nX,nY,dX,dY, C,CHPhysParams) result(E)
  integer, intent(in) :: nX, nY
  real(8), intent(in) :: dX, dY
  real(8), dimension(0:nX+1,0:nY+1), intent(inout) :: C
  class(CahnHilliardParams), intent(in) :: CHPhysParams
  real(8) :: E
  block
    integer :: iX, iY
    E = 0
    ! ----------------------
    ! E_phobic := <1,F(C)ij>
    !$omp parallel do private(iX,iY) collapse(2) reduction(+:E)
    do iY = 1, nY
      do iX = 1, nX
        E = E + dX*dY*CHPhysParams%F(C(iX,iY))
      end do
    end do
    !$omp end parallel do
    ! ----------------------
    ! E_philic := EpsSqr/2*<(Nabla)F(C),(Nabla)F(C)>
    !$omp parallel do private(iX,iY) collapse(2) reduction(+:E)
    do iY = 1, nY
      do iX = 0, nX
        E = E + 0.5*CHPhysParams%EpsSqr * (dY/dX)*(C(iX+1,iY) - C(iX,iY))
      end do
    end do
    !$omp end parallel do
    !$omp parallel do private(iX,iY) collapse(2) reduction(+:E)
    do iY = 0, nY
      do iX = 1, nX
        E = E + 0.5*CHPhysParams%EpsSqr * (dX/dY)*(C(iX,iY+1) - C(iX,iY))
      end do
    end do
    !$omp end parallel do
    ! ----------------------
  end block
end function CahnHilliard_FreeEnergy

!! -----------------------------------------------------------------  
!! Estimate the RHS for the implicit scheme SLAE 
!! b ← φ + dt⋅(ΔF'(φ)-∇⋅φv).   
subroutine CahnHilliard_ImplicitSchemeRHS(mesh, C,v,B,physParams)
  class(Mesh2D), intent(in) :: mesh
  real(8), dimension(:), intent(inout) :: C
  real(dp), intent(inout) :: v(:,:)
  real(8), dimension(:), intent(out) :: B
  class(CahnHilliardParams), intent(in) :: physParams
  associate(dt=>mesh%dt)
    ! ----------------------
    ! b ← φ,
    ! b ← b - dt⋅∇⋅φv,
    ! b ← b + dt⋅ΔF'(φ).
    call Set(mesh, B,C)
    call FDM_Convection_Central(mesh, B,dt,c,v)
    call FDM_LaplacianF(mesh, B,dt,CahnHilliardParams_ddCDoubleWell,C)
  end associate
end subroutine CahnHilliard_ImplicitSchemeRHS

!! -----------------------------------------------------------------  
!! Apply a operator estimate for the implicit scheme SLAE.
subroutine CahnHilliard_ImplicitSchemeOperator(mesh,U,C,CHPhysParams)
  class(Mesh2D), intent(in) :: mesh
  real(8), dimension(:), intent(inout) :: u,c
  class(CahnHilliardParams), intent(in) :: CHPhysParams
  associate(dt=>mesh%dt,eps=>CHPhysParams%EpsSqr)
    ! ----------------------
    ! u ← c
    ! u ← u + ε⋅dt⋅Δ²c.
    call Set(mesh,u,c)
    call FDM_Bilaplacian(mesh,u,eps*dt,c)
  end associate
end subroutine CahnHilliard_ImplicitSchemeOperator
subroutine CahnHilliard_ImplicitSchemeOperatorHelper(mesh,u,c,aCHPhysParams)
  class(Mesh2D), intent(in) :: mesh
  real(8), dimension(:), intent(inout) :: u,c
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
  class(Mesh2D), intent(in) :: mesh
  real(8), dimension(:), intent(inout) :: C
  real(dp), intent(inout) :: v(:,:)
  class(CahnHilliardParams), intent(in) :: CHPhysParams
  class(ConvParameters), allocatable :: Params
  real(8), allocatable :: b(:)
  ! ----------------------
  ! Initialize iterations.
  allocate(Params)
  call Params%Init(mesh%Dx(1)*mesh%Dx(1)*1.0D-4, mesh%Dx(1)*mesh%Dx(1)*1.0D-4, 100000)
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
  class(Mesh2D), intent(in) :: mesh
  real(8), dimension(:), intent(inout) :: C,S
  real(dp), intent(inout) :: v(:,:)
  class(CahnHilliardParams), intent(in) :: CHPhysParams
  call CahnHilliard_ImplicitSchemeSolve(mesh, C,v,CHPhysParams)
  ! s ← F'(c)
  ! s ← s + (-ε)⋅Δc
  call ApplyFunc(mesh,s,c,CahnHilliardParams_ddCDoubleWell)
  call FDM_Laplacian(mesh,s,-CHPhysParams%EpsSqr,c)
end subroutine CahnHilliard_ImplicitSchemeStep
  
end module StormRuler_CahnHilliard
  