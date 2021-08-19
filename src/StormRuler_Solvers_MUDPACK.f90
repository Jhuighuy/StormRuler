!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
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
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!
module StormRuler_Solvers_MUDPACK

#$use 'StormRuler_Params.fi'
#$if HAS_MUDPACK

use StormRuler_Parameters, only: dp, ip
use StormRuler_Mesh, only: tMesh
use StormRuler_ConvParams, only: tConvParams

use StormRuler_Helpers

use, intrinsic :: iso_fortran_env, only: error_unit

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

abstract interface
  function tFunc_mud2sa(x, y)
    import dp
    real(dp), intent(in) :: x, y
    real(dp) :: tFunc_mud2sa
  end function tFunc_mud2sa
  subroutine tBoundary_mud2sa(kbdy, xory, alfa, gbdy)
    import dp
    integer kbdy
    real(dp) xory, alfa, gbdy
  end subroutine tBoundary_mud2sa
end interface

interface
  subroutine mud2sa(iparam, dparam, work,  &
      & sigma_x, sigma_y, lambda_wrapper, bndyc, rhs, phi, mgopt, ierror)
    import tFunc_mud2sa, tBoundary_mud2sa
    ! <<<<<<<<<<<<<<<<<<<<<<
    integer, intent(inout) :: iparam(17), mgopt(4), ierror
    double precision, intent(inout) :: dparam(6)
    double precision :: work(*), rhs(*), phi(*)
    procedure(tFunc_mud2sa) :: sigma_x, sigma_y, lambda_wrapper
    procedure(tBoundary_mud2sa) :: bndyc
    ! >>>>>>>>>>>>>>>>>>>>>>
  end subroutine mud2sa
end interface

private :: mud2sa!, mud3sa

interface Solve_DivWGrad_MUDPACK
#$do rank = 0, 0*NUM_RANKS
  module procedure Solve_DivWGrad_MUDPACK$rank
#$end do
end interface Solve_DivWGrad_MUDPACK

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Split size into the coefficient and exponent: 
!! coef*( 2**(exp-1) ) + 1 == size, coef > 1.
!! ----------------------------------------------------------------- !!
subroutine SplitSize_MUDPACK(size, coef, exp)
  ! <<<<<<<<<<<<<<<<<<<<<<
  integer, intent(in) :: size
  integer, intent(out) :: coef, exp
  ! >>>>>>>>>>>>>>>>>>>>>>

  ! TODO: optimize me.
  do exp = 16, 2, -1
    associate(p => 2**(exp-1))
      if (mod(size - 1, p) == 0) then
        coef = (size - 1)/p
        if (coef /= 1) then
          return
        end if
      end if
    end associate
  end do
  coef = size - 1; exp = 1
end subroutine SplitSize_MUDPACK

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Solve second-order approximate variable weight Laplacian
!! equation: ∇⋅(w∇u) - λu = b on the rectangular domain
!! using the MUDPACK library.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, 0*NUM_RANKS
subroutine Solve_DivWGrad_MUDPACK$rank(mesh, u, lambda, w, b, params)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: lambda, w(@:,:), b(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  type(tConvParams), intent(inout) :: params
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer(ip) :: iCell
  integer(ip) :: iCell_x, iCell_y
  integer(ip) :: iNode_x, iNode_y

  real(dp), allocatable :: w_MD_cc(:,:), b_MD_cc(:,:)
  real(dp), allocatable :: u_MD_nc(:,:), b_MD_nc(:,:)

  integer :: ierror, iparam(17), mgopt(4)
  real(dp) :: dparam(6)
  real(dp), allocatable :: work(:)

  ! ----------------------
  ! Initialize and configure MUDPACK.
  ! ----------------------
  iparam(:) = 0; dparam(:) = 0.0_dp
  associate( &
    &                 phase => iparam( 1), & ! intl
    &    domainLength_x_min => dparam( 1), & ! xa
    &    domainLength_x_max => dparam( 2), & ! xb
    &    domainLength_y_min => dparam( 3), & ! yc
    &    domainLength_y_max => dparam( 4), & ! yd
    &             bcs_x_min => iparam( 2), & ! nxa
    &             bcs_x_max => iparam( 3), & ! nxb
    &             bcs_y_min => iparam( 4), & ! nyc
    &             bcs_y_max => iparam( 5), & ! nyd
    &            numNodes_x => iparam(10), & ! nx
    &            numNodes_y => iparam(11), & ! ny
    &       numNodes_x_coef => iparam( 6), & ! ixp > 1
    &       numNodes_y_coef => iparam( 7), & ! jyq > 1
    &        numNodes_x_exp => iparam( 8), & ! iex
    &        numNodes_y_exp => iparam( 9), & ! jey
    & initialGuessSpecified => iparam(12), & ! iguess
    &          maxNumCycles => iparam(13), & ! maxcy
    &      relaxationMethod => iparam(14), & ! method
    &    multigridCycleType =>  mgopt( 1), & ! kcycle
    &         workspaceSize => iparam(15), & ! length
    & computedWorkspaceSize => iparam(16), & !
    &     relativeTolerance => dparam( 5) )  ! tolmax
    ! ----------------------
    ! Specify domain and BCs.
    ! ----------------------
    domainLength_x_min = 0.0_dp
    domainLength_y_min = 0.0_dp
    domainLength_x_max = mesh%dl(1)*mesh%MDIndexBounds(1)
    domainLength_y_max = mesh%dl(3)*mesh%MDIndexBounds(2)
    ! value = 0: periodic boundary conditions,
    ! value = 1: prescribed boundary conditions, specified in the solution field,
    ! value = 2: mixed derivative boundary conditions.
    ! TODO: set the real boundary conditions.
    bcs_x_min = 0; bcs_x_max = 0
    bcs_y_min = 0; bcs_y_max = 0

    ! ----------------------
    ! Set up mesh size.
    ! ----------------------
    ! The following equalities should hold:
    ! numNodes_x = numNodes_x_coef*( 2**(numNodes_x_exp-1) ) + 1
    ! numNodes_y = numNodes_x_coef*( 2**(numNodes_y_exp-1) ) + 1
    numNodes_x = mesh%MDIndexBounds(1) + 1
    numNodes_y = mesh%MDIndexBounds(2) + 1
    call SplitSize_MUDPACK(numNodes_x, numNodes_x_coef, numNodes_x_exp)
    call SplitSize_MUDPACK(numNodes_y, numNodes_y_coef, numNodes_y_exp)

    ! ----------------------
    ! Set up initial guess and method.
    ! ----------------------
    initialGuessSpecified = 0
    ! value = 0: point relaxation
    ! value = 1: line relaxation in X direction,
    ! value = 2: line relaxation in Y direction,
    ! value = 3: line relaxation in both X and Y directions.
    relaxationMethod = 0
    maxNumCycles = params%MaxNumIterations
    relativeTolerance = params%RelativeTolerance
    ! value = 0: default cycle (W).
    ! value = 1: V cycle.
    ! value = 2: W cycle.
    multigridCycleType = 0

    ! ----------------------
    ! Allocate the workspace.
    ! ----------------------
    ! Set estimated size to 0 to force the MUDPACK
    ! to fail with error 9 and read the real workspace size.
    workspaceSize = 0
    phase = 0; call mud2sa(iparam, dparam, work, &
      & w_wrapper, w_wrapper, lambda_wrapper, boundary, b_MD_nc, u_MD_nc, mgopt, ierror)
    if (ierror /= 9) then
      write(error_unit, *) &
        & 'MUD2SA FAILED WITH ERROR CODE OTHER THAN `9`', ierror
    end if
    workspaceSize = computedWorkspaceSize
    allocate( work(workspaceSize) )
  end associate

  ! ----------------------
  ! Unwrap and interpolate the fields from cell to nodes.
  ! ----------------------
  associate( &
    & numCells_x => mesh%MDIndexBounds(1), &
    & numCells_y => mesh%MDIndexBounds(2) )
    allocate(w_MD_cc(numCells_x, numCells_y))
    allocate(b_MD_cc(numCells_x, numCells_y))
  end associate
  !$omp parallel do schedule(static) &
  ! 'w' is unwraped for fast lookup,
  ! 'b' is unrapped to be interpolated to nodes.
  !$omp & default(private) shared(mesh, w, w_MD_cc, b, b_MD_cc)
  do iCell = 1_ip, mesh%NumCells
    associate(iCellMD => mesh%CellMDIndex(:,iCell))
      w_MD_cc(iCellMD(1), iCellMD(2)) = w(iCell)
      b_MD_cc(iCellMD(1), iCellMD(2)) = b(iCell)
    end associate
  end do
  !$omp end parallel do
  ! Interpolate 'b' to nodes.
  associate( &
    & numNodes_x => mesh%MDIndexBounds(1) + 1, &
    & numNodes_y => mesh%MDIndexBounds(2) + 1 )
    allocate( &
      & u_MD_nc(numNodes_x, numNodes_y), &
      & b_MD_nc(numNodes_x, numNodes_y) )
    call Interpolate2D_CellToNode0(b_MD_nc, b_MD_cc)
    deallocate(b_MD_cc)
  end associate

  ! ----------------------
  ! Perform the discretization phase.
  ! ----------------------
  associate(phase => iparam(1))
    phase = 0; call mud2sa(iparam, dparam, work, &
      & w_wrapper, w_wrapper, lambda_wrapper, boundary, b_MD_nc, u_MD_nc, mgopt, ierror)
    if (ierror < 0) then
      write(error_unit, *) &
        & 'MUD2SA DISCRETIZATION PHASE WARNING, IERROR=', ierror
    else if (ierror > 0) then
      write(error_unit, *) &
        & 'MUD2SA DISCRETIZATION PHASE FAILED, IERROR=', ierror
      error stop 1
    end if
  end associate

  ! ----------------------
  ! Perform the approximation phase.
  ! ----------------------
  associate(phase => iparam(1))
    phase = 1; call mud2sa(iparam, dparam, work, &
      & w_wrapper, w_wrapper, lambda_wrapper, boundary, b_MD_nc, u_MD_nc, mgopt, ierror)
    if (ierror < 0) then
      write(error_unit, *) &
        & 'MUD2SA APPROXIMATION PHASE WARNING, IERROR=', ierror
    else if (ierror > 0) then
      write(error_unit, *) &
        & 'MUD2SA APPROXIMATION PHASE FAILED, IERROR=', ierror
      error stop 1
    end if
  end associate

  ! ----------------------
  ! Interpolate the solution from nodes to cells.
  ! ----------------------
  !$omp parallel do schedule(static) &
  !$omp & default(private) shared(mesh, u, u_MD_nc)
  do iCell = 1_ip, mesh%NumCells
    associate(iCellMD => mesh%CellMDIndex(:,iCell))
      associate(iNode2D => iCellMD)
        u(iCell) = &
          & ( u_MD_nc(iNode2D(1) + 0, iNode2D(2) + 0) + &
          &   u_MD_nc(iNode2D(1) + 0, iNode2D(2) + 1) + &
          &   u_MD_nc(iNode2D(1) + 1, iNode2D(2) + 0) + &
          &   u_MD_nc(iNode2D(1) + 1, iNode2D(2) + 1) )*0.25_dp
      end associate
    end associate
  end do
  !$omp end parallel do

contains
  
  real(dp) function w_wrapper(x, y)
    ! <<<<<<<<<<<<<<<<<<<<<<
    real(dp), intent(in) :: x, y
    ! >>>>>>>>>>>>>>>>>>>>>>
    w_wrapper = 1.0
  end function w_wrapper

  real(dp) function lambda_wrapper(x, y)
    ! <<<<<<<<<<<<<<<<<<<<<<
    real(dp), intent(in) :: x, y
    ! >>>>>>>>>>>>>>>>>>>>>>
    lambda_wrapper = max(0.0e-8_dp, lambda)
  end function lambda_wrapper

  subroutine boundary(kbdy, xory, alfa, gbdy)
    integer kbdy
    real(dp) xory, alfa, gbdy
  end subroutine boundary
end subroutine Solve_DivWGrad_MUDPACK$rank
#$end do

#$endif

end module StormRuler_Solvers_MUDPACK
