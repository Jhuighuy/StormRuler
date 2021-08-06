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
module StormRuler_KrylovSolvers

#$use 'StormRuler_Parameters.f90'

use StormRuler_Parameters, only: dp
use StormRuler_ConvParams, only: tConvParams
use StormRuler_Mesh, only: tMesh

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

abstract interface
#$do rank = 0, NUM_RANKS
  module subroutine tMeshOperator$rank(mesh, Au, u, opParams)
    class(tMesh), intent(in) :: mesh
    real(dp), intent(in), pointer :: Au(@:,:), u(@:,:)
    class(*), intent(in) :: opParams
  end subroutine tMeshOperator$rank
#$end do
end interface

interface Solve_CG
#$do rank = 0, NUM_RANKS
  module subroutine Solve_CG$rank(mesh, u, b, LOp, opParams, params)
    class(tMesh), intent(in) :: mesh
    real(dp), intent(in), pointer :: u(@:,:), b(@:,:)
    procedure(tMeshOperator$rank) :: LOp
    class(*), intent(in) :: opParams
    type(tConvParams), intent(inout) :: params
  end subroutine Solve_CG$rank
#$end do
end interface Solve_CG

interface Solve_BiCGStab
#$do rank = 0, NUM_RANKS
  module subroutine Solve_BiCGStab$rank(mesh, u, b, LOp, opParams, params)
    class(tMesh), intent(in) :: mesh
    real(dp), intent(in), pointer :: u(@:,:), b(@:,:)
    procedure(tMeshOperator$rank) :: LOp
    class(*), intent(in) :: opParams
    type(tConvParams), intent(inout) :: params
  end subroutine Solve_BiCGStab$rank
#$end do
end interface Solve_BiCGStab

#$if HAS_MKL
interface Solve_CG_MKL
#$do rank = 0, NUM_RANKS
  module subroutine Solve_CG_MKL$rank(mesh, u, b, LOp, opParams, params)
    class(tMesh), intent(in) :: mesh
    real(dp), intent(in), pointer :: u(@:,:), b(@:,:)
    procedure(tMeshOperator$rank) :: LOp
    class(*), intent(in) :: opParams
    type(tConvParams), intent(inout) :: params
  end subroutine Solve_CG_MKL$rank
#$end do
end interface Solve_CG_MKL

interface Solve_FGMRES_MKL
#$do rank = 0, NUM_RANKS
  module subroutine Solve_FGMRES_MKL$rank(mesh, u, b, LOp, opParams, params)
    class(tMesh), intent(in) :: mesh
    real(dp), intent(in), pointer :: u(@:,:), b(@:,:)
    procedure(tMeshOperator$rank) :: LOp
    class(*), intent(in) :: opParams
    type(tConvParams), intent(inout) :: params
  end subroutine Solve_FGMRES_MKL$rank
#$end do
end interface Solve_FGMRES_MKL
#$end if

end module StormRuler_KrylovSolvers
