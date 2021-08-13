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
module StormRuler_BLAS

#$use 'StormRuler_Parameters.f90'

use StormRuler_Parameters, only: dp, ip
use StormRuler_Helpers, only: &
  & @{tMapFunc$$, tSMapFunc$$@|@0, NUM_RANKS}@, &
  & operator(.inner.), operator(.outer.)
use StormRuler_Mesh, only: tMesh

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

interface Fill
#$do rank = 0, NUM_RANKS
  module procedure Fill$rank
#$end do
end interface Fill

interface Set
#$do rank = 0, NUM_RANKS
  module procedure Set$rank
#$end do
end interface Set

interface Dot
#$do rank = 0, NUM_RANKS
  module procedure Dot$rank
#$end do
end interface Dot

interface Add
#$do rank = 0, NUM_RANKS
  module procedure Add$rank
#$end do
end interface Add

interface Sub
#$do rank = 0, NUM_RANKS
  module procedure Sub$rank
#$end do
end interface Sub

interface Mul
#$do rank = 0, NUM_RANKS
  module procedure Mul$rank
#$end do
end interface Mul

interface Mul_Inner
#$do rank = 0, NUM_RANKS-1
  module procedure Mul_Inner$rank
#$end do
end interface Mul_Inner

interface Mul_Outer
#$do rank = 0, NUM_RANKS-1
  module procedure Mul_Outer$rank
#$end do
end interface Mul_Outer

interface FuncProd
#$do rank = 0, NUM_RANKS
  module procedure FuncProd$rank
#$end do 
end interface FuncProd

interface SFuncProd
#$do rank = 0, NUM_RANKS
  module procedure SFuncProd$rank
#$end do 
end interface SFuncProd

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Fill vector components: u ← α.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine Fill$rank(mesh, u, alpha)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: alpha
  real(dp), intent(inout) :: u(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer(ip) :: iCell
  
  ! ----------------------
  !$omp parallel do schedule(static) 
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = alpha
  end do
  !$omp end parallel do
end subroutine Fill$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! u ← v
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine Set$rank(mesh, u, v)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: v(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: iCell
  
  ! ----------------------
  !$omp parallel do schedule(static) &
  !$omp default(private) shared(mesh, u, v)
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = v(@:,iCell)
  end do
  !$omp end parallel do
end subroutine Set$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute dot product: d ← <u⋅v>.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
function Dot$rank(mesh, u, v) result(d)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: u(@:,:), v(@:,:)
  real(dp) :: d
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: iCell
  
  d = 0.0_dp
  ! ----------------------
  !$omp parallel do reduction(+:d) schedule(static) &
  !$omp & default(private) shared(mesh, u, v)
  do iCell = 1, mesh%NumCells
#$if rank == 0
    d = d + u(iCell)*v(iCell)
#$else if rank == 1
    d = d + dot_product(u(:,iCell), v(:,iCell))
#$else
    d = d + sum(u(@:,iCell)*v(@:,iCell))
#$end if
  end do
  !$omp end parallel do
end function Dot$rank
#$end do

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute linear combination: u ← βv + αw.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine Add$rank(mesh, u, v, w, alpha, beta)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: v(@:,:), w(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  real(dp), intent(in), optional :: alpha, beta
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: iCell
  real(dp) :: a, b
  a = 1.0_dp; if (present(alpha)) a = alpha
  b = 1.0_dp; if (present(beta)) b = beta
  
  ! ----------------------
  !$omp parallel do schedule(static) &
  !$omp & default(private) shared(mesh, a, b, u, v, w)
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = b*v(@:,iCell) + a*w(@:,iCell)
  end do
  !$omp end parallel do
end subroutine Add$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute linear combination: u ← βv - αw.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine Sub$rank(mesh, u, v, w, alpha, beta)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: v(@:,:), w(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  real(dp), intent(in), optional :: alpha, beta
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: iCell
  real(dp) :: a, b
  a = 1.0_dp; if (present(alpha)) a = alpha
  b = 1.0_dp; if (present(beta)) b = beta
  
  ! ----------------------
  !$omp parallel do schedule(static) &
  !$omp & default(private) shared(mesh, a, b, u, v, w)
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = b*v(@:,iCell) - a*w(@:,iCell)
  end do
  !$omp end parallel do
end subroutine Sub$rank
#$end do

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute product: u̅ ← vw̅.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine Mul$rank(mesh, u, v, w, power)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: v(:), w(@:,:)
  real(dp), intent(inout) :: u(@:,:)
  integer(ip), intent(in), optional :: power
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: iCell, p
  p = 1; if (present(power)) p = power
  
  ! ----------------------
  !$omp parallel do schedule(static) &
  !$omp & default(private) shared(mesh, u, v, w, p)
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = (v(iCell)**p)*w(@:,iCell)
  end do
  !$omp end parallel do
end subroutine Mul$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute an inner product: u ← v̅⋅w̅.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS-1
subroutine Mul_Inner$rank(mesh, u, vBar, wBar)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: vBar(:,:), wBar(:,@:,:)
  real(dp), intent(inout) :: u(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: iCell
  
  ! ----------------------
  !$omp parallel do schedule(static) &
  !$omp & default(private) shared(mesh, u, vBar, wBar)
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = vBar(:,iCell).inner.wBar(:,@:,iCell)
  end do
  !$omp end parallel do
end subroutine Mul_Inner$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute an outer product: û ← v̅⊗w̅.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS-1
subroutine Mul_Outer$rank(mesh, uHat, vBar, wBar)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: vBar(:,:), wBar(@:,:)
  real(dp), intent(inout) :: uHat(:,@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: iCell
  
  ! ----------------------
  !$omp parallel do schedule(static) &
  !$omp & default(private) shared(mesh, uHat, vBar, wBar)
  do iCell = 1, mesh%NumCells
    uHat(:,@:,iCell) = vBar(:,iCell).outer.wBar(@:,iCell)
  end do
  !$omp end parallel do
end subroutine Mul_Outer$rank
#$end do

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a function product: v ← f(u).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine FuncProd$rank(mesh, v, u, f)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: u(@:,:)
  real(dp), intent(inout) :: v(@:,:)
  procedure(tMapFunc$rank) :: f
  ! >>>>>>>>>>>>>>>>>>>>>>
  
  integer(ip) :: iCell
  
  ! ----------------------
#$if not NAG_COMPILER
  !$omp parallel do schedule(static) &
  !$omp & default(private) shared(mesh, u, v)
#$end if
  do iCell = 1, mesh%NumCells
    v(@:,iCell) = f(u(@:,iCell))
  end do
#$if not NAG_COMPILER
  !$omp end parallel do
#$end if
end subroutine FuncProd$rank
#$end do

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute a function product: v ← f(u).
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
#$do rank = 0, NUM_RANKS
subroutine SFuncProd$rank(mesh, v, u, f)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(in) :: mesh
  real(dp), intent(in) :: u(@:,:)
  real(dp), intent(inout) :: v(@:,:)
  procedure(tSMapFunc$rank) :: f
  ! >>>>>>>>>>>>>>>>>>>>>>

  integer(ip) :: iCell

  ! ----------------------
#$if not NAG_COMPILER
  !$omp parallel do schedule(static) &
  !$omp & default(private) shared(mesh, u, v)
#$end if
  do iCell = 1, mesh%NumCells
    v(@:,iCell) = f(mesh%dl(::2)*mesh%CellMDIndex(:,iCell), u(@:,iCell))
  end do
#$if not NAG_COMPILER
  !$omp end parallel do
#$end if
end subroutine SFuncProd$rank
#$end do

end module StormRuler_BLAS
