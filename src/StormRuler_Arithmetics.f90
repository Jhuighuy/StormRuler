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
module StormRuler_Arithmetics

use StormRuler_Mesh
use StormRuler_Parameters
#@use 'StormRuler_Parameters.f90'

implicit none

interface Fill
#@do rank = 0, NUM_RANKS
  module procedure Fill$rank
#@end do
end interface Fill

interface Set
#@do rank = 0, NUM_RANKS
  module procedure Set$rank
#@end do
end interface Set

interface Dot
#@do rank = 0, NUM_RANKS
  module procedure Dot$rank
#@end do
end interface Dot

interface Add
#@do rank = 0, NUM_RANKS
  module procedure Add$rank
#@end do
end interface Add

interface Sub
#@do rank = 0, NUM_RANKS
  module procedure Sub$rank
#@end do
end interface Sub

interface Mul
#@do rank = 0, NUM_RANKS
  module procedure Mul$rank
#@end do
end interface Mul

interface Mul_Inner
#@do rank = 0, NUM_RANKS-1
  module procedure Mul_Inner$rank
#@end do
end interface Mul_Inner

interface Mul_Outer
#@do rank = 0, NUM_RANKS-1
  module procedure Mul_Outer$rank
#@end do
end interface Mul_Outer

abstract interface
  pure function MathFunc$0(x) result(f)
    import dp
    real(dp), intent(in) :: x
    real(dp) :: f
  end function MathFunc$0
end interface
#@do rank = 1, NUM_RANKS
abstract interface
  pure function MathFunc$rank(x) result(f)
    import dp
    real(dp), intent(in) :: x(@:)
    real(dp) :: f(@{size(x,dim=$$+1)}@)
  end function MathFunc$rank
end interface
#@end do

interface ApplyFunc
#@do rank = 0, NUM_RANKS
  module procedure ApplyFunc$rank
#@end do
end interface ApplyFunc

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
contains
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Fill vector components: u ← α.
#@do rank = 0, NUM_RANKS
subroutine Fill$rank(mesh,u,alpha)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in), optional :: alpha
  real(dp), intent(out) :: u(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell
  real(dp) :: a
  a = 0.0_dp; if (present(alpha)) a = alpha
  !$omp parallel do
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = a
  end do
  !$omp end parallel do
end subroutine Fill$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! u = v
#@do rank = 0, NUM_RANKS
subroutine Set$rank(mesh,u,v)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: v(@:,:)
  real(dp), intent(out) :: u(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell
  !$omp parallel do
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = v(@:,iCell)
  end do
  !$omp end parallel do
end subroutine Set$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute dot product: d ← <u⋅v>.
#@do rank = 0, NUM_RANKS
function Dot$rank(mesh,u,v) result(d)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: u(@:,:), v(@:,:)
  real(dp) :: d
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell
  d = 0.0_dp
  associate(dv=>product(mesh%Dx))
    !$omp parallel do reduction(+:d)
    do iCell = 1, mesh%NumCells
#@if rank == 0
      d += dv*u(iCell)*v(iCell)
#@else
      d += dv*sum(u(@:,iCell)*v(@:,iCell))
#@end if
    end do
    !$omp end parallel do
  end associate
end function Dot$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute linear combination: u ← βv + αw.
#@do rank = 0, NUM_RANKS
subroutine Add$rank(mesh,u,v,w,alpha,beta)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: v(@:,:), w(@:,:)
  real(dp), intent(in), optional :: alpha, beta
  real(dp), intent(out) :: u(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell
  real(dp) :: a, b
  a = 1.0_dp; if (present(alpha)) a = alpha
  b = 1.0_dp; if (present(beta)) b = beta 
  !$omp parallel do
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = b*v(@:,iCell) + a*w(@:,iCell)
  end do
  !$omp end parallel do
end subroutine Add$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute linear combination: u ← βv - αw.
#@do rank = 0, NUM_RANKS
subroutine Sub$rank(mesh,u,v,w,alpha,beta)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: v(@:,:), w(@:,:)
  real(dp), intent(in), optional :: alpha, beta
  real(dp), intent(out) :: u(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell
  real(dp) :: a, b
  a = 1.0_dp; if (present(alpha)) a = alpha 
  b = 1.0_dp; if (present(beta)) b = beta 
  !$omp parallel do
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = b*v(@:,iCell) - a*w(@:,iCell)
  end do
  !$omp end parallel do
end subroutine Sub$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute "u ← vw".
#@do rank = 0, NUM_RANKS
subroutine Mul$rank(mesh,u,v,w)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: v(:), w(@:,:)
  real(dp), intent(out) :: u(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell
  !$omp parallel do
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = v(iCell)*w(@:,iCell)
  end do
  !$omp end parallel do
end subroutine Mul$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute "u ← v⋅w".
#@do rank = 0, NUM_RANKS-1
subroutine Mul_Inner$rank(mesh,u,v,w)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: v(:,:),w(:,@:,:)
  real(dp), intent(out) :: u(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell
  !$omp parallel do
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = Inner(v(:,iCell),w(:,@:,iCell))
  end do
  !$omp end parallel do
end subroutine Mul_Inner$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute "u ← v⊗w".
#@do rank = 0, NUM_RANKS-1
subroutine Mul_Outer$rank(mesh,u,v,w)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: v(:,:),w(@:,:)
  real(dp), intent(out) :: u(:,@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell
  !$omp parallel do
  do iCell = 1, mesh%NumCells
    u(:,@:,iCell) = Outer(v(:,iCell),w(@:,iCell))
  end do
  !$omp end parallel do
end subroutine Mul_Outer$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Apply alpha function.
#@do rank = 0, NUM_RANKS
subroutine ApplyFunc$rank(mesh,Fu,u,f)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: u(@:,:)
  real(dp), intent(out) :: Fu(@:,:)
  procedure(MathFunc$rank) :: f
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell
  !$omp parallel do
  do iCell = 1, mesh%NumCells
    Fu(@:,iCell) = f(u(@:,iCell))
  end do
  !$omp end parallel do
end subroutine ApplyFunc$rank
#@end do
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!

end module StormRuler_Arithmetics
