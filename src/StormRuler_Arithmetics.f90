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
module StormRuler_Arithmetics

use StormRuler_Mesh
use StormRuler_Parameters
#@use 'StormRuler_Parameters.f90'

implicit none

interface Zero
#@do rank = 0, NumRanks
  module procedure Zero$rank
#@end do
end interface Zero

interface Set
#@do rank = 0, NumRanks
  module procedure Set$rank
#@end do
end interface Set

interface Dot
#@do rank = 0, NumRanks
  module procedure Dot$rank
#@end do
end interface Dot

interface Add
#@do rank = 0, NumRanks
  module procedure Add$rank
#@end do
end interface Add

interface Sub
#@do rank = 0, NumRanks
  module procedure Sub$rank
#@end do
end interface Sub

interface Mul
#@do rank = 0, NumRanks
  module procedure Mul$rank
#@end do
end interface Mul

interface ApplyFunc
#@do rank = 0, NumRanks
  module procedure ApplyFunc$rank
#@end do
end interface ApplyFunc

abstract interface
  pure function MathFunc$0(x) result(f)
    import dp
    real(dp), intent(in) :: x
    real(dp) :: f
  end function MathFunc$0
end interface
#@do rank = 1, NumRanks
abstract interface
  pure function MathFunc$rank(x) result(f)
    import dp
    real(dp), intent(in) :: x(@:)
    real(dp) :: f(@{size(x, dim=$$+1)}@)
  end function MathFunc$rank
end interface
#@end do

contains

!! -----------------------------------------------------------------  
!! u = 0
#@do rank = 0, NumRanks
subroutine Zero$rank(mesh, u)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(inout) :: u(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell
  !$omp parallel do private(iCell)
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = 0
  end do
  !$omp end parallel do
end subroutine Zero$rank
#@end do

!! -----------------------------------------------------------------  
!! u = v
#@do rank = 0, NumRanks
subroutine Set$rank(mesh, u,v)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(inout) :: u(@:,:)
  real(dp), intent(in) :: v(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell
  !$omp parallel do private(iCell)
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = v(@:,iCell)
  end do
  !$omp end parallel do
end subroutine Set$rank
#@end do

!! -----------------------------------------------------------------  
!! Compute a dot product.
#@do rank = 0, NumRanks
function Dot$rank(mesh, u,v) result(d)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: u(@:,:), v(@:,:)
  real(dp) :: d
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell
  d = 0.0_dp
  associate(dv=>product(mesh%Dx))
    !$omp parallel do private(iCell) reduction(+:d)
    do iCell = 1, mesh%NumCells
#@if rank == 0
      d += dv*u(iCell)*v(iCell)
#@else if rank == 1
      d += dv*dot_product(u(:,iCell),v(:,iCell))
#@else
      d += dv*sum(u(@:,iCell)*v(@:,iCell))
#@end if
    end do
    !$omp end parallel do
  end associate
end function Dot$rank
#@end do

!! -----------------------------------------------------------------  
!! Compute "u ← c⋅v + a⋅w".
#@do rank = 0, NumRanks
subroutine Add$rank(mesh, u,v,w,a,c)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(out) :: u(@:,:)
  real(dp), intent(inout) :: v(@:,:), w(@:,:)
  real(dp), intent(in), optional :: a, c
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell
  real(dp) :: aa, cc
  aa = 1.0_dp; if (present(a)) aa = a 
  cc = 1.0_dp; if (present(c)) cc = c 
  !$omp parallel do private(iCell)
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = cc*v(@:,iCell) + aa*w(@:,iCell)
  end do
  !$omp end parallel do
end subroutine Add$rank
#@end do

!! -----------------------------------------------------------------  
!! Compute "u ← c⋅v - a⋅w".
#@do rank = 0, NumRanks
subroutine Sub$rank(mesh, u,v,w, a,c)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(out) :: u(@:,:)
  real(dp), intent(inout) :: v(@:,:), w(@:,:)
  real(dp), intent(in), optional :: a, c
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell
  real(dp) :: aa, cc
  aa = 1.0_dp; if (present(a)) aa = a 
  cc = 1.0_dp; if (present(c)) cc = c 
  !$omp parallel do private(iCell)
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = cc*v(@:,iCell) - aa*w(@:,iCell)
  end do
  !$omp end parallel do
end subroutine Sub$rank
#@end do

!! -----------------------------------------------------------------  
!! Compute "u ← v⋅w".
#@do rank = 0, NumRanks
subroutine Mul$rank(mesh, u,v,w)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(out) :: u(@:,:)
  real(dp), intent(inout) :: v(:), w(@:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell
  !$omp parallel do private(iCell)
  do iCell = 1, mesh%NumCells
    u(@:,iCell) = v(iCell)*w(@:,iCell)
  end do
  !$omp end parallel do
end subroutine Mul$rank
#@end do

!! -----------------------------------------------------------------  
!! Apply a function.
#@do rank = 0, NumRanks
subroutine ApplyFunc$rank(mesh, Fu,u, f)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(inout) :: u(@:,:), Fu(@:,:)
  procedure(MathFunc$rank) :: f
  ! >>>>>>>>>>>>>>>>>>>>>>
  integer :: iCell
  !$omp parallel do private(iCell)
  do iCell = 1, mesh%NumCells
    Fu(@:,iCell) = f(u(@:,iCell))
  end do
  !$omp end parallel do
end subroutine ApplyFunc$rank
#@end do

end module StormRuler_Arithmetics
