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

implicit none

interface Zero
  module procedure Zero1,Zero2
end interface Zero
interface Set
  module procedure Set1,Set2
end interface Set
interface Add
  module procedure Add1,Add2
end interface Add
interface Sub
  module procedure Sub1,Sub2
end interface Sub
interface Mul
  module procedure Mul2
end interface Mul
interface Dot
  module procedure Dot1,Dot2
end interface Dot

contains

!! -----------------------------------------------------------------  
!! Y = X
subroutine Zero1(mesh, Y)
  class(Mesh2D), intent(in) :: mesh
  real(8), dimension(:), intent(inout) :: Y
  block
    integer :: iCell
    !$omp parallel do private(iCell)
    do iCell = 1, mesh%NumCells
      Y(iCell) = 0
    end do
    !$omp end parallel do
  end block
end subroutine Zero1
subroutine Zero2(mesh, Y)
  class(Mesh2D), intent(in) :: mesh
  real(8), dimension(:,:), intent(inout) :: Y
  block
    integer :: iCell
    !$omp parallel do private(iCell)
    do iCell = 1, mesh%NumCells
      Y(:,iCell) = 0
    end do
    !$omp end parallel do
  end block
end subroutine Zero2

!! -----------------------------------------------------------------  
!! Y = X
subroutine Set1(mesh, Y,X)
  class(Mesh2D), intent(in) :: mesh
  real(8), dimension(:), intent(inout) :: Y
  real(8), dimension(:), intent(in) :: X
  block
    integer :: iCell
    !$omp parallel do private(iCell)
    do iCell = 1, mesh%NumCells
      Y(iCell) = X(iCell)
    end do
    !$omp end parallel do
  end block
end subroutine Set1
subroutine Set2(mesh, Y,X)
  class(Mesh2D), intent(in) :: mesh
  real(8), dimension(:,:), intent(inout) :: Y
  real(8), dimension(:,:), intent(in) :: X
  block
    integer :: iCell
    !$omp parallel do private(iCell)
    do iCell = 1, mesh%NumCells
      Y(:,iCell) = X(:,iCell)
    end do
    !$omp end parallel do
  end block
end subroutine Set2

!! -----------------------------------------------------------------  
!! Compute a dot product.
function Dot1(mesh, u,v) result(d)
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: u(:), v(:)
  real(dp) :: d
  block
    integer :: iCell
    d = 0
    !$omp parallel do private(iCell) reduction(+:d)
    do iCell = 1, mesh%NumCells
      d = d + (mesh%Dx(1)**2)*u(iCell)*v(iCell)
    end do
    !$omp end parallel do
  end block
end function Dot1
function Dot2(mesh, u,v) result(d)
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(in) :: u(:,:), v(:,:)
  real(dp) :: d
  block
    integer :: iCell
    d = 0
    !$omp parallel do private(iCell) reduction(+:d)
    do iCell = 1, mesh%NumCells
      d = d + (mesh%Dx(1)**2)*dot_product(u(:,iCell),v(:,iCell))
    end do
    !$omp end parallel do
  end block
end function Dot2

!! -----------------------------------------------------------------  
!! Compute "u = c⋅v + a⋅w".
subroutine Add1(mesh, u,v,w,a,c)
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(out) :: u(:)
  real(dp), intent(inout) :: v(:), w(:)
  real(dp), intent(in), optional :: a, c
  real(dp) :: b, d
  b = 1.0_dp; if (present(a)) b = a 
  d = 1.0_dp; if (present(c)) d = c 
  block
    integer :: iCell
    !$omp parallel do private(iCell)
    do iCell = 1, mesh%NumCells
      u(iCell) = d*v(iCell) + b*w(iCell)
    end do
    !$omp end parallel do
  end block
end subroutine Add1
subroutine Add2(mesh, u,v,w,a,c)
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(out) :: u(:,:)
  real(dp), intent(inout) :: v(:,:), w(:,:)
  real(dp), intent(in), optional :: a, c
  real(dp) :: b, d
  b = 1.0_dp; if (present(a)) b = a 
  d = 1.0_dp; if (present(c)) d = c 
  block
    integer :: iCell
    !$omp parallel do private(iCell)
    do iCell = 1, mesh%NumCells
      u(:,iCell) = d*v(:,iCell) + b*w(:,iCell)
    end do
    !$omp end parallel do
  end block
end subroutine Add2

!! -----------------------------------------------------------------  
!! Compute "u = c⋅v - a⋅w".
subroutine Sub1(mesh, u,v,w, a,c)
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(out) :: u(:)
  real(dp), intent(inout) :: v(:), w(:)
  real(dp), intent(in), optional :: a, c
  real(dp) :: b, d
  b = 1.0_dp; if (present(a)) b = a 
  d = 1.0_dp; if (present(c)) d = c 
  block
    integer :: iCell
    !$omp parallel do private(iCell)
    do iCell = 1, mesh%NumCells
      u(iCell) = d*v(iCell) - b*w(iCell)
    end do
    !$omp end parallel do
  end block
end subroutine Sub1
subroutine Sub2(mesh, u,v,w,c,a)
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(out) :: u(:,:)
  real(dp), intent(inout) :: v(:,:), w(:,:)
  real(dp), intent(in), optional :: a, c
  real(dp) :: b, d
  b = 1.0_dp; if (present(a)) b = a 
  d = 1.0_dp; if (present(c)) d = c 
  block
    integer :: iCell
    !$omp parallel do private(iCell)
    do iCell = 1, mesh%NumCells
      u(:,iCell) = d*v(:,iCell) - b*w(:,iCell)
    end do
    !$omp end parallel do
  end block
end subroutine Sub2

subroutine Mul2(mesh, u,v,w)
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(out) :: u(:,:)
  real(dp), intent(inout) :: v(:), w(:,:)
  block
    integer :: iCell
    !$omp parallel do private(iCell)
    do iCell = 1, mesh%NumCells
      u(:,iCell) = v(iCell)*w(:,iCell)
    end do
    !$omp end parallel do
  end block
end subroutine Mul2

!! -----------------------------------------------------------------  
!! Apply a function.
subroutine ApplyFunc(mesh, Fu,u,f)
  class(Mesh2D), intent(in) :: mesh
  real(dp), intent(inout) :: u(:), Fu(:)
  procedure(MathFunc) :: f
  block
    integer :: iCell
    !$omp parallel do private(iCell)
    do iCell = 1, mesh%NumCells
      Fu(iCell) = f(u(iCell))
    end do
    !$omp end parallel do
  end block
end subroutine ApplyFunc

end module StormRuler_Arithmetics
