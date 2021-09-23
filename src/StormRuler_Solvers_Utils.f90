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
module StormRuler_Solvers_Utils

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp
use StormRuler_Mesh, only: tMesh
use StormRuler_BLAS, only: Set, Sub, Norm_2
#$for type_, _ in SCALAR_TYPES
use StormRuler_BLAS, only: tMatVecFunc$type_$1
use StormRuler_Solvers_Precond, only: tPrecondFunc$type_
#$end for

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Compute squared residual norm: ‖𝒃 - 𝓐[𝓟]𝒚‖₂.
!! ----------------------------------------------------------------- !!
real(dp) function ResidualNorm(mesh, y, b, MatVec, env, Precond, precond_env)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(:,:), y(:,:)
  procedure(tMatVecFuncR$1) :: MatVec
  class(*), intent(inout) :: env
  procedure(tPrecondFuncR), optional :: Precond
  class(*), allocatable, intent(inout), optional :: precond_env
  ! >>>>>>>>>>>>>>>>>>>>>>

  real(dp), pointer :: r(:,:), x(:,:)
  allocate(r, x, mold=y)

  ! ----------------------
  ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  !   𝒙 ← 𝓟𝒚, 𝒓 ← 𝓐𝒙, 
  ! 𝗲𝗹𝘀𝗲: 𝒓 ← 𝓐𝒚, 𝗲𝗻𝗱 𝗶𝗳 
  ! 𝒓 ← 𝒃 - 𝒓,
  ! 𝗼𝘂𝘁 ← ‖𝒓‖₂.
  ! ----------------------
  if (present(Precond)) then
    call Precond(mesh, x, y, MatVec, env, precond_env)
    call MatVec(mesh, r, x, env)
  else
    call MatVec(mesh, r, y, env)
  end if
  call Sub(mesh, r, b, r)
  ResidualNorm = Norm_2(mesh, r)

end function ResidualNorm

!! ----------------------------------------------------------------- !!
!! Compute squared residual norm: ‖(𝓐[𝓟])ᵀ(𝒃 - 𝓐[𝓟]𝒚)‖₂.
!! ----------------------------------------------------------------- !!
real(dp) function ResidualNorm_Squared(mesh, y, b, &
    & MatVec, env, MatVec_T, env_T, Precond, precond_env, Precond_T, precond_env_T)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: b(:,:), y(:,:)
  procedure(tMatVecFuncR$1) :: MatVec, MatVec_T
  class(*), intent(inout) :: env, env_T
  procedure(tPrecondFuncR), optional :: Precond, Precond_T
  class(*), allocatable, intent(inout), optional :: precond_env, precond_env_T
  ! >>>>>>>>>>>>>>>>>>>>>>

  real(dp), pointer :: r(:,:), x(:,:)
  allocate(r, x, mold=y)

  ! ----------------------
  ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  !   𝒙 ← 𝓟𝒚, 𝒓 ← 𝓐𝒙, 
  ! 𝗲𝗹𝘀𝗲: 𝒓 ← 𝓐𝒚, 𝗲𝗻𝗱 𝗶𝗳 
  ! 𝒓 ← 𝒃 - 𝒓.
  ! ----------------------
  if (present(Precond)) then
    call Precond(mesh, x, y, MatVec, env, precond_env)
    call MatVec(mesh, r, x, env)
  else
    call MatVec(mesh, r, y, env)
  end if
  call Sub(mesh, r, b, r)

  ! ----------------------
  ! 𝗶𝗳 𝓟 ≠ 𝗻𝗼𝗻𝗲:
  !   𝒙 ← 𝓐ᵀ𝒓, 𝒓 ← 𝓟ᵀ𝒙, 
  ! 𝗲𝗹𝘀𝗲: 𝒙 ← 𝒓, 𝒓 ← 𝓐ᵀ𝒙, 𝗲𝗻𝗱 𝗶𝗳
  ! 𝗼𝘂𝘁 ← ‖𝒓‖₂.
  ! ----------------------
  if (present(Precond)) then
    call MatVec_T(mesh, x, r, env_T)
    call Precond_T(mesh, r, x, MatVec_T, env_T, precond_env_T)
  else
    call Set(mesh, x, r)
    call MatVec_T(mesh, r, x, env_T)
  end if
  ResidualNorm_Squared = Norm_2(mesh, r)

end function ResidualNorm_Squared

end module StormRuler_Solvers_Utils
