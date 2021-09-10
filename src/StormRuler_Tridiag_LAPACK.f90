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
module StormRuler_Tridiag_LAPACK

#$use 'StormRuler_Params.fi'
#$if HAS_LAPACK

use StormRuler_Parameters, only: dp, ip
use StormRuler_Tridiag, only: tTridiagMatrix

use, intrinsic :: iso_fortran_env, only: error_unit

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

#$if HAS_MKL
include 'mkl.fi'
#$else
external dsterf, dsteqr, dstein
#$end if

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute all eigenvalues 𝜃ⱼ of the symmetric tridiagonal matrix:
!! 𝜽 = {𝜃ⱼ}, 𝓣𝒚ⱼ = 𝜃ⱼ𝒚ⱼ, using the root-free QR algorithm (from LAPACK).
!! See: https://intel.ly/38Wdyqb
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine ComputeEigenvalues_Symm_LAPACK(T, theta)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tTridiagMatrix), intent(in) :: T
  real(dp), intent(inout) :: theta(:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  real(dp), allocatable :: subdiag(:)
  integer(ip) :: info

  ! ----------------------
  ! Prepair the input arguments.
  ! ----------------------
  theta = T%Diag(:)
  subdiag = T%Subdiag(2:)

  ! ----------------------
  ! Evaluate eigenvalues with LAPACK.
  ! ----------------------
  call dsterf(T%Dim, theta, subdiag, info)
  if (info /= 0) then
    write(error_unit, *) 'LAPACK DSTREF failed, INFO=', info
    error stop 1
  end if

end subroutine ComputeEigenvalues_Symm_LAPACK

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Compute all eigenpairs (𝜃ⱼ,𝒚ⱼ) of the symmetric tridiagonal matrix:
!! 𝒀 = {𝒚ⱼ}, 𝜽 = {𝜃ⱼ}, 𝓣𝒚ⱼ = 𝜃ⱼ𝒚ⱼ, using the QR algorithm (from LAPACK).
!! See: https://intel.ly/37N95pe
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine ComputeEigenpairs_Symm_LAPACK(T, theta, y)
  ! <<<<<<<<<<<<<<<<<<<<<<
  class(tTridiagMatrix), intent(in) :: T
  real(dp), intent(inout) :: theta(:)
  real(dp), intent(inout), optional :: y(:,:)
  ! >>>>>>>>>>>>>>>>>>>>>>

  character(len=1) :: compute_y
  real(dp), allocatable :: subdiag(:), work(:)
  integer(ip) :: dim_y, info

  ! ----------------------
  ! Prepair the input arguments.
  ! ----------------------
  theta = T%Diag(:)
  subdiag = T%Subdiag(2:)
  if (present(y)) then
    compute_y = 'I'; dim_y = T%Dim
    allocate(work(2*T%Dim-2))
  else
    compute_y = 'N'; dim_y = 1
    allocate(work(1))
  end if

  ! ----------------------
  ! Evaluate eigenpairs with LAPACK.
  ! ----------------------
  call dsteqr(compute_y, T%Dim, theta, subdiag, y, dim_y, work, info)
  if (info /= 0) then
    write(error_unit, *) 'LAPACK DSTEQR failed, INFO=', info
    error stop 1
  end if

end subroutine ComputeEigenpairs_Symm_LAPACK

end module StormRuler_Tridiag_LAPACK

#$end if
