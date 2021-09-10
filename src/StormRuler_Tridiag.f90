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
module StormRuler_Tridiag

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Real tridiagonal square matrix 𝓣.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type :: tTridiagMatrix
  ! ----------------------
  ! Matrix dimension, 𝘥𝘪𝘮(𝓣).
  ! ----------------------
  integer(ip) :: Dim
  ! ----------------------
  ! Is matrix symmetric, 𝓣 = 𝓣ᵀ.
  ! ----------------------
  logical :: Symmetric

  ! ----------------------
  ! Matrix diagonal, 𝘥𝘪𝘢𝘨(𝓣)ᵢ = 𝓣ᵢ,ᵢ.
  ! Shape is [1, Dim].
  ! ----------------------
  real(dp), pointer :: Diag(:)
  ! ----------------------
  ! Matrix subdiagonal, 𝘴𝘶𝘣𝘥𝘪𝘢𝘨(𝓣)ᵢ = 𝓣ᵢ,ᵢ₋₁.
  ! Shape is [2, Dim].
  ! ----------------------
  real(dp), pointer :: Subdiag(:)
  ! ----------------------
  ! Matrix superdiagonal, 𝘴𝘶𝘱𝘦𝘳𝘥𝘪𝘢𝘨(𝓣)ᵢ = 𝓣ᵢ,ᵢ₊₁.
  ! Shape is [1, Dim-1].
  ! ----------------------
  real(dp), pointer :: Superdiag(:)
end type tTridiagMatrix

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

end module StormRuler_Tridiag
