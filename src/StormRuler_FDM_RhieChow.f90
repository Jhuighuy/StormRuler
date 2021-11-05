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
module StormRuler_FDM_RhieChow

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip, i8, gCylCoords
use StormRuler_Helpers, only: Flip

use StormRuler_Mesh, only: tMesh
use StormRuler_Array, only: tArrayR

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! The central FDM-approximate divergence with 
!! Rhie-Chow type correction: 𝒗 ← 𝒗 - 𝜆∇⋅𝒖⃗ - 𝜏𝓡𝓒(𝒑,𝛒).
!!
!! Rhie-Chow correction is used to eleminate the checkerboard 
!! pressure phenomenon that may lead to the pressure-velocity 
!! decoupling in the incompressible simulations.
!!
!! Shape of 𝒖⃗ is [1,NumDims]×[1,NumVars]×[1, NumAllCells],
!! shape of 𝒗 is [1,NumVars]×[1, NumAllCells].
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine FDM_Divergence_RhieChow(mesh, vArr, lambda, uVecArr, tau, pArr, rhoAny)
  class(tMesh), intent(inout) :: mesh
  class(tArrayR), intent(in) :: uVecArr, pArr
  class(*), intent(in) :: rhoAny
  class(tArrayR), intent(inout) :: vArr
  real(dp), intent(in) :: lambda, tau

  error stop 'RHIR CHOW!!!'
end subroutine FDM_Divergence_RhieChow

end module StormRuler_FDM_RhieChow
