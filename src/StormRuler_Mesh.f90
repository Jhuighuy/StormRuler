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
module StormRuler_Mesh

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip, i8

use StormRuler_Helpers, only: ErrorStop, PrintWarning, &
  & Flip, IndexOf, InsertTo, I2S, R2S, IntToPixel
use StormRuler_IO, only: IOList, IOListItem, @{IOListItem$$@|@0, 2}@

use, intrinsic :: iso_fortran_env, only: error_unit

#$if HAS_OpenMP
use :: omp_lib
#$endif

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Semi-structured multidimensional mesh.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !! 
type :: tMesh
  ! ----------------------
  ! Number of spatial dimension.
  ! ----------------------
  integer(ip) :: NumDims

  ! ----------------------
  ! Number of faces (edges in 2D) per each cell.
  ! NumCellFaces = 4 for 2D meshes, NumCellFaces = 6 for 3D meshes.
  ! ----------------------
  integer(ip) :: NumCellFaces
  ! ----------------------
  ! Number of the connection per each cell.
  ! Typically, NumCellConns = NumCellFaces.
  ! For LBM simulations, extended connectivity is required, here 
  ! NumCellConns = 9 for 2D meshes, NumCellConns = 19 for 3D meshes (D2Q9 and D3Q19).
  ! ----------------------
  integer(ip) :: NumCellConns

  ! ----------------------
  ! Number of interior cells.
  ! This value should be used for field operations.
  ! ----------------------
  integer(ip) :: NumCells
  ! ----------------------
  ! Number of boundary cells 
  ! (domain-exterior cells, first in the ghost layers).
  ! Cells may be counted multiple times..
  ! ----------------------
  integer(ip) :: NumBndCells
  ! ----------------------
  ! Total number of cells (including interior, boundary and ghost cells).
  ! This value should be used for field allocation.
  ! ----------------------
  integer(ip) :: NumAllCells
  ! ----------------------
  ! Cell connectivity table.
  ! Shape is [1, NumCellConns]×[1, NumAllCells].
  ! ----------------------
  integer(ip), allocatable :: CellToCell(:,:)
  ! ----------------------
  ! Logical table for the periodic face connections.
  ! Shape is [1, NumCellConns]×[1, NumAllCells].
  ! ----------------------
  logical, allocatable, private :: mIsCellConnPeriodic(:,:)

  ! ----------------------
  ! Multidimensional index bounds table.
  ! Shape is [1, NumDims].
  ! ----------------------
  integer(ip), allocatable :: CellIndexBounds(:)
  ! ----------------------
  ! Cell multidimensional index table.
  ! Shape is [1, NumDims]×[1, NumAllCells].
  ! ----------------------
  integer(ip), allocatable :: CellIndex(:,:)

  ! ----------------------
  ! Number of boundary condition marks.
  ! ----------------------
  integer(ip) :: NumBndMarks
  ! ----------------------
  ! Addresses of the boundary cell indices per each mark.
  ! Shape is [1, NumBndMarks+1].
  ! ----------------------
  integer(ip), allocatable :: BndCellAddrs(:)
  ! ----------------------
  ! Indices of the boundary cell connections per each mark.
  ! Shape is [BndCellAddrs(1), BndCellAddrs(NumBndMarks + 1) - 1].
  ! ----------------------
  integer(i8), allocatable :: BndCellConns(:)

  ! ----------------------
  ! Distance between centers of the adjacent cells per face.
  ! Shape is [1, NumCellConns].
  ! ----------------------
  real(dp), allocatable :: dl(:)
  ! ----------------------
  ! Difference between centers of the adjacent cells per face.
  ! Shape is [1, NumDims]×[1, NumCellConns].
  ! ----------------------
  real(dp), allocatable :: dr(:,:)

contains

  ! ----------------------
  ! Initializers.
  ! ----------------------

  ! ----------------------
  ! Parallel mesh walkthough subroutines.
  ! ----------------------
  procedure :: RunCellKernel => RunMeshCellSimpleKernel
  procedure :: RunCellKernel_Block => RunMeshBlockCellKernel
  procedure :: RunCellKernel_Sum => RunMeshSumCellKernel
  procedure :: RunCellKernel_Min => RunMeshMinCellKernel
  procedure :: RunCellKernel_Max => RunMeshMaxCellKernel

  ! ----------------------
  ! Field wrappers.
  ! ----------------------
  generic :: CellCenter => CellCenterVec, CellCenterCoord
  procedure :: CellCenterVec => GetMeshCellCenterVec
  procedure :: CellCenterCoord => GetMeshCellCenterCoord
  procedure :: IsCellConnPeriodic => IsMeshCellConnPeriodic

  ! ----------------------
  ! Ordering routines. 
  ! ----------------------
  procedure :: ApplyOrdering => ApplyMeshOrdering

  ! ----------------------
  ! Printers.
  ! ----------------------
  procedure :: PrintTo_Neato => PrintMeshToNeato
  procedure :: PrintTo_LegacyVTK => PrintMeshToLegacyVTK
  
end type tMesh

abstract interface
  subroutine tKernelFunc(cell)
    import ip
    integer(ip), intent(in) :: cell
  end subroutine tKernelFunc
end interface

abstract interface
  subroutine tBlockKernelFunc(firstCell, lastCell)
    import ip
    integer(ip), intent(in) :: firstCell, lastCell
  end subroutine tBlockKernelFunc
end interface

abstract interface
  function tReduceKernelFunc(cell) result(r)
    import ip, dp
    integer(ip), intent(in) :: cell
    real(dp) :: r
  end function tReduceKernelFunc
end interface

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Launch a cell kernel.
!! ----------------------------------------------------------------- !!
subroutine RunMeshCellSimpleKernel(mesh, Kernel)
  class(tMesh), intent(in) :: mesh
  procedure(tKernelFunc) :: Kernel

  integer :: cell

  !$omp parallel do default(none) shared(mesh) schedule(static)
  do cell = 1, mesh%NumCells
    call Kernel(cell)
  end do
  !$omp end parallel do

end subroutine RunMeshCellSimpleKernel

!! ----------------------------------------------------------------- !!
!! Launch a block-kernel.
!! ----------------------------------------------------------------- !!
subroutine RunMeshBlockCellKernel(mesh, BlockKernel)
  class(tMesh), intent(in) :: mesh
  procedure(tBlockKernelFunc) :: BlockKernel

#$if not HAS_OpenMP

  call BlockKernel(1, mesh%NumCells)

#$else

  integer(ip) :: i, thread
  integer(ip), allocatable :: ranges(:)

  ! ----------------------
  ! Compute block range for each thread.
  ! ----------------------
  allocate(ranges(omp_get_max_threads()))
  associate(size => mesh%NumCells)
    ranges(:) = size/omp_get_max_threads()
    associate(remainder => ranges(:mod(size, omp_get_max_threads())))
      remainder(:) = remainder(:) + 1
    end associate
  end associate
  ranges = [1, ranges]
  do i = 1, omp_get_max_threads()
    ranges(i + 1) = ranges(i) + ranges(i + 1) 
  end do

  ! ----------------------
  ! Lauch threads.
  ! ----------------------
  !$omp parallel default(none) shared(ranges) private(thread)
  thread = omp_get_thread_num() + 1
  call BlockKernel(ranges(thread), ranges(thread + 1) - 1)
  !$omp end parallel

#$endif
end subroutine RunMeshBlockCellKernel

!! ----------------------------------------------------------------- !!
!! Launch a SUM-reduction cell kernel.
!! ----------------------------------------------------------------- !!
function RunMeshSumCellKernel(mesh, Kernel) result(sum)
  class(tMesh), intent(in) :: mesh
  procedure(tReduceKernelFunc) :: Kernel
  real(dp) :: sum

  integer :: cell

  sum = 0.0_dp
  !$omp parallel do default(none) shared(mesh) &
  !$omp & schedule(static) reduction(+:sum)
  do cell = 1, mesh%NumCells
    sum = sum + Kernel(cell)
  end do
  !$omp end parallel do

end function RunMeshSumCellKernel

!! ----------------------------------------------------------------- !!
!! Launch a min-reduction cell kernel.
!! ----------------------------------------------------------------- !!
function RunMeshMinCellKernel(mesh, Kernel) result(minValue)
  class(tMesh), intent(in) :: mesh
  procedure(tReduceKernelFunc) :: Kernel
  real(dp) :: minValue

  integer :: cell

  minValue = +huge(minValue)
  !$omp parallel do default(none) shared(mesh) &
  !$omp & schedule(static) reduction(min:minValue)
  do cell = 1, mesh%NumCells
    minValue = min(minValue, Kernel(cell))
  end do
  !$omp end parallel do

end function RunMeshMinCellKernel

!! ----------------------------------------------------------------- !!
!! Launch a max-reduction cell kernel.
!! ----------------------------------------------------------------- !!
function RunMeshMaxCellKernel(mesh, Kernel) result(maxValue)
  class(tMesh), intent(in) :: mesh
  procedure(tReduceKernelFunc) :: Kernel
  real(dp) :: maxValue

  integer :: cell

  maxValue = -huge(maxValue)
  !$omp parallel do default(none) shared(mesh) &
  !$omp & schedule(static) reduction(max:maxValue)
  do cell = 1, mesh%NumCells
    maxValue = max(maxValue, Kernel(cell))
  end do
  !$omp end parallel do

end function RunMeshMaxCellKernel

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! ----------------------------------------------------------------- !!
!! Get cell center.
!! ----------------------------------------------------------------- !!
pure function GetMeshCellCenterVec(mesh, cell) result(cellCenter)
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: cell
  real(dp) :: cellCenter(mesh%NumDims)

  !! TODO: add initial offset as mesh class field.
  associate(cellIndex => mesh%CellIndex(:,cell))
    cellCenter = [0.0_dp,0.0_dp] + mesh%dl(:2*mesh%NumDims:2)*(cellIndex - 0.5_dp)
  end associate

end function GetMeshCellCenterVec
pure function GetMeshCellCenterCoord(mesh, dim, cell) result(cellCenterCoord)
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: dim, cell
  real(dp) :: cellCenterCoord

  associate(cellCenter => mesh%CellCenterVec(cell))
    cellCenterCoord = cellCenter(dim)
  end associate

end function GetMeshCellCenterCoord

!! ----------------------------------------------------------------- !!
!! Check if cell connects it to the BC cell or 
!! to the interioir periodic. 
!! ----------------------------------------------------------------- !!
pure logical function IsMeshCellConnPeriodic(mesh, cellConn, cell) result(isPeriodicFace)
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: cell, cellConn
  
  isPeriodicFace = allocated(mesh%mIsCellConnPeriodic)
  if (isPeriodicFace) then
    isPeriodicFace = mesh%mIsCellConnPeriodic(cellConn, cell)
  end if

end function IsMeshCellConnPeriodic

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Initialize the mesh from image pixels.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine InitMeshFromImage(mesh, dl, stencil, image, &
    & fluidColor, colorToMark, numBndLayers, separatedBndCells)
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: dl(:)
  integer(ip), intent(in) :: image(0:,0:)
  integer(ip), intent(in) :: fluidColor, colorToMark(:)
  character(len=*), intent(in) :: stencil
  integer(ip), intent(in) :: numBndLayers
  logical, intent(in) :: separatedBndCells

  integer(ip) :: x, y, xx, yy, incAddr, boundaryLayer
  integer(ip) :: cell, cellCell, cellConn, bndCell, bndCellCell, mark
  integer(ip), allocatable :: cache(:,:)
  integer(ip), allocatable :: connIncr(:,:)

  !! TODO: Unified 2D/3D case.
  !! TODO: Check image for correct boundaries. 

  if (numBndLayers > 1.and.(.not.separatedBndCells)) then
    call ErrorStop('Non-separated boundary cells '// &
                 & 'support only a single boundary layer.')
  end if

  ! ----------------------
  ! Set mesh stencil.
  ! ----------------------
  select case(stencil)
    case('D1Q2')
      mesh%NumDims = 1
      mesh%NumCellFaces = 2
      mesh%NumCellConns = 2
      allocate(connIncr(1,2))
      connIncr(:,1) = [+1]; connIncr(:,2) = [-1]

    case('D1Q3')
      mesh%NumDims = 1
      mesh%NumCellFaces = 2
      mesh%NumCellConns = 3
      allocate(connIncr(1,3))
      connIncr(:,3) = [ 0]
      connIncr(:,1) = [+1]; connIncr(:,2) = [-1]

    case('D2Q4')
      mesh%NumDims = 2
      mesh%NumCellFaces = 4
      mesh%NumCellConns = 4
      allocate(connIncr(2,4))
      connIncr(:,1) = [+1, 0]; connIncr(:,2) = [-1, 0]
      connIncr(:,3) = [ 0,+1]; connIncr(:,4) = [ 0,-1]

    case('D2Q5')
      mesh%NumDims = 2
      mesh%NumCellFaces = 4
      mesh%NumCellConns = 5
      allocate(connIncr(2,5))
      connIncr(:,5) = [ 0, 0]
      connIncr(:,1) = [+1, 0]; connIncr(:,2) = [-1, 0]
      connIncr(:,3) = [ 0,+1]; connIncr(:,4) = [ 0,-1]

    case('D2Q9')
      mesh%NumDims = 2
      mesh%NumCellFaces = 4
      mesh%NumCellConns = 9
      allocate(connIncr(2,9))
      connIncr(:,9) = [ 0, 0]
      connIncr(:,1) = [+1, 0]; connIncr(:,2) = [-1, 0]
      connIncr(:,3) = [ 0,+1]; connIncr(:,4) = [ 0,-1]
      connIncr(:,5) = [+1,+1]; connIncr(:,6) = [-1,-1]
      connIncr(:,7) = [+1,-1]; connIncr(:,8) = [-1,+1]

    case('D3Q6')
      mesh%NumDims = 3
      mesh%NumCellFaces = 6
      mesh%NumCellConns = 6
      allocate(connIncr(3,6))
      connIncr(:,1) = [+1, 0, 0]; connIncr(:,2) = [-1, 0, 0]
      connIncr(:,3) = [ 0,+1, 0]; connIncr(:,4) = [ 0,-1, 0]
      connIncr(:,5) = [ 0, 0,+1]; connIncr(:,6) = [ 0,-1,-1]
      error stop 'not implemented'

    case('D3Q7')
      mesh%NumDims = 3
      mesh%NumCellFaces = 6
      mesh%NumCellConns = 7
      allocate(connIncr(3,7))
      connIncr(:,7) = [ 0, 0, 0]
      connIncr(:,1) = [+1, 0, 0]; connIncr(:,2) = [-1, 0, 0]
      connIncr(:,3) = [ 0,+1, 0]; connIncr(:,4) = [ 0,-1, 0]
      connIncr(:,5) = [ 0, 0,+1]; connIncr(:,6) = [ 0,-1,-1]
      error stop 'not implemented'

    case('D3Q19')
      mesh%NumDims = 3
      mesh%NumCellFaces = 6
      mesh%NumCellConns = 19
      allocate(connIncr(3,19))
      connIncr(:,19) = [ 0, 0, 0]
      connIncr(:, 1) = [+1, 0, 0]; connIncr(:, 2) = [-1, 0, 0]
      connIncr(:, 3) = [ 0,+1, 0]; connIncr(:, 4) = [ 0,-1, 0]
      connIncr(:, 5) = [ 0, 0,+1]; connIncr(:, 6) = [ 0,-1,-1]
      connIncr(:, 7) = [+1,+1, 0]; connIncr(:, 8) = [-1,-1, 0]
      connIncr(:, 9) = [+1,-1, 0]; connIncr(:,10) = [-1,+1, 0]
      connIncr(:,11) = [ 0,+1,+1]; connIncr(:,12) = [ 0,-1,-1]
      connIncr(:,13) = [ 0,+1,-1]; connIncr(:,14) = [ 0,-1,+1]
      connIncr(:,15) = [+1, 0,+1]; connIncr(:,16) = [-1, 0,-1]
      connIncr(:,17) = [+1, 0,-1]; connIncr(:,18) = [-1, 0,+1]
      error stop 'not implemented'

    case('D3Q27')
      mesh%NumDims = 3
      mesh%NumCellFaces = 6
      mesh%NumCellConns = 27
      allocate(connIncr(3,27))
      connIncr(:,27) = [ 0, 0, 0]
      connIncr(:, 1) = [+1, 0, 0]; connIncr(:, 2) = [-1, 0, 0]
      connIncr(:, 3) = [ 0,+1, 0]; connIncr(:, 4) = [ 0,-1, 0]
      connIncr(:, 5) = [ 0, 0,+1]; connIncr(:, 6) = [ 0,-1,-1]
      connIncr(:, 7) = [+1,+1, 0]; connIncr(:, 8) = [-1,-1, 0]
      connIncr(:, 9) = [+1,-1, 0]; connIncr(:,10) = [-1,+1, 0]
      connIncr(:,11) = [ 0,+1,+1]; connIncr(:,12) = [ 0,-1,-1]
      connIncr(:,13) = [ 0,+1,-1]; connIncr(:,14) = [ 0,-1,+1]
      connIncr(:,15) = [+1, 0,+1]; connIncr(:,16) = [-1, 0,-1]
      connIncr(:,17) = [+1, 0,-1]; connIncr(:,18) = [-1, 0,+1]
      connIncr(:,19) = [+1,+1,+1]; connIncr(:,20) = [-1,-1,-1]
      connIncr(:,21) = [+1,+1,-1]; connIncr(:,22) = [-1,-1,+1]
      connIncr(:,23) = [+1,-1,+1]; connIncr(:,24) = [-1,+1,-1]
      connIncr(:,25) = [+1,+1,-1]; connIncr(:,26) = [-1,-1,+1]
      error stop 'not implemented'

    case default
      call ErrorStop('Invalid stencil, expected one of: `D1Q2`, '// &
        & '`D1Q3`, `D2Q4`, `D2Q5`, `D2Q9`, `D3Q6`, `D3Q7`, `D3Q19`, `D3Q27`.')
  end select
  allocate(mesh%dl(mesh%NumCellConns))
  allocate(mesh%dr(mesh%NumDims,mesh%NumCellConns))
  do cellConn = 1, mesh%NumCellConns
    mesh%dr(:,cellConn) = 1.0_dp*connIncr(:,cellConn)
    mesh%dl(cellConn) = norm2(dl(:)*connIncr(:,cellConn))
  end do

  ! ----------------------
  ! Cache cell indices and count the cells.
  ! ----------------------
  mesh%NumCells = 0
  mesh%CellIndexBounds = shape(image) - 2
  allocate(cache(0:(mesh%CellIndexBounds(1) + 1), &
               & 0:(mesh%CellIndexBounds(2) + 1))) 
  cache(:,:) = 0
  mesh%NumBndMarks = size(colorToMark)
  allocate(mesh%BndCellAddrs(mesh%NumBndMarks + 1))
  mesh%BndCellAddrs(:) = 0
  do y = 1, mesh%CellIndexBounds(2)
    do x = 1, mesh%CellIndexBounds(1)
      if (image(x,y) /= fluidColor) cycle

      ! ----------------------
      ! Count and cache the interior cell.
      ! ----------------------
      mesh%NumCells = mesh%NumCells + 1
      cell = mesh%NumCells; cache(x,y) = cell

      ! ----------------------
      ! Go through the connected pixels 
      ! and count the BC cells per each mark.
      ! ----------------------
      do cellConn = 1, mesh%NumCellConns
        xx = x + connIncr(1,cellConn); yy = y + connIncr(2,cellConn)
        if ((image(xx,yy) /= fluidColor).and. &
            & (separatedBndCells.or.(cache(xx,yy) == 0))) then

          mark = IndexOf(image(xx,yy), colorToMark)
          if (mark == 0) then
            call ErrorStop('Unexpected image color at ('//I2S([xx,yy])//')')
          end if
          cache(xx,yy) = -mark
          mesh%BndCellAddrs(mark + 1) = mesh%BndCellAddrs(mark + 1) + 1

        end if
      end do
    end do
  end do

  ! ----------------------
  ! Precompute the boundary cells indices
  ! and count the total amount of the cells.
  ! ----------------------
  mesh%BndCellAddrs(1) = mesh%NumCells + 1
  do mark = 1, mesh%NumBndMarks

    if (mesh%BndCellAddrs(mark + 1) == 0) then
      call PrintWarning('Unused boundary mark '//I2S(mark)//' detected.')      
    end if
    mesh%BndCellAddrs(mark + 1) = &
      & mesh%BndCellAddrs(mark + 1) + mesh%BndCellAddrs(mark)

  end do
  mesh%NumBndCells = mesh%BndCellAddrs(mesh%NumBndMarks + 1) - mesh%NumCells - 1
  mesh%NumAllCells = mesh%NumCells + numBndLayers*mesh%NumBndCells

  ! ----------------------
  ! Fill the interiour cell indices and connectivity,
  ! generate the boundary cells, fill indices and partially connectivity.
  ! ----------------------
  allocate(mesh%CellIndex(mesh%NumDims,mesh%NumAllCells))
  allocate(mesh%CellToCell(mesh%NumCellConns,mesh%NumAllCells))
  allocate(mesh%BndCellConns( &
    & mesh%BndCellAddrs(1):(mesh%BndCellAddrs(mesh%NumBndMarks + 1) - 1)))
  mesh%NumAllCells = mesh%NumCells + mesh%NumBndCells
  do y = 1, mesh%CellIndexBounds(2)
    do x = 1, mesh%CellIndexBounds(1)
      if (image(x,y) /= fluidColor) cycle
      cell = cache(x,y)

      ! ----------------------
      ! Fill cell index.
      ! ----------------------
      mesh%CellIndex(:,cell) = [x,y]

      ! ----------------------
      ! Fill cell connectivity.
      ! ----------------------
      mesh%CellToCell(:,cell) = 0
      do cellConn = 1, mesh%NumCellConns
        xx = x + connIncr(1,cellConn); yy = y + connIncr(2,cellConn)
        if (cache(xx,yy) > 0) then

          ! ----------------------
          ! Connect to the existing cell.
          ! If the case of boundary cell, also connect to
          ! the current cell.
          ! ----------------------
          cellCell = cache(xx,yy)
          mesh%CellToCell(cellConn,cell) = cellCell
          if (cellCell > mesh%NumCells) then
            bndCell = cellCell
            mesh%CellToCell(Flip(cellConn),bndCell) = cell
            !mesh%BndCellConns(bndCell) = min(mesh%BndCellConns(bndCell), cellConn)
          end if

        else if (cache(xx,yy) < 0) then

          ! ----------------------
          ! Connect to the non-existing boundary cell.
          ! ----------------------
          mark = -cache(xx,yy)
          bndCell = mesh%BndCellAddrs(mark)
          mesh%BndCellConns(mesh%BndCellAddrs(mark)) = cellConn
          mesh%BndCellAddrs(mark) = mesh%BndCellAddrs(mark) + 1

          ! ----------------------
          ! Cache the cell cell to avoid duplicates.
          ! ----------------------
          if (.not.separatedBndCells) cache(xx,yy) = bndCell

          ! ----------------------
          ! Fill the boundary cell index
          ! and connection to the interior cell.
          ! ----------------------
          mesh%CellIndex(:,bndCell) = [xx,yy]
          mesh%CellToCell(cellConn,cell) = bndCell
          mesh%CellToCell(Flip(cellConn),bndCell) = cell

          ! ----------------------
          ! Generate boundary layers.
          ! ----------------------
          do boundaryLayer = 2, numBndLayers
            mesh%NumAllCells = mesh%NumAllCells + 1
            bndCellCell = mesh%NumAllCells

            mesh%CellToCell(cellConn,bndCell) = bndCellCell
            mesh%CellToCell(Flip(cellConn),bndCellCell) = bndCell
            xx = xx + connIncr(1,cellConn); yy = yy + connIncr(2,cellConn)
            mesh%CellIndex(:,bndCellCell) = [xx,yy]

            bndCell = bndCellCell
          end do

        else

          call ErrorStop('Invalid value `0` in cache at ('//I2S([xx,yy])//').')

        end if
      end do
    end do
  end do

  ! ----------------------
  ! Fix the boundary mark addresses.
  ! ----------------------
  mesh%BndCellAddrs = cshift(mesh%BndCellAddrs, -1)
  mesh%BndCellAddrs(1) = mesh%NumCells + 1

  ! ----------------------
  ! Print statistics.
  ! ----------------------
  print *
  print *, '-=-=-=-=-=-=-=-'
  print *, 'Mesh statistics (InitFromImage):'
  print *, '-=-=-=-=-=-=-=-'
  print *, ' * Stencil:                  '//stencil
  print *, ' * Index bounds:             '//I2S(mesh%CellIndexBounds)
  print *, ' * Number of cells:          '//I2S(mesh%NumCells)
  print *, ' * Number of boundary cells: '//I2S(mesh%NumBndCells)
  print *, ' * Number of all cells:      '//I2S(mesh%NumAllCells)
  print *

end subroutine InitMeshFromImage

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Apply the specified cell ordering.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine ApplyMeshOrdering(mesh, iperm)
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(in) :: iperm(:)

  integer(ip) :: cell, cellConn
  integer(ip), allocatable :: cellToCellTemp(:,:)
  integer(ip), allocatable :: cellIndexTemp(:,:)

  ! ----------------------
  ! Order cell-cell graph columns.
  ! ----------------------
  do cell = 1, mesh%NumAllCells
    do cellConn = 1, mesh%NumCellFaces
      associate(cellCell => mesh%CellToCell(cellConn, cell))

        if (0 < cellCell.and.cellCell <= mesh%NumCells) cellCell = iperm(cellCell)

      end associate
    end do
  end do

  ! ----------------------
  ! Order cell-cell graph rows.
  ! ----------------------
  cellToCellTemp = mesh%CellToCell(:,:mesh%NumCells)
  cellIndexTemp = mesh%CellIndex(:,:mesh%NumCells)
  do cell = 1, mesh%NumCells
    associate(cellCell => iperm(cell))

      mesh%CellToCell(:,cellCell) = cellToCellTemp(:,cell)
      mesh%CellIndex(:,cellCell) = cellIndexTemp(:,cell)

    end associate
  end do

end subroutine ApplyMeshOrdering

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Print mesh connectivity in '.dot' format.
!! Output may be visualized with: 'neato -n -Tpng c2c.dot > c2c.png' 
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine PrintMeshToNeato(mesh, file)
  class(tMesh), intent(inout) :: mesh
  character(len=*), intent(in) :: file

  integer(ip) :: unit
  integer(ip) :: cell, cellCell, cellConn
  integer(ip) :: mark, bndCell, bndCellConn, gstCell, gstCellCell

  integer(ip), parameter :: DPI = 72

  character(len=10), parameter :: SHAPES(*) = &
    & [ character(len=10) :: 'box', 'polygon', 'egg', 'triangle', &
    &   'diamond', 'trapezium', 'parallelogram', 'house', 'pentagon', &
    &   'haxagon', 'septagon', 'octagon', 'star' ]
  
  character(len=10), parameter :: PALETTE(*) = &
    & [ character(len=10) :: 'red', 'green', 'blue', 'fuchsia', &
    &   'yellow', 'wheat' ]

  print *
  print *, '-=-=-=-=-=-=-=-'
  print *, 'Print to Neato: '//"'"//file//"'"
  print *, '-=-=-=-=-=-=-=-'
  print *
  
  if (mesh%NumDims /= 2) then
    call ErrorStop('Only 2D meshes can be printed to Neato')
  end if

  open(newunit=unit, file=file, status='replace')
  
  write(unit, "(A)") '# Generated by StormRuler/Mesh2Neato'
  write(unit, "(A)") 'digraph Mesh {'
  
  ! ----------------------
  ! Write cell MD indices of the cells as graph nodes.
  ! ----------------------
  ! Difference node shapes are use for different face
  ! directions of the BC/ghost cells faces.
  do cell = 1, mesh%NumCells
    associate(cellCoord => DPI*mesh%CellIndex(:,cell))

      write(unit, "('  ', A, '[label=', A, ' ', 'pos=', A, ']')") &
        & 'C'//I2S(cell), '"C"', & ! ID and label
        & '"'//I2S(cellCoord(1))//','//I2S(cellCoord(2))//'!"' ! pos

    end associate
  end do

  do mark = 1, mesh%NumBndMarks
    do bndCell = mesh%BndCellAddrs(mark), mesh%BndCellAddrs(mark + 1) - 1
      bndCellConn = mesh%BndCellConns(bndCell)
      gstCell = bndCell
      do while(gstCell /= 0)
        associate(cellCoord => DPI*mesh%CellIndex(:,gstCell))

          write(unit, "('  ', A, '[label=', A, ' pos=', A, ' shape=', A, ']')") &
            & 'C'//I2S(gstCell), & ! ID
            & '"B'//I2S(mark)//'"', & ! label
            & '"'//I2S(cellCoord(1))//','//I2S(cellCoord(2))//'!"', & ! pos
            & trim(SHAPES((bndCellConn+1)/2)) ! shape

        end associate
        gstCell = mesh%CellToCell(bndCellConn, gstCell)
      end do
    end do
  end do
  
  ! ----------------------
  ! Write cell to cell connectivity as graph edges.
  ! ----------------------
  ! BC/ghost cells are colored by the BC mark.
  do cell = 1, mesh%NumCells
    do cellConn = 1, mesh%NumCellConns
      ! Do not dump the periodic connections.
      if (mesh%IsCellConnPeriodic(cellConn, cell)) cycle

      cellCell = mesh%CellToCell(cellConn, cell)
      write(unit, "('  ', A, '->', A)") &
        & 'C'//I2S(cell), 'C'//I2S(cellCell) ! source and dest ID

    end do
  end do

  do mark = 1, mesh%NumBndMarks
    do bndCell = mesh%BndCellAddrs(mark), mesh%BndCellAddrs(mark + 1) - 1
      bndCellConn = mesh%BndCellConns(bndCell)
      gstCell = bndCell
      do
        associate(color => PALETTE(mark))

          gstCellCell = mesh%CellToCell(Flip(bndCellConn), gstCell) 
          write(unit, "('  ', A, '->', A, ' [color=', A, ']')") &
            & 'C'//I2S(gstCell), & ! source ID
            & 'C'//I2S(gstCellCell), color ! dest ID and color

          gstCellCell = mesh%CellToCell(bndCellConn, gstCell) 
          if (gstCellCell == 0) exit
          write(unit, "('  ', A, '->', A, ' [color=', A, ']')") &
            & 'C'//I2S(gstCell), & ! source ID
            & 'C'//I2S(gstCellCell), color ! dest ID and color

        end associate
        gstCell = gstCellCell
      end do
    end do
  end do

  write(unit, "(A)") '}'
  
  ! ----------------------
  ! Close file and exit.
  ! ----------------------
  close(unit)

end subroutine PrintMeshToNeato

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Print mesh in Legacy VTK '.vtk' format.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine PrintMeshToLegacyVTK(mesh, file, fields)
  class(tMesh), intent(inout) :: mesh
  character(len=*), intent(in) :: file
  type(IOList), intent(in), optional :: fields
  
  integer(ip) :: unit
  integer(ip) :: numVtkCells, vtkCell
  integer(ip) :: cell, cellCell, cellConn, cellCellCell, cellCellFace
  logical, allocatable :: cellIsFullInternal(:)
  class(IOListItem), pointer :: item

  !! TODO: it feels like this thing can be seriously optimized..

  if ((mesh%NumDims /= 2).and.(mesh%NumDims /= 3)) then
    error stop 'Only 2D/3D meshes can be printed to Legacy VTK'
  end if

  print *, 'Print to Legacy VTK: ', file
  
  open(newunit=unit, file=file, status='replace')

  write(unit, "(A)") '# vtk DataFile Version 2.0'
  write(unit, "(A)") '# Generated by StormRuler/Mesh2VTK'
  write(unit, "(A)") 'ASCII'
  write(unit, "(A)") 'DATASET UNSTRUCTURED_GRID'
  write(unit, "(A)") ''
  
  ! ----------------------
  ! Write cell centers as VTK nodes.
  ! ----------------------
  write(unit, "('POINTS ', A, ' double')") I2S(mesh%NumCells)
  do cell = 1, mesh%NumCells
    associate(r => mesh%CellCenter(cell))
      if (mesh%NumDims == 2) then
        write(unit, "(A, ' ', A, ' ', A)") R2S(r(1)), R2S(r(2)), '0.0'
      else
        write(unit, "(A, ' ', A, ' ', A)") R2S(r(1)), R2S(r(2)), R2S(r(3))
      end if
    end associate
  end do
  write(unit, "(A)") ''
  
  ! ----------------------
  ! Precompute if cell is far away from boundaries.
  ! ----------------------
  allocate(cellIsFullInternal(mesh%NumCells))
  !$omp parallel do
  do cell = 1, mesh%NumCells
    cellIsFullInternal(cell) = &
      & all(mesh%CellToCell(:,cell) <= mesh%NumCells)
  end do
  !$omp end parallel do
  
  ! ----------------------
  ! Write VTK cells.
  ! VTK cell connectivity is generated on the fly via parity algorithm.
  ! ----------------------
#$do writePass = 1, 2

  numVtkCells = 0
  !$omp parallel do reduction(+:numVtkCells) &
  !$omp private(cellConn, cellCell, cellCellFace, cellCellCell) if($writePass == 1)
  do cell = 1, mesh%NumCells

    ! ----------------------
    ! Locate the second cell (VTK cell node).
    ! ----------------------
    do cellConn = 1, mesh%NumCellFaces
      if (mesh%IsCellConnPeriodic(cellConn, cell)) cycle
      cellCell = mesh%CellToCell(cellConn, cell)
      if (cellCell > mesh%NumCells) cycle

      ! ----------------------
      ! Locate the third cell (VTK cell node).
      ! ----------------------
      do cellCellFace = 1, mesh%NumCellFaces
        if (cellCellFace == cellConn) cycle
        if (mesh%IsCellConnPeriodic(cellCellFace, cell)) cycle
        cellCellCell = mesh%CellToCell(cellCellFace, cell)
        if (cellCellCell > mesh%NumCells) cycle

        ! ----------------------
        ! Generate the VTK cell (triangle).
        ! ----------------------
        ! Denote cells with same index components parity as primary.
        ! Primary cells are used to generate VTK cells around them.
        ! Secondary cells are used to generate VTK cells near boundaries.
        associate(cellIndex => mesh%CellIndex(:,cell))
          if (mod(cellIndex(1) - cellIndex(2), 2) /= 0) then
            ! Skip secondary cells far away from boundaries.
            if ( cellIsFullInternal(cellCell).or. &
              &  cellIsFullInternal(cellCellCell) ) cycle
          end if
        end associate

        numVtkCells = numVtkCells + 1
#$if writePass == 2
        write(unit, "('3 ', A, ' ', A, ' ', A)") &
          & I2S(cell-1), I2S(cellCell-1), I2S(cellCellCell-1)
#$end if

      end do
    end do
  end do
  !$omp end parallel do

  ! ----------------------
  ! Write a number of the VTK cells as the result of the first pass.
  ! ----------------------
#$if writePass == 1
  associate(numVTKCellNodes => numVtkCells*(mesh%NumDims + 2))
    write(unit, "('CELLS ', A, ' ', A)") I2S(numVtkCells), I2S(numVTKCellNodes)
  end associate
#$end if

#$end do

  ! ----------------------
  ! Write VTK cell types.
  ! ----------------------
  write(unit, "('CELL_TYPES ', A)") I2S(numVtkCells)
  do vtkCell = 1, numVtkCells
    write(unit, "(A)") '5' ! 5 is triangle.
  end do
  write(unit, "(A)") ''

  ! ----------------------
  ! Write fields.
  ! ----------------------
  write(unit, "('POINT_DATA ', A)") I2S(mesh%NumCells)
  if (present(fields)) then
    item => fields%first
    do while(associated(item))
      select type(item)
        ! ----------------------
        ! Scalar field.
        ! ----------------------
        class is(IOListItem$0)
          write(unit, "('SCALARS ', A, ' double 1')") item%name 
          write(unit, "(A)") 'LOOKUP_TABLE default'
          do cell = 1, mesh%NumCells
            write(unit, "(A)") R2S(item%values(cell))
          end do

        ! ----------------------
        ! Vector field.
        ! ----------------------
        class is(IOListItem$1)
          write(unit, "('VECTORS ', A, ' double')") item%name 
          do cell = 1, mesh%NumCells
            associate(v => item%values(:,cell))
              write(unit, '(A, " ", A, " ", A)') R2S(v(1)), R2S(v(2)), '0.0'
            end associate
          end do

        end select
      write(unit, "(A)") ''
      item => item%next
    end do
  end if

  ! ----------------------
  ! Close file and exit.
  ! ----------------------
  close(unit)

end subroutine PrintMeshToLegacyVTK

end module StormRuler_Mesh
