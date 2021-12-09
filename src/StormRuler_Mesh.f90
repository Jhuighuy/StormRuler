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

use StormRuler_Consts, only: bp, ip, dp

use StormRuler_Helpers, only: &
  & ErrorStop, PrintLog, PrintWarning, &
  & Flip, IndexOf, InsertTo, I2S, R2S, IntToRgb

use, intrinsic :: iso_fortran_env, only: error_unit

#$use 'StormRuler_Macros.fi'

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
  ! ----------------------
  character(len=:), allocatable :: Stencil 
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
  integer(ip), allocatable :: IndexBounds(:)
  ! ----------------------
  ! Multidimensional index to cell table.
  ! Shape is [1, IndexBounds(1)]×…×[1, IndexBounds(NumDims)].
  !! TODO: 3D case.
  ! ----------------------
  integer(ip), allocatable :: IndexToCell(:,:)
  ! ----------------------
  ! Cell multidimensional index table.
  ! Shape is [1, NumDims]×[1, NumAllCells].
  ! ----------------------
  integer(ip), allocatable :: CellToIndex(:,:)

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
  integer(bp), allocatable :: BndCellConns(:)

  ! ----------------------
  ! Distance between centers of the adjacent cells per face.
  ! Shape is [1, NumCellConns].
  ! ----------------------
  real(dp), allocatable :: dl(:)
  ! ----------------------
  ! Unit difference between centers of the adjacent cells per face.
  ! Shape is [1, NumDims]×[1, NumCellConns].
  ! ----------------------
  integer(ip), allocatable :: dr(:,:)

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
  procedure :: IsCellRed => IsMeshCellRed
  procedure :: IsCellInternal => IsMeshCellInternal
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
!! Check if mesh cell is red or black.
!! ----------------------------------------------------------------- !!
elemental logical function IsMeshCellRed(mesh, cell) result(isRed)
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: cell

  associate(cellIndex => mesh%CellToIndex(:,cell))
    isRed = mod(cellIndex(1) - cellIndex(2), 2) == 0!all(mod(cellIndex, 2) == 1)
  end associate

end function IsMeshCellRed

!! ----------------------------------------------------------------- !!
!! Check if mesh cell is internal or ghost.
!! ----------------------------------------------------------------- !!
elemental logical function IsMeshCellInternal(mesh, cell) result(isInternal)
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: cell

  isInternal = (1 <= cell).and.(cell <= mesh%NumCells)

end function IsMeshCellInternal

!! ----------------------------------------------------------------- !!
!! Get cell center.
!! ----------------------------------------------------------------- !!
pure function GetMeshCellCenterVec(mesh, cell) result(cellCenter)
  class(tMesh), intent(in) :: mesh
  integer(ip), intent(in) :: cell
  real(dp) :: cellCenter(mesh%NumDims)

  !! TODO: add initial offset as mesh class field.
  associate(cellIndex => mesh%CellToIndex(:,cell))
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
!! Initialize the mesh stencil.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine InitMeshStencil(mesh, spacings, stencil)
  class(tMesh), intent(inout) :: mesh
  real(dp), intent(in) :: spacings(:)
  character(len=*), intent(in) :: stencil

  integer(ip) :: cellConn

  mesh%Stencil = stencil
  select case(mesh%Stencil)

    case('D1Q2')
      mesh%NumDims = 1
      mesh%NumCellFaces = 2
      mesh%NumCellConns = 2
      allocate(mesh%dr(1,2))
      mesh%dr(:,1) = [+1]; mesh%dr(:,2) = [-1]

    case('D1Q3')
      mesh%NumDims = 1
      mesh%NumCellFaces = 2
      mesh%NumCellConns = 3
      allocate(mesh%dr(1,3))
      mesh%dr(:,3) = [ 0]
      mesh%dr(:,1) = [+1]; mesh%dr(:,2) = [-1]

    case('D2Q4')
      mesh%NumDims = 2
      mesh%NumCellFaces = 4
      mesh%NumCellConns = 4
      allocate(mesh%dr(2,4))
      mesh%dr(:,1) = [+1, 0]; mesh%dr(:,2) = [-1, 0]
      mesh%dr(:,3) = [ 0,+1]; mesh%dr(:,4) = [ 0,-1]

    case('D2Q5')
      mesh%NumDims = 2
      mesh%NumCellFaces = 4
      mesh%NumCellConns = 5
      allocate(mesh%dr(2,5))
      mesh%dr(:,5) = [ 0, 0]
      mesh%dr(:,1) = [+1, 0]; mesh%dr(:,2) = [-1, 0]
      mesh%dr(:,3) = [ 0,+1]; mesh%dr(:,4) = [ 0,-1]

    case('D2Q9')
      mesh%NumDims = 2
      mesh%NumCellFaces = 4
      mesh%NumCellConns = 9
      allocate(mesh%dr(2,9))
      mesh%dr(:,9) = [ 0, 0]
      mesh%dr(:,1) = [+1, 0]; mesh%dr(:,2) = [-1, 0]
      mesh%dr(:,3) = [ 0,+1]; mesh%dr(:,4) = [ 0,-1]
      mesh%dr(:,5) = [+1,+1]; mesh%dr(:,6) = [-1,-1]
      mesh%dr(:,7) = [+1,-1]; mesh%dr(:,8) = [-1,+1]

    case('D3Q6')
      mesh%NumDims = 3
      mesh%NumCellFaces = 6
      mesh%NumCellConns = 6
      allocate(mesh%dr(3,6))
      mesh%dr(:,1) = [+1, 0, 0]; mesh%dr(:,2) = [-1, 0, 0]
      mesh%dr(:,3) = [ 0,+1, 0]; mesh%dr(:,4) = [ 0,-1, 0]
      mesh%dr(:,5) = [ 0, 0,+1]; mesh%dr(:,6) = [ 0,-1,-1]

    case('D3Q7')
      mesh%NumDims = 3
      mesh%NumCellFaces = 6
      mesh%NumCellConns = 7
      allocate(mesh%dr(3,7))
      mesh%dr(:,7) = [ 0, 0, 0]
      mesh%dr(:,1) = [+1, 0, 0]; mesh%dr(:,2) = [-1, 0, 0]
      mesh%dr(:,3) = [ 0,+1, 0]; mesh%dr(:,4) = [ 0,-1, 0]
      mesh%dr(:,5) = [ 0, 0,+1]; mesh%dr(:,6) = [ 0,-1,-1]

    case('D3Q19')
      mesh%NumDims = 3
      mesh%NumCellFaces = 6
      mesh%NumCellConns = 19
      allocate(mesh%dr(3,19))
      mesh%dr(:,19) = [ 0, 0, 0]
      mesh%dr(:, 1) = [+1, 0, 0]; mesh%dr(:, 2) = [-1, 0, 0]
      mesh%dr(:, 3) = [ 0,+1, 0]; mesh%dr(:, 4) = [ 0,-1, 0]
      mesh%dr(:, 5) = [ 0, 0,+1]; mesh%dr(:, 6) = [ 0,-1,-1]
      mesh%dr(:, 7) = [+1,+1, 0]; mesh%dr(:, 8) = [-1,-1, 0]
      mesh%dr(:, 9) = [+1,-1, 0]; mesh%dr(:,10) = [-1,+1, 0]
      mesh%dr(:,11) = [ 0,+1,+1]; mesh%dr(:,12) = [ 0,-1,-1]
      mesh%dr(:,13) = [ 0,+1,-1]; mesh%dr(:,14) = [ 0,-1,+1]
      mesh%dr(:,15) = [+1, 0,+1]; mesh%dr(:,16) = [-1, 0,-1]
      mesh%dr(:,17) = [+1, 0,-1]; mesh%dr(:,18) = [-1, 0,+1]

    case('D3Q27')
      mesh%NumDims = 3
      mesh%NumCellFaces = 6
      mesh%NumCellConns = 27
      allocate(mesh%dr(3,27))
      mesh%dr(:,27) = [ 0, 0, 0]
      mesh%dr(:, 1) = [+1, 0, 0]; mesh%dr(:, 2) = [-1, 0, 0]
      mesh%dr(:, 3) = [ 0,+1, 0]; mesh%dr(:, 4) = [ 0,-1, 0]
      mesh%dr(:, 5) = [ 0, 0,+1]; mesh%dr(:, 6) = [ 0,-1,-1]
      mesh%dr(:, 7) = [+1,+1, 0]; mesh%dr(:, 8) = [-1,-1, 0]
      mesh%dr(:, 9) = [+1,-1, 0]; mesh%dr(:,10) = [-1,+1, 0]
      mesh%dr(:,11) = [ 0,+1,+1]; mesh%dr(:,12) = [ 0,-1,-1]
      mesh%dr(:,13) = [ 0,+1,-1]; mesh%dr(:,14) = [ 0,-1,+1]
      mesh%dr(:,15) = [+1, 0,+1]; mesh%dr(:,16) = [-1, 0,-1]
      mesh%dr(:,17) = [+1, 0,-1]; mesh%dr(:,18) = [-1, 0,+1]
      mesh%dr(:,19) = [+1,+1,+1]; mesh%dr(:,20) = [-1,-1,-1]
      mesh%dr(:,21) = [+1,+1,-1]; mesh%dr(:,22) = [-1,-1,+1]
      mesh%dr(:,23) = [+1,-1,+1]; mesh%dr(:,24) = [-1,+1,-1]
      mesh%dr(:,25) = [+1,+1,-1]; mesh%dr(:,26) = [-1,-1,+1]

    case default
      call ErrorStop('Invalid stencil, expected one of: `D1Q2`, '// &
        & '`D1Q3`, `D2Q4`, `D2Q5`, `D2Q9`, `D3Q6`, `D3Q7`, `D3Q19`, `D3Q27`.')

  end select

  allocate(mesh%dl(mesh%NumCellConns))
  do cellConn = 1, mesh%NumCellConns
    mesh%dl(cellConn) = norm2(spacings(:)*mesh%dr(:,cellConn))
  end do

end subroutine InitMeshStencil

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Initialize the mesh from image pixels.
!! ( this function is more suited for FDM computations. )
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine InitMeshFromImage(mesh, image, fluidColor, colorToMark, &
                           & numBndLayers, separatedBndCells)
  class(tMesh), intent(inout) :: mesh
  integer(ip), intent(in) :: image(0:,0:)
  integer(ip), intent(in) :: fluidColor, colorToMark(:)
  integer(ip), intent(in) :: numBndLayers
  logical, intent(in) :: separatedBndCells

  integer(ip) :: x, y, xx, yy, incAddr, bndLayer
  integer(ip) :: cell, cellCell, cellConn, bndCell, bndCellCell, mark

  !! TODO: Unified 2D/3D case.
  !! TODO: Check image for correct boundaries.

  if (numBndLayers > 1.and.(.not.separatedBndCells)) then
    call ErrorStop('Non-separated boundary cells '// &
                 & 'support only a single boundary layer.')
  end if

  ! ----------------------
  ! Assign interior cell indices, count the boundary cells per mark.
  ! ----------------------
  mesh%NumCells = 0
  mesh%IndexBounds = shape(image) - 2
  allocate(mesh%IndexToCell(0:(mesh%IndexBounds(1) + 1), &
                          & 0:(mesh%IndexBounds(2) + 1))) 
  mesh%IndexToCell(:,:) = 0
  mesh%NumBndMarks = size(colorToMark)
  allocate(mesh%BndCellAddrs(mesh%NumBndMarks + 1))
  mesh%BndCellAddrs(:) = 0
  do y = 1, mesh%IndexBounds(2)
    do x = 1, mesh%IndexBounds(1)
      if (image(x,y) /= fluidColor) cycle

      ! ----------------------
      ! Assign index to the interior cell.
      ! ----------------------
      mesh%NumCells = mesh%NumCells + 1
      cell = mesh%NumCells; mesh%IndexToCell(x,y) = cell

      ! ----------------------
      ! Count the boundary cells.
      ! ----------------------
      do cellConn = 1, mesh%NumCellConns
        xx = x + mesh%dr(1,cellConn); yy = y + mesh%dr(2,cellConn)
        if ((image(xx,yy) /= fluidColor).and. &
            & (separatedBndCells.or.(mesh%IndexToCell(xx,yy) == 0))) then

          mark = IndexOf(image(xx,yy), colorToMark)
          if (mark == 0) then
            call ErrorStop('Unexpected image color at ('//I2S([xx,yy])//')')
          end if

          mesh%IndexToCell(xx,yy) = -mark
          mesh%BndCellAddrs(mark + 1) = mesh%BndCellAddrs(mark + 1) + 1

        end if
      end do
    end do
  end do

  ! ----------------------
  ! Assign the boundary cells indices and count the total amount of the cells.
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
  ! Fill the cell indices and connectivity.
  ! ----------------------
  allocate(mesh%CellToIndex(mesh%NumDims,mesh%NumAllCells))
  allocate(mesh%CellToCell(mesh%NumCellConns,mesh%NumAllCells))
  associate(lbound => mesh%BndCellAddrs(1), & 
          & rbound => mesh%BndCellAddrs(mesh%NumBndMarks + 1) - 1)
    allocate(mesh%BndCellConns(lbound:rbound))
  end associate
  mesh%NumAllCells = mesh%NumCells + mesh%NumBndCells
  do y = 1, mesh%IndexBounds(2)
    do x = 1, mesh%IndexBounds(1)
      if (image(x,y) /= fluidColor) cycle

      ! ----------------------
      ! Fill the cell index.
      ! ----------------------
      cell = mesh%IndexToCell(x,y)
      mesh%CellToIndex(:,cell) = [x,y]

      ! ----------------------
      ! Fill the cell connectivity.
      ! ----------------------
      mesh%CellToCell(:,cell) = 0
      do cellConn = 1, mesh%NumCellConns
        xx = x + mesh%dr(1,cellConn); yy = y + mesh%dr(2,cellConn)
        if (mesh%IndexToCell(xx,yy) > 0) then

          ! ----------------------
          ! Neighbour cell has the assigned index.
          ! ----------------------
          cellCell = mesh%IndexToCell(xx,yy)
          mesh%CellToCell(cellConn,cell) = cellCell

          ! ----------------------
          ! If the current cell is interior and the neighbour cell is boundary, 
          ! also connect it to the current one.
          ! ----------------------
          if ((cell <= mesh%NumCells).and.(cellCell > mesh%NumCells)) then
            bndCell = cellCell
            mesh%CellToCell(Flip(cellConn),bndCell) = cell
          end if

        else if (mesh%IndexToCell(xx,yy) < 0) then

          ! ----------------------
          ! Neighbour cell is boundary cell without the assigned index.
          ! ----------------------
          mark = -mesh%IndexToCell(xx,yy)
          bndCell = mesh%BndCellAddrs(mark)
          mesh%BndCellConns(bndCell) = cellConn
          mesh%BndCellAddrs(mark) = mesh%BndCellAddrs(mark) + 1

          ! ----------------------
          ! Fill the boundary cell index and connection to the interior cell.
          ! ----------------------
          mesh%CellToIndex(:,bndCell) = [xx,yy]
          mesh%CellToCell(cellConn,cell) = bndCell
          mesh%CellToCell(Flip(cellConn),bndCell) = cell

          ! ----------------------
          ! Assign the cell index to avoid duplicates in non-separated case.
          ! ----------------------
          if (.not.separatedBndCells) mesh%IndexToCell(xx,yy) = bndCell

          ! ----------------------
          ! Optionally generate boundary layers.
          ! ----------------------
          do bndLayer = 2, numBndLayers
            mesh%NumAllCells = mesh%NumAllCells + 1
            bndCellCell = mesh%NumAllCells
            mesh%CellToCell(cellConn,bndCell) = bndCellCell
            mesh%CellToCell(Flip(cellConn),bndCellCell) = bndCell
            xx = xx + mesh%dr(1,cellConn); yy = yy + mesh%dr(2,cellConn)
            mesh%CellToIndex(:,bndCellCell) = [xx,yy]
            bndCell = bndCellCell
          end do

        else

          call ErrorStop('Invalid value `0` in mesh%IndexToCell at ('//I2S([xx,yy])//').')

        end if
      end do
    end do
  end do

  ! ----------------------
  ! Fix the cells addresses.
  ! ----------------------
  mesh%BndCellAddrs = cshift(mesh%BndCellAddrs, -1)
  mesh%BndCellAddrs(1) = mesh%NumCells + 1

  ! ----------------------
  ! Print mesh statistics.
  ! ----------------------
  call PrintLog('')
  call PrintLog('-=-=-=-=-=-=-=-')
  call PrintLog('Mesh statistics (InitFromImage):')
  call PrintLog('-=-=-=-=-=-=-=-')
  call PrintLog(' * Stencil:                  '//mesh%Stencil)
  call PrintLog(' * Index bounds:             '//I2S(mesh%IndexBounds))
  call PrintLog(' * Number of cells:          '//I2S(mesh%NumCells))
  call PrintLog(' * Number of boundary cells: '//I2S(mesh%NumBndCells))
  call PrintLog(' * Number of all cells:      '//I2S(mesh%NumAllCells))
  call PrintLog('')

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
      associate(cellCell => mesh%CellToCell(cellConn,cell))

        if (0 < cellCell.and.cellCell <= mesh%NumCells) cellCell = iperm(cellCell)

      end associate
    end do
  end do

  ! ----------------------
  ! Order cell-cell graph rows.
  ! ----------------------
  cellToCellTemp = mesh%CellToCell(:,:mesh%NumCells)
  cellIndexTemp = mesh%CellToIndex(:,:mesh%NumCells)
  do cell = 1, mesh%NumCells
    associate(cellCell => iperm(cell))

      mesh%CellToCell(:,cellCell) = cellToCellTemp(:,cell)
      mesh%CellToIndex(:,cellCell) = cellIndexTemp(:,cell)

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

  character(len=10), parameter :: PALETTE(*) = &
    & [ character(len=10) :: 'red', 'green', 'blue', 'fuchsia', &
    &   'yellow', 'wheat' ]

  character(len=20), parameter :: SHAPES(*) = &
    & [ character(len=20) :: 'box', 'polygon', 'egg', 'triangle', &
    &   'diamond', 'trapezium', 'parallelogram', 'house', 'pentagon', &
    &   'haxagon', 'septagon', 'octagon', 'star' ]

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
    associate(cellCoord => DPI*mesh%CellToIndex(:,cell))

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
        associate(cellCoord => DPI*mesh%CellToIndex(:,gstCell))

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

end module StormRuler_Mesh
