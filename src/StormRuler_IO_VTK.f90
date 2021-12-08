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
module StormRuler_IO_VTK

use StormRuler_Consts, only: bp, ip, dp, endl

use StormRuler_Helpers, only: &
  & ErrorStop, PrintLog, PrintWarning, Flip, I2S, R2S

use StormRuler_Mesh, only: tMesh

use StormRuler_IO, only: IOList, IOListItem, @{IOListItem$$@|@0, 2}@
use StormRuler_IO_Encoder, only: tEncoder, tAsciiEncoder, tBinaryEncoder
!use StormRuler_IO_Base64, only: tBase64OutputStream

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! ...
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
type :: tVtkMesh
  ! ----------------------
  ! Type of the VTK cell.
  ! ----------------------
  integer(bp) :: VtkCellType

  ! ----------------------
  ! Number of the VTK cells.
  ! ----------------------
  integer(ip) :: NumVtkCells
  ! ----------------------
  ! Number of nodes per VTK cell.
  ! ----------------------
  integer(ip) :: NumVtkCellNodes

  ! ----------------------
  ! VTK cell-VTK node connectivity table. 
  ! Shape is [1, NumVtkCellNodes]×[1, NumVtkCells].
  ! ----------------------
  integer(ip), allocatable :: VtkCellToNode(:,:)
end type tVtkMesh

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! ----------------------------------------------------------------- !!
!! Initialize VTK mesh.
!! ----------------------------------------------------------------- !!
subroutine IO_InitVtkMesh(vtkMesh, mesh)
  class(tVtkMesh), intent(inout) :: vtkMesh
  class(tMesh), intent(in) :: mesh

  integer(ip) :: cell, cellCells(4), cellFace, orthCellFaces(2)

  vtkMesh%VtkCellType = merge(5_bp, 10_bp, mesh%NumDims == 2)
  vtkMesh%NumVtkCells = 0
  vtkMesh%NumVtkCellNodes = mesh%NumDims + 1
  allocate(vtkMesh%VtkCellToNode(vtkMesh%NumVtkCellNodes, 4*mesh%NumCells))

  ! ----------------------
  ! Generate the VTK mesh connectivity.
  ! ----------------------
  do cell = 1, mesh%NumCells
    if (mesh%IsCellRed(cell)) cycle

    ! ----------------------
    ! Locate the first internal neighbour.
    ! ----------------------
    do cellFace = 1, mesh%NumCellFaces
      cellCells(1) = mesh%CellToCell(cellFace,cell)
      if (mesh%IsCellInternal(cellCells(1))) exit
    end do

    if (mesh%NumDims == 2) then

      orthCellFaces(1) = 2*mod((cellFace - 1)/2 + 1, 2) + 1
      cellCells(2) = mesh%CellToCell(orthCellFaces(1),cell)
      cellCells(3) = mesh%CellToCell(Flip(orthCellFaces(1)),cell)
      if (any(mesh%IsCellInternal(cellCells(2:3)))) then

        ! ----------------------
        ! This is a regular cell, 
        ! generate the VTK triangles around it.
        ! ----------------------
        call MakeVtkCell(cell, cellCells)
        cellCells(1) = mesh%CellToCell(Flip(cellFace),cell)
        if (mesh%IsCellInternal(cellCells(1))) then
          call MakeVtkCell(cell, cellCells)
        end if

      else

        ! ----------------------
        ! This is a hanging cell, 
        ! generate the VTK triangles around the neighbour.
        ! ----------------------
        cellCells(2) = mesh%CellToCell(orthCellFaces(1),cellCells(1))
        cellCells(3) = mesh%CellToCell(Flip(orthCellFaces(1)),cellCells(1))
        call MakeVtkCell(cell, cellCells)

      end if

    else if (mesh%NumDims == 2) then

      orthCellFaces(1) = 2*mod((cellFace - 1)/2 + 1, 3) + 1
      orthCellFaces(2) = 2*mod((cellFace - 1)/2 + 2, 3) + 1
      cellCells(2) = mesh%CellToCell(orthCellFaces(1),cell)
      cellCells(3) = mesh%CellToCell(orthCellFaces(2),cell)
      cellCells(4) = mesh%CellToCell(Flip(orthCellFaces(1)),cell)
      cellCells(5) = mesh%CellToCell(Flip(orthCellFaces(2)),cell)
      if (any(mesh%IsCellInternal(cellCells(2:5)))) then

        ! ----------------------
        ! This is a regular cell, 
        ! generate the VTK tetrahedrons around it.
        ! ----------------------
        call MakeVtkCell(cell, cellCells)
        cellCells(1) = mesh%CellToCell(Flip(cellFace),cell)
        if (mesh%IsCellInternal(cellCells(1))) then
          call MakeVtkCell(cell, cellCells)
        end if

      else

        ! ----------------------
        ! This is a hanging cell, 
        ! generate the VTK tetrahedrons around the neighbour.
        ! ----------------------
        cellCells(2) = mesh%CellToCell(orthCellFaces(1),cellCells(1))
        cellCells(3) = mesh%CellToCell(orthCellFaces(2),cellCells(1))
        cellCells(4) = mesh%CellToCell(Flip(orthCellFaces(1)),cellCells(1))
        cellCells(5) = mesh%CellToCell(Flip(orthCellFaces(2)),cellCells(1))
        call MakeVtkCell(cell, cellCells)

      end if

    end if

  end do

  ! ----------------------
  ! Print the VTK mesh statistics.
  ! ----------------------
  call PrintLog('')
  call PrintLog('-=-=-=-=-=-=-=-')
  call PrintLog('VTK mesh statistics:')
  call PrintLog('-=-=-=-=-=-=-=-')
  call PrintLog(' * Number of VTK cells: '//I2S(vtkMesh%NumVtkCells))
  call PrintLog('')

contains
  subroutine PushVtkCell(vtkCellNodes)
    integer(ip), intent(in) :: vtkCellNodes(:)

    vtkMesh%NumVtkCells = vtkMesh%NumVtkCells + 1
    vtkMesh%VtkCellToNode(:,vtkMesh%NumVtkCells) = vtkCellNodes

  end subroutine PushVtkCell
  subroutine MakeVtkCell(cell, cellCells)
    integer(ip), intent(in) :: cell, cellCells(:)

    if (mesh%NumDims == 2) then

      if (mesh%IsCellInternal(cellCells(2))) then
        call PushVtkCell([cell,cellCells(1),cellCells(2)])
      end if
      if (mesh%IsCellInternal(cellCells(3))) then
        call PushVtkCell([cell,cellCells(1),cellCells(3)])
      end if

    else if (mesh%NumDims == 3) then

      if (all(mesh%IsCellInternal(cellCells(2:3)))) then
        call PushVtkCell([cell, cellCells(1), cellCells(2:3)])
      end if
      if (all(mesh%IsCellInternal(cellCells(3:4)))) then
        call PushVtkCell([cell, cellCells(1), cellCells(3:4)])
      end if
      if (all(mesh%IsCellInternal(cellCells(4:5)))) then
        call PushVtkCell([cell, cellCells(1), cellCells(4:5)])
      end if
      if (all(mesh%IsCellInternal(cellCells(2:5:3)))) then
        call PushVtkCell([cell, cellCells(1), cellCells(2:5:3)])
      end if

    end if

  end subroutine MakeVtkCell
end subroutine IO_InitVtkMesh

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Print mesh in Legacy VTK '.vtk' format.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine IO_WriteToUnstructuredVTK(mesh, file, fields)
  class(tMesh), intent(inout) :: mesh
  type(IOList), intent(in), optional :: fields
  character(len=*), intent(in) :: file
  
  integer(ip) :: unit
  integer(ip) :: cell, vtkCell
  class(IOListItem), pointer :: item

  logical, parameter :: binary = .true., singleReals = .true.

  class(tVtkMesh), allocatable, save :: vtkMesh
  class(tEncoder), allocatable :: encoder

  if ((mesh%NumDims /= 2).and.(mesh%NumDims /= 3)) then
    error stop 'Only 2D/3D meshes can be printed to Legacy VTK'
  end if

  if (.not.allocated(vtkMesh)) then
    allocate(vtkMesh)
    call IO_InitVtkMesh(vtkMesh, mesh)
  end if

  call PrintLog('')
  call PrintLog('-=-=-=-=-=-=-=-')
  call PrintLog('Write to unstructured VTK file: '//file)
  call PrintLog('-=-=-=-=-=-=-=-')
  call PrintLog('')

  ! ----------------------
  ! Open file, initialize encoder and print the header.
  ! ----------------------
  open(newunit=unit, file=file, access='stream', status='replace')
  write(unit) '# vtk DataFile Version 3.0', endl
  write(unit) '# StormRuler unstructured VTK writer.', endl
  if (binary) then
    write(unit) 'BINARY', endl
    encoder = tBinaryEncoder(endianness='big', singleReals=singleReals)
  else
    write(unit) 'ASCII', endl
    encoder = tAsciiEncoder()
  end if
  write(unit) 'DATASET UNSTRUCTURED_GRID', endl
  write(unit) endl

  ! ----------------------
  ! Write VTK points (as centers of the mesh cells).
  ! ----------------------
  write(unit) 'POINTS ', I2S(mesh%NumCells), ' float', endl
  do cell = 1, mesh%NumCells
    associate(r => mesh%CellCenter(cell))
      write(unit) encoder%Encode(r, paddedSize=3)
    end associate
  end do
  write(unit) encoder%Leftover(), endl, endl

  ! ----------------------
  ! Write VTK cells (as in the VTK mesh).
  ! ----------------------
  write(unit) 'CELLS ', &
    & I2S(vtkMesh%NumVtkCells), ' ', &
    & I2S(vtkMesh%NumVtkCells*(vtkMesh%NumVtkCellNodes + 1)), endl
  do vtkCell = 1, vtkMesh%NumVtkCells
    write(unit) encoder%Encode(vtkMesh%NumVtkCellNodes)
    write(unit) encoder%Encode(vtkMesh%VtkCellToNode(:,vtkCell) - 1)
  end do
  write(unit) encoder%Leftover(), endl, endl

  write(unit) 'CELL_TYPES ', I2S(vtkMesh%NumVtkCells), endl
  do vtkCell = 1, vtkMesh%NumVtkCells
    write(unit) encoder%Encode(int(vtkMesh%VtkCellType, kind=ip))
  end do
  write(unit) encoder%Leftover(), endl, endl

  ! ----------------------
  ! Write fields.
  ! ----------------------
  if (present(fields)) then
    write(unit) 'POINT_DATA ', I2S(mesh%NumCells), endl
    item => fields%first
    do while(associated(item))
      select type(item)
        ! ----------------------
        ! Scalar field.
        ! ----------------------
        class is(IOListItem$0)
          write(unit) 'SCALARS ', item%name, ' float 1', endl
          write(unit) 'LOOKUP_TABLE default', endl
          write(unit) encoder%Encode(item%values(:mesh%NumCells))

        ! ----------------------
        ! Vector field.
        ! ----------------------
        class is(IOListItem$1)
          write(unit) 'VECTORS ', item%name, ' float', endl
          write(unit) encoder%Encode(item%values(:,:mesh%NumCells), paddedSize=3)

        end select
      write(unit) encoder%Leftover(), endl, endl
      item => item%next
    end do
  end if

  ! ----------------------
  ! Close file and exit.
  ! ----------------------
  close(unit)

end subroutine IO_WriteToUnstructuredVTK

end module StormRuler_IO_VTK
