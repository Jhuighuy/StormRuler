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

#$use 'StormRuler_Params.fi'

use StormRuler_Parameters, only: dp, ip, i8

use StormRuler_Helpers, only:  &
  & ErrorStop, PrintLog, PrintWarning, & 
  & Flip, I2S, R2S

use StormRuler_Mesh, only: tMesh

use StormRuler_IO, only: IOList, IOListItem, @{IOListItem$$@|@0, 2}@
use StormRuler_IO_Stream, only: tOutputStream, &
  & tTextFileOutputStream, tBinaryFileOutputStream
use StormRuler_IO_Base64, only: tBase64OutputStream

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
  integer(i8) :: VtkCellType

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

  vtkMesh%VtkCellType = merge(5_i8, 10_i8, mesh%NumDims == 2)
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
        if (mesh%IsCellInternal(cellCells(2))) then
          call PushVtkCell([cell,cellCells(1),cellCells(2)])
        end if
        if (mesh%IsCellInternal(cellCells(3))) then
          call PushVtkCell([cell,cellCells(1),cellCells(3)])
        end if
        cellCells(1) = mesh%CellToCell(Flip(cellFace),cell)
        if (mesh%IsCellInternal(cellCells(1))) then
          if (mesh%IsCellInternal(cellCells(2))) then
            call PushVtkCell([cell,cellCells(1),cellCells(2)])
          end if
          if (mesh%IsCellInternal(cellCells(3))) then
            call PushVtkCell([cell,cellCells(1),cellCells(3)])
          end if
        end if

      else

        ! ----------------------
        ! This is a hanging cell, 
        ! generate the VTK triangles around the neighbour.
        ! ----------------------
        cellCells(2) = mesh%CellToCell(orthCellFaces(1),cellCells(1))
        cellCells(3) = mesh%CellToCell(Flip(orthCellFaces(1)),cellCells(1))
        if (mesh%IsCellInternal(cellCells(2))) then
          call PushVtkCell([cell,cellCells(1),cellCells(2)])
        end if
        if (mesh%IsCellInternal(cellCells(3))) then
          call PushVtkCell([cell,cellCells(1),cellCells(3)])
        end if

      end if

    else if (mesh%NumDims == 2) then

      orthCellFaces(1) = 2*mod((cellFace - 1)/2 + 1, 3) + 1
      orthCellFaces(2) = 2*mod((cellFace - 1)/2 + 2, 3) + 1
      cellCells(2) = mesh%CellToCell(orthCellFaces(1),cell)
      cellCells(3) = mesh%CellToCell(Flip(orthCellFaces(1)),cell)
      cellCells(4) = mesh%CellToCell(orthCellFaces(2),cell)
      cellCells(5) = mesh%CellToCell(Flip(orthCellFaces(2)),cell)

      error stop 'not implemented..'

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
  character, parameter :: endl = new_line('A')
  class(IOListItem), pointer :: item
  logical, parameter :: binary = .true.
  
  class(tVtkMesh), allocatable, save :: vtkMesh
  class(tOutputStream), allocatable :: stream

  if ((mesh%NumDims /= 2).and.(mesh%NumDims /= 3)) then
    error stop 'Only 2D/3D meshes can be printed to Legacy VTK'
  end if

  if (.not.allocated(vtkMesh)) then
    allocate(vtkMesh)
    call IO_InitVtkMesh(vtkMesh, mesh)
  end if
  if (binary) then
    allocate(tBinaryFileOutputStream :: stream)
  else
    allocate(tTextFileOutputStream :: stream)
  end if

  call PrintLog('')
  call PrintLog('-=-=-=-=-=-=-=-')
  call PrintLog('Write to unstructured VTK file: '//file)
  call PrintLog('-=-=-=-=-=-=-=-')
  call PrintLog('')

  ! ----------------------
  ! Open file, initialize stream and print the header.
  ! ----------------------
  open(newunit=unit, file=file, access='stream', status='replace')
  write(unit) '# vtk DataFile Version 3.0', endl
  write(unit) '# StormRuler unstructured VTK writer.', endl
  select type(stream)
    class is(tTextFileOutputStream)
      write(unit) 'ASCII', endl
      call stream%Init(unit)
    class is(tBinaryFileOutputStream)
      write(unit) 'BINARY', endl
      call stream%Init(unit)
  end select
  write(unit) 'DATASET UNSTRUCTURED_GRID', endl
  write(unit) endl

  ! ----------------------
  ! Write VTK nodes (as cell centers in the mesh).
  ! ----------------------
  write(unit) 'POINTS ', I2S(mesh%NumCells), ' double', endl
  do cell = 1, mesh%NumCells
    associate(r => mesh%CellCenter(cell))
      if (mesh%NumDims == 2) then
        call stream%Write([r, 0.0_dp])
      else
        call stream%Write(r)
      end if
    end associate
  end do
  call stream%Finalize()
  write(unit) endl, endl

  ! ----------------------
  ! Write VTK cells (as in the VTK mesh).
  ! ----------------------
  write(unit) 'CELLS ', I2S(vtkMesh%NumVtkCells), ' ', &
    & I2S(vtkMesh%NumVtkCells*(vtkMesh%NumVtkCellNodes + 1)), endl
  do vtkCell = 1, vtkMesh%NumVtkCells
    call stream%Write([vtkMesh%NumVtkCellNodes])
    call stream%Write(vtkMesh%VtkCellToNode(:,vtkCell) - 1)
  end do
  call stream%Finalize()
  write(unit) endl, endl
  write(unit) 'CELL_TYPES ', I2S(vtkMesh%NumVtkCells), endl
  do vtkCell = 1, vtkMesh%NumVtkCells
    call stream%Write([int(vtkMesh%VtkCellType, kind=ip)])
  end do
  call stream%Finalize()
  write(unit) endl, endl

  ! ----------------------
  ! Write fields.
  ! ----------------------
  if (present(fields)) then
    write(unit) 'POINT_DATA', ' ', I2S(mesh%NumCells), endl
    item => fields%first
    do while(associated(item))
      select type(item)
        ! ----------------------
        ! Scalar field.
        ! ----------------------
        class is(IOListItem$0)
          write(unit) 'SCALARS ', item%name, ' double 1', endl
          write(unit) 'LOOKUP_TABLE default', endl
          do cell = 1, mesh%NumCells
            call stream%Write([item%values(cell)])
          end do

        ! ----------------------
        ! Vector field.
        ! ----------------------
        class is(IOListItem$1)
          write(unit) 'VECTORS ', item%name, ' double', endl
          do cell = 1, mesh%NumCells
            call stream%Write([item%values(:,cell), 0.0_dp])
          end do

        end select
      write(unit) endl, endl
      item => item%next
    end do
  end if

  ! ----------------------
  ! Close file and exit.
  ! ----------------------
  close(unit)

end subroutine IO_WriteToUnstructuredVTK

end module StormRuler_IO_VTK
