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
  & ErrorStop, PrintLog, PrintWarning, &
  & MergeString, Flip, I2S, R2S

use StormRuler_Mesh, only: tMesh

use StormRuler_IO_Stream, only: tOutputStream, tUnitOutputStream
use StormRuler_IO_Stream_Base64, only: tBase64OutputStream
use StormRuler_IO_Writer, only: tWriter, tTextWriter, tBinaryWriter

use StormRuler_IO, only: IOList, IOListItem, @{IOListItem$$@|@0, 2}@

#$use 'StormRuler_Macros.fi'

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Print mesh in the simple VTK ('.vtk') format.
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine IO_WriteDenseStructuredVTK(mesh, file, fields)
  class(tMesh), intent(inout) :: mesh
  type(IOList), intent(in), optional :: fields
  character(len=*), intent(in) :: file

  integer(ip) :: unit
  integer(ip) :: x, y, cell, numCells
  class(IOListItem), pointer :: item

  logical, parameter :: binary = .true., singleReals = .true.

  class(tOutputStream), allocatable :: stream
  class(tWriter), allocatable :: writer

  call PrintLog('Write to VTS file: '//file)

  ! ----------------------
  ! Open file, initialize stream and writer.
  ! ----------------------
  open(newunit=unit, file=file, access='stream', status='replace')
  stream = tUnitOutputStream(unit)
  if (binary) then
    stream = tBase64OutputStream(stream)
    writer = tBinaryWriter(singleReals=singleReals)
  else
    writer = tTextWriter(separator=endl)
  end if

  ! ----------------------
  ! Write the header.
  ! ----------------------
  write(unit) '<?xml version="1.0"?>', endl
  write(unit) '<VTKFile type="StructuredGrid"'
  write(unit) ' version="0.1" '!, ' header_type="UInt64"'
#$if BIG_ENDIAN
  write(unit) ' byte_order="BigEndian">', endl
#$else
  write(unit) ' byte_order="LittleEndian">', endl
#$end if

  ! ----------------------
  ! Write the structured grid.
  ! ----------------------
  write(unit) '<StructuredGrid WholeExtent="'
  write(unit) '0 ', I2S(mesh%IndexBounds(1)), ' 0 ', I2S(mesh%IndexBounds(2)), ' 0 0'
  write(unit) '">', endl
  write(unit) '<Piece Extent="'
  write(unit) '0 ', I2S(mesh%IndexBounds(1)), ' 0 ', I2S(mesh%IndexBounds(2)), ' 0 0'
  write(unit) '">', endl

  write(unit) '<Points>', endl
  write(unit) '<DataArray'
  write(unit) ' type="', MergeString('Float32', 'Float64', singleReals), '"'
  write(unit) ' format="', MergeString('binary', 'ascii', binary), '"'
  write(unit) ' NumberOfComponents="3">'
  call stream%BeginWrite()
  do y = 0, mesh%IndexBounds(2)
    do x = 0, mesh%IndexBounds(1)
      call writer%Write(stream, 1.0_dp*[x, y, 1])
    end do
  end do
  call stream%EndWrite()
  write(unit) '</DataArray>', endl
  write(unit) '</Points>', endl

  ! ----------------------
  ! Write cell visibility.
  ! ----------------------
  write(unit) '<CellData>', endl
  write(unit) '<DataArray Name="vtkGhostLevels" type="UInt8"'
  write(unit) ' format="', MergeString('binary', 'ascii', binary), '">'
  call stream%BeginWrite()
  do y = 1, mesh%IndexBounds(2)
    do x = 1, mesh%IndexBounds(1)
      cell = mesh%IndexToCell(x,y)
      call writer%Write(stream, merge(0_bp, 17_bp, mesh%IsCellInternal(cell)))
    end do
  end do
  call stream%EndWrite()
  write(unit) '</DataArray>', endl

  ! ----------------------
  ! Write fields.
  ! ----------------------
  if (present(fields)) then
    item => fields%first
    do while(associated(item))
      write(unit) '<DataArray Name="', item%name, '"'
      write(unit) ' type="', MergeString('Float32', 'Float64', singleReals), '"'
      write(unit) ' format="', MergeString('binary', 'ascii', binary), '"'
      select type(item)
        ! ----------------------
        ! Scalar field.
        ! ----------------------
        class is(IOListItem$0)
          write(unit) '>', endl
          call stream%BeginWrite()
          do y = 1, mesh%IndexBounds(2)
            do x = 1, mesh%IndexBounds(1)
              cell = mesh%IndexToCell(x,y)
              if (mesh%IsCellInternal(cell)) then
                call writer%Write(stream, item%values(cell))
              else
                call writer%Write(stream, 0.0_dp)
              end if
            end do
          end do
          call stream%EndWrite()

        ! ----------------------
        ! Vector field.
        ! ----------------------
        class is(IOListItem$1)
          write(unit) ' NumberOfComponents="3">'
          call stream%BeginWrite()
          do y = 1, mesh%IndexBounds(2)
            do x = 1, mesh%IndexBounds(1)
              cell = mesh%IndexToCell(x,y)
              if (mesh%IsCellInternal(cell)) then
                call writer%Write(stream, item%values(:,cell), paddedSize=3)
              else
                call writer%Write(stream, [0.0_dp], paddedSize=3)
              end if
            end do
          end do
          call stream%EndWrite()

        end select
      write(unit) '</DataArray>', endl
      item => item%next
    end do
  end if

  ! ----------------------
  ! Write footer and close the file.
  ! ----------------------
  write(unit) '</CellData>', endl
  write(unit) '</Piece>', endl
  write(unit) '</StructuredGrid>', endl
  write(unit) '</VTKFile>', endl
  close(unit)

end subroutine IO_WriteDenseStructuredVTK

end module StormRuler_IO_VTK
