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

#$use 'StormRuler_Macros.fi'

use StormRuler_Consts, only: bp, ip, dp, endl

use StormRuler_Helpers, only: &
  & ErrorStop, PrintLog, PrintWarning, &
  & MergeString, Flip, I2S, R2S

use StormRuler_Mesh, only: tMesh

use StormRuler_IO_Writer, only: tWriter, tTextWriter, tBinaryWriter
use StormRuler_IO_Stream, only: tOutputStream, tUnitOutputStream
use StormRuler_IO_Stream_Base64, only: tBase64OutputStream
#$if HAS_ZLIB
use StormRuler_IO_Stream_ZLib, only: tZLibOutputStream
#$end if

use StormRuler_IO, only: IOList, IOListItem, @{IOListItem$$@|@0, 2}@

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

implicit none

!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< !!
!! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> !!

contains

!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
!! Write the mesh and fields in the VTK Image Data format ('.vti').
!! -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- !!
subroutine IO_WriteVtkImageData(mesh, file, fields)
  class(tMesh), intent(inout) :: mesh
  type(IOList), intent(in), optional :: fields
  character(len=*), intent(in) :: file

  integer(ip) :: unit
  integer(ip) :: x, y, cell, numPixels
  class(IOListItem), pointer :: item

  logical, parameter :: binary = .true., compressed = .true., singleReals = .false.

  class(tOutputStream), allocatable :: stream
  class(tWriter), allocatable :: writer

  call PrintLog('Write to VTI file: '//file)

  numPixels = mesh%IndexBounds(1)*mesh%IndexBounds(2)
  if (mesh%NumDims == 3) numPixels = numPixels*mesh%IndexBounds(3)

  ! ----------------------
  ! Open file, initialize stream and writer.
  ! ----------------------
  open(newunit=unit, file=file, access='stream', status='replace')
  stream = tUnitOutputStream(unit)
  if (binary) then
    stream = tBase64OutputStream(stream)
#$if HAS_ZLIB
    if (compressed) stream = tZLibOutputStream(stream)
#$end if
    writer = tBinaryWriter(singleReals=singleReals)
  else
    writer = tTextWriter(separator=endl)
  end if

  ! ----------------------
  ! Write the header.
  ! ----------------------
  write(unit) '<?xml version="1.0"?>', endl
  write(unit) '<VTKFile type="ImageData"'
  write(unit) ' version="1.0"'
  if (binary) then
    !! TODO: 64-bit header type.
    write(unit) ' header_type="UInt32"'
  end if
#$if BIG_ENDIAN
  write(unit) ' byte_order="BigEndian"'
#$else
  write(unit) ' byte_order="LittleEndian"'
#$end if
#$if HAS_ZLIB
  if (binary.and.compressed) then
    write(unit) ' compressor="vtkZLibDataCompressor"'
  end if
#$end if
  write(unit) '>', endl

  ! ----------------------
  ! Write the mesh header.
  ! ----------------------
  !! TODO: write true extents, origin and spacing.
  write(unit) '<ImageData WholeExtent="'
  write(unit) '0 ', I2S(mesh%IndexBounds(1)), ' 0 ', I2S(mesh%IndexBounds(2)), ' 0 0'
  write(unit) '" Origin="0 0 0" Spacing="1.0e-2 1.0e-2 1.0">', endl
  write(unit) '<Piece Extent="'
  write(unit) '0 ', I2S(mesh%IndexBounds(1)), ' 0 ', I2S(mesh%IndexBounds(2)), ' 0 0'
  write(unit) '">', endl

  ! ----------------------
  ! Write cell visibility.
  ! ----------------------
  write(unit) '<CellData>', endl
  write(unit) '<DataArray Name="vtkGhostType" type="UInt8"'
  write(unit) ' format="', MergeString('binary', 'ascii', binary), '">', endl
  call stream%BeginWrite()
  if (binary.and.(.not.compressed)) then
    call writer%Write(stream, numPixels)
  end if
  !! TODO: 3D case.
  do y = 1, mesh%IndexBounds(2)
    do x = 1, mesh%IndexBounds(1)
      cell = mesh%IndexToCell(x,y)
      call writer%Write(stream, merge(0_bp, 17_bp, mesh%IsCellInternal(cell)))
    end do
  end do
  call stream%EndWrite()
  write(unit) endl, '</DataArray>', endl

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
          if (binary.and.(.not.compressed)) then
            call writer%Write(stream, merge(4, 8, singleReals)*numPixels)
          end if
          !! TODO: 3D case.
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
          if (binary.and.(.not.compressed)) then
            call writer%Write(stream, 3*merge(4, 8, singleReals)*numPixels)
          end if
          !! TODO: 3D case.
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
      write(unit) endl, '</DataArray>', endl
      item => item%next
    end do
  end if

  ! ----------------------
  ! Write footer and close the file.
  ! ----------------------
  write(unit) '</CellData>', endl
  write(unit) '</Piece>', endl
  write(unit) '</ImageData>', endl
  write(unit) '</VTKFile>', endl
  close(unit)

end subroutine IO_WriteVtkImageData

end module StormRuler_IO_VTK
