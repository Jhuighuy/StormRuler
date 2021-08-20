module helpers
use StormRuler_Mesh
implicit none
real(dp), parameter :: pi = 4*atan(1.0D0)
integer(ip), Parameter :: Nx = 100, Ny = 100  
Real(dp), Parameter :: Dx = 2*pi/Nx, Dy = 2*pi/Ny, Dt = Dx*Dx
contains

Function to_str(k)
  integer(ip), Intent(In) :: k
  Character(Len=20) :: to_str
  Write (to_str, *) k
  to_str = adjustl(to_str)
End Function to_str  

end module helpers
  
#$if True
program nsch
  use helpers
  use StormRuler_Helpers
  use StormRuler_Mesh
  use StormRuler_CahnHilliard
  use StormRuler_NavierStokes
  use StormRuler_IO
  use StormRuler_ConvParams

  implicit none

  integer(ip) :: L,M,iCell
  integer(ip), allocatable :: pixels(:,:)
  Real(dp), Dimension(:), Allocatable, target :: C,p,s
  Real(dp), Dimension(:,:), Allocatable, target :: w,v,f
  type(CahnHilliardParams) :: CH
  class(tMesh), allocatable :: mesh, mesh1
  type(IOList) :: fields
  integer(ip), allocatable :: colorToBCM(:)
  type(tConvParams) :: params

  ! ----------------------
  ! Print an epic banner.
  print *, ''
  print *, '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
  print *, '    _____ __                       ____        __           '
  print *, '   / ___// /_____  _________ ___  / __ \__  __/ ___  _____  '
  print *, '   \__ \/ __/ __ \/ ___/ __ `__ \/ /_/ / / / / / _ \/ ___/  '
  print *, '  ___/ / /_/ /_/ / /  / / / / / / _, _/ /_/ / /  __/ /      '
  print *, ' /____/\__/\____/_/  /_/ /_/ /_/_/ |_|\__,_/_/\___/_/       '
  print *, '                                                            '
  print *, '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
  print *, ''

  if (.false.) then  
    allocate(mesh1)
    call Load_PPM('test/Domain-100.ppm',pixels)
    colorToBCM = [PixelToInt([255,255,255]),PixelToInt([255,0,0])]
    call mesh1%InitFromImage2D(pixels,0,colorToBCM,3)
    call mesh1%PrintTo_Neato('test/c2c.dot')
    call mesh1%PrintTo_LegacyVTK('test/c2c.vtk')

    allocate(c(1:mesh1%NumAllCells))
    allocate(s(1:mesh1%NumAllCells))
    s(:) = 1.0_dp
    c(:) = 1.0_dp

    call fields%Add('c',c)
    call fields%Add('s',s)

    call SolvePoisson(mesh1,s,c)
    call mesh1%PrintTo_LegacyVTK('out/fld-'//I2S(0)//'.vtk',fields)
  end if

  if (.true.) then
    allocate(mesh)
    mesh%dt = dt
    call mesh%InitRect(dx,nx,.true.,dy,ny,.true.,20)
    call mesh%SetRange()

    allocate(c(1:mesh%NumAllCells))
    allocate(s(1:mesh%NumAllCells))
    allocate(p(1:mesh%NumAllCells))
    allocate(w(1:2,1:mesh%NumAllCells))
    allocate(v(1:2,1:mesh%NumAllCells))
    allocate(f(1:2,1:mesh%NumAllCells))

    call fields%Add('velocity',v)
    call fields%Add('pressure',p)
    call fields%Add('phase',c)

    do iCell = 1, mesh%NumAllCells
      associate(r => mesh%CellCenter(iCell)-[pi,pi])
        p(iCell) = 1.0_dp
        c(iCell) = 1.0_dp
        v(:,iCell) = [0.0,0.0]
        f(:,iCell) = [0.0,-1.0]
        if ((abs(r(1)) < 1).and.((abs(r(2)) < 1))) c(iCell) = -1
      end associate
    end Do

    CH%EpsSqr = 0.01D0
    call mesh%PrintTo_LegacyVTK('out/fld-'//I2S(0)//'.vtk',fields)
    Do L=1,200
      Do M=1,10
        Call CahnHilliard_ImplicitSchemeStep(mesh, c,s,v, CH)
        Call NavierStokes_PredictVelocity(mesh,v,p,c,s,f,0.1_dp,1.0_dp,0.0_dp)
      End Do
      call mesh%PrintTo_LegacyVTK('out/fld-'//I2S(L)//'.vtk',fields)
    End Do
  end if

end program nsch
#$end if
