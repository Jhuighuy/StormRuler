module helpers
  use StormRuler_Mesh
  implicit none
  real(8), parameter :: pi = 4*atan(1.0D0)
  Integer, Parameter :: Nx = 100, Ny = 100  
  Real(8), Parameter :: Dx = 2*pi/Nx, Dy = 2*pi/Ny, Dt = Dx*Dx
  contains
  Function to_str(k)
    Integer, Intent(In) :: k
    Character(Len=20) :: to_str
    Write (to_str, *) k
    to_str = adjustl(to_str)
  End Function to_str  
  
  Subroutine print_grid2(C, l)
    Integer, Intent(In) :: l
    Real(8), Dimension(0:Nx+1,0:Ny+1), Intent(InOut) :: C
    Integer :: output
    Write(*,*) l
    Open(NewUnit=output, file='out/fields-'//Trim(to_str(l))//'.csv', Status='replace')
    Write (output, *) 'x,y,z,c'
    Block
      Integer :: Ix,Iy
      Do Ix=1,Nx
        Do Iy=1,Ny
          Write(output, '(E12.6,A,E12.6,A,E12.6,A,E12.6,A)') &
            Dx*Ix, ',', Dy*Iy, ',', 0.0, ',', C(Ix,Iy), ','
        End Do
      End Do
    End Block
    Close(output)
  End Subroutine print_grid2

  Subroutine print_mesh2(mesh, u, l)
    class(Mesh2D), intent(in) :: mesh
    Real(8), Dimension(:), Intent(in) :: u
    Integer, Intent(In) :: l
    Integer :: output
    Write(*,*) l
    Open(NewUnit=output, file='out/fields-'//Trim(to_str(l))//'.csv', Status='replace')
    Write (output, *) 'x,y,z,c'
    Block
      integer :: iCell
      do iCell = 1, mesh%NumAllCells
        associate(r => mesh%CellCenter(iCell,:))
          Write(output, '(E12.6,A,E12.6,A,E12.6,A,E12.6,A)') &
            r(1), ',', r(2), ',', r(3), ',', u(iCell), ','
        end associate
      end do
    End Block
    Close(output)
  End Subroutine print_mesh2
  Subroutine print_mesh3(mesh, u,p,c, l)
    class(Mesh2D), intent(in) :: mesh
    Real(8), Dimension(:,:), Intent(in) :: u
    Real(8), Dimension(:), Intent(in) :: p,c
    Integer, Intent(In) :: l
    Integer :: output
    Write(*,*) l
    Open(NewUnit=output, file='out/fields-'//Trim(to_str(l))//'.csv', Status='replace')
    Write (output, *) 'x,y,z,vx,vy,p,c'
    Block
      integer :: iCell
      do iCell = 1, mesh%NumAllCells
        associate(r => mesh%CellCenter(iCell,:))
          Write(output, '(E12.6,A,E12.6,A,E12.6,A,E12.6,A,E12.6,A,E12.6,A,E12.6,A)') &
            r(1), ',', r(2), ',', r(3), ',', u(1,iCell), ',', u(2,iCell), ',', p(iCell), ',', c(iCell), ','
        end associate
      end do
    End Block
    Close(output)
  End Subroutine print_mesh3
end module helpers
  
program nsch
  use helpers
  use StormRuler_Mesh
  use StormRuler_CahnHilliard
  use StormRuler_NavierStokes
  implicit none

  Integer :: L,M,iCell
  Real(8), Dimension(:), Allocatable :: C,p,s
  Real(8), Dimension(:,:), Allocatable :: w,v,f
  type(CahnHilliardParams) :: CH
  class(Mesh2D), allocatable :: mesh
  
  
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
                                                      

  allocate(mesh)
  mesh%dt = dt
  call mesh%InitRect(dx,nx,.true.,dy,ny,.true.,20)

  !allocate(u(1:mesh%NumAllCells))
  !allocate(v(1:mesh%NumAllCells))
  allocate(c(1:mesh%NumAllCells))
  allocate(s(1:mesh%NumAllCells))
  allocate(p(1:mesh%NumAllCells))
  allocate(w(1:2,1:mesh%NumAllCells))
  allocate(v(1:2,1:mesh%NumAllCells))
  allocate(f(1:2,1:mesh%NumAllCells))

  !do iCell = 1, mesh%NumAllCells
  !  associate(r => mesh%CellCenter(iCell,:))
  !    !if (r(1)>=0.25.and.r(1)<=0.75) v(:,iCell)=[1,0]
  !    p(iCell) = -0.25_dp*(cos(2*r(1))+cos(2*r(2)))
  !    v(:,iCell) = [+cos(r(1))*sin(r(2))&
  !           ,-sin(r(1))*cos(r(2))]
  !  end associate
  !End Do

  do iCell = 1, mesh%NumAllCells
    associate(r => mesh%CellCenter(iCell,:)-[pi,pi])
      p(iCell) = 1
      c(iCell) = 1
      v(:,iCell) = [0.0,0.0]
      f(:,iCell) = [0.0,-1.0]
      if ((abs(r(1)) < 1).and.((abs(r(2)) < 1))) c(iCell) = -1
    end associate
  End Do

  CH%EpsSqr = 0.01D0
  Call print_mesh3(mesh, v,p,c, 0)
  Do L=1,200
    Do M=1,10
      Call CahnHilliard_ImplicitSchemeStep(mesh, c,s,v, CH)
      Call NavierStokes_PredictVelocity(mesh,v,p,c,s,f,0.1_dp,1.0_dp,0.0_dp)
    End Do
    Call print_mesh3(mesh, v,p,c, L)
  End Do

  !CH%EpsSqr = 0.01D0
  !Block
  !  do iCell = 1, mesh%NumAllCells
  !    Call RANDOM_NUMBER(C(iCell))
  !    C(iCell) = 2*C(iCell) - 1
  !  End Do
  !End Block
  !Call print_mesh2(mesh, C, 0)
  !Do L=1,200
  !  Do M=1,1
  !    Call CahnHilliard_ImplicitSchemeStep(mesh, C,s,v, CH)
  !  End Do
  !  Call print_mesh2(mesh, C, L)
  !End Do

end program nsch
