!---------------------------------------------------------------------
! SCLD
!	- calculates semiclassical TST
!	  density and sum of states for an ion-molecule 
!	  induced dipole
!---------------------------------------------------------------------
! numE		: int, number of energy bins 
! delE		: real*8, energy stepsize
! mss		: real*8, reduced mass
! al		: 2D real*8, polarizability tensor
! chrg		: int, ion's charge
! dpl		: real*8, dipole moment vector
program SCLD
  use TS
  implicit none
  integer :: numE,chrg
  real(kind=8) :: delE,mss,dip,mn,mi
  real(kind=8), dimension(3) :: dpl
  real(kind=8), dimension(3,3) :: al
  real(kind=8) :: NE,rhoE,bc,cm2au,amu2me
  real(kind=8) :: enr,g,pi
  integer :: i

  cm2au = 1.d0/2.194746313710d9
  amu2me = 1836.1526675d0
  pi = 3.1415926535897932

  open(file='SCLD.dat',unit=100,status='old')
  read(100,*) mi
  read(100,*) mn
  read(100,*) chrg
  read(100,*) al(1,1:3)
  read(100,*) al(2,1:3)
  read(100,*) al(3,1:3)
  read(100,*) dpl(1:3)
  read(100,*) numE,delE
  close(unit=100)

  !convert units
  mss = mi*mn/(mi + mn)
  mss = mss*amu2me

  !calculate g
  g = chrg**2.d0*(al(1,1) + al(2,2) + al(3,3))/3.d0
  
  !do i=0,numE-1
  !special for energy = 0. 
  open(file='SCLD.out',unit=200,status='replace')
  open(file='CLLD.out',unit=300,status='replace')
  do i=1,numE
    enr = i*delE*cm2au

    !calculate bc
    bc = TS_bc_nodpl(enr,al,chrg,mss)
    NE = P_E_nodpl(enr,bc,g,mss)

    write(200,*) i,i*delE,NE,bc
    write(300,*) i,i*delE,2.d0*pi*sqrt(2.d0*g/enr),(2.d0*g/enr)**0.25d0

  end do   
  close(unit=200)
  close(unit=300)

end program SCLD
!---------------------------------------------------------------------
