!---------------------------------------------------------------------
! TS	
!	- contains subroutine for TS properties
!---------------------------------------------------------------------
module TS
  use P
  implicit none

contains
!---------------------------------------------------------------------
! TS_V_nodpl
!	- calculates the effective potential given an energy and 
!	  seperation, with a no dplole potential
!---------------------------------------------------------------------
! enr		: real*8, energy
! b		: real*8, initial vertical seperation
! al		: 2D real*8, polarizability tensor
! chrg		: int, charge of ion
real(kind=8) function TS_V_nodpl(enr,b,al,chrg)
  implicit none
  integer, intent(in) :: chrg
  real(kind=8), intent(in) :: enr,b
  real(kind=8), dimension(3,3), intent(in) :: al
  real(kind=8) :: g
  g = chrg**2.0d0*(al(1,1) + al(2,2) + al(3,3))/3.0d0
  TS_V_nodpl = 0.5*enr**2.0d0*b**4.0d0/(g)
end function TS_V_nodpl

!---------------------------------------------------------------------
! TS_omg_nodpl
!	- calculates the effective potential given an energy and 
!	  seperation, with a nodplole potential
!---------------------------------------------------------------------
! enr		: real*8, energy
! b		: real*8, initial vertical seperation
! al		: 2D real*8, polarizability tensor
! chrg		: int, charge of ion
! mss		: real*8, mass
real(kind=8) function TS_om_nodpl(enr,b,al,chrg,mss)
  implicit none
  integer, intent(in) :: chrg
  real(kind=8), intent(in) :: enr,b,mss
  real(kind=8), dimension(3,3), intent(in) :: al
  real(kind=8) :: g
  g = chrg**2.d0*(al(1,1) + al(2,2) + al(3,3))/3.d0
  TS_om_nodpl = 2.d0*sqrt(enr**3.d0 * b**6.d0/(mss*g**2.d0))
end function TS_om_nodpl

!---------------------------------------------------------------------
! TS_bc_nodpl
!	- given energy and parameters, calculate an effective max
!	  b
!---------------------------------------------------------------------
! enr		: real*8, energy
! al		: 2D real*8, polarizability tensor
! chrg		: int, charge of ion
! mss		: real*8, mass
! Pmin		: real*8, minimum probability to consider
real(kind=8) function TS_bc_nodpl(enr,al,chrg,mss)
  implicit none
  integer, intent(in) :: chrg
  real(kind=8), intent(in) :: enr,mss
  real(kind=8), dimension(3,3), intent(in) :: al
  real(kind=8) :: g,Pmin,bc,delb,Peb
  Pmin = 1.0d-6
  g = chrg**2.d0*(al(1,1) + al(2,2) + al(3,3))/3.d0
  bc = 0.d0 
  delb = 1.d-1
  Peb = P_Eb_nodpl(enr,bc,g,mss)
  open(file='test',unit=100,status='replace')
  !go to above bc
  do while (Peb .gt. Pmin)
    bc = bc + delb
    Peb = P_Eb_nodpl(enr,bc,g,mss)
    write(100,*) bc,Peb
  end do 
  !get to bc
  delb = 1.d-3
  do while (Peb .lt. Pmin)
   bc = bc - delb
   Peb = P_Eb_nodpl(enr,bc,g,mss)
    write(100,*) bc,Peb
  end do 
  TS_bc_nodpl = bc
  close(unit=100)
end function TS_bc_nodpl
!---------------------------------------------------------------------
end module TS
!---------------------------------------------------------------------
