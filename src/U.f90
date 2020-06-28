!---------------------------------------------------------------------
! U
!	- module containing subrotuines for evaluating attractive
!	  potential
!---------------------------------------------------------------------
module U
  implicit none

contains 
!---------------------------------------------------------------------
! U_rA
!	- calculates attractive potential as a function of distance
!	  and angle 
!---------------------------------------------------------------------
! r		: real*8, radius 
! chrg		: int, charge of ion
! al		: 2D real*8, polarizability tensor
! ang		: real*8, angle
! dip		: 1D real*8, dipole vector
real(kind=8) function U_rA(r,chrg,al,ang,dip)
  implicit none
  integer, intent(in) :: chrg
  real(kind=8), intent(in) :: r,al,ang,dip 
  real(kind=8), dimension(3), intent(in) :: dip
  real(kind=8), dimension(3,3), intent(in) :: al
  real(kind=8) :: temp,a
  temp = 0.d0
  !just have the no dipole case for now
  a = (al(1,1) + al(2,2) + al(3,3))/3.0d0 
  temp = -0.5d0*chrg**2*a/r**4.0d0
  U_rA = temp
end function U_rA
!---------------------------------------------------------------------

end module U
