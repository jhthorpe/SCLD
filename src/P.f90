!---------------------------------------------------------------------
! P
!	- calculates probability	
!---------------------------------------------------------------------
module P
  implicit none

contains
!---------------------------------------------------------------------
! P_Eb_nodpl
!	- calculates transmission probability
!---------------------------------------------------------------------
! enr		: real*8, energy
! b		: real*8, initial offset
! g		: real*8, gamma
! mss		: real*8, reduced mass
real(kind=8) function P_Eb_nodpl(enr,b,g,mss)
  implicit none
  real(kind=8) :: pi
  real(kind=8), intent(in) :: enr,b,g,mss
  pi = 3.1415926535897932
!  P_Eb_nodpl = 1.d0/(1.d0 + exp(-1.d0*pi*sqrt(mss)* &
!                        (g/(enr**0.5d0*b**3.d0) &
!                         - 0.5d0*enr**0.5d0*b) ) )
  P_Eb_nodpl = 1.d0/(1.d0 + exp(pi*sqrt(mss)*&
               (0.5d0*b*enr**0.5d0 - g/(b**3.d0*enr**0.5d0))))
!  P_Eb_nodpl = 1.d0/(1.d0 + exp(-1.d0*pi*sqrt(mss)* &
!                        (g**2.d0/(enr**0.5d0*b**3.d0) &
!                         - 0.5d0*enr**0.5d0*b*g) ) )
end function P_Eb_nodpl
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! P_E_nodpl
!	- calculates transmission probability and integrates out b
!---------------------------------------------------------------------
! enr		: real*8, energy
! bc		: real*8, critical distance
! g		: real*8, gamma
! mss		: real*8, reduced mass
real(kind=8) function P_E_nodpl(enr,bc,g,mss)
  implicit none
  real(kind=8) :: pi
  real(kind=8), intent(in) :: enr,bc,g,mss
  real(kind=8) :: temp,b,db
  integer :: nb,i
  pi = 3.1415926535897932
  nb = 1d5
  temp = 0.d0
  db = bc/nb
  do i=0,nb-1
    b = i*db
    temp = temp + 2.d0*pi*b*P_Eb_nodpl(enr,b,g,mss)*db 
  end do 
  P_E_nodpl = temp
end function P_E_nodpl
!---------------------------------------------------------------------
end module P
!---------------------------------------------------------------------
