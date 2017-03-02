! copyright info:
!
!                             @Copyright 2005
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek
! Brigham Young University - Hao Wang

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang 
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! initNH.f90
! Program Description
! ===========================================================================
!       This writes out the NH hamiltonian--should be "somewhat" constant
!
! ===========================================================================
! Code written by:
! Keith
! ===========================================================================
!
! Program Declaration
! ===========================================================================
subroutine writeHNose(time,dt,natoms,xmass,ratom,vatom,etot) 
use noseHoover
implicit none
integer, intent(in):: natoms
real, intent(in) :: time,dt,etot
real, intent(in), dimension(natoms) :: xmass
real, intent(in), dimension(3,natoms) :: ratom,vatom
real enk,hn
integer i

call get_enk(natoms,xmass,vatom,enk)     
hn = etot + enk
do i = 1, nnos
 if (i.eq.1) then
  hn = hn + 0.5*v_xi(i)*v_xi(i)*Q_i(i) + gkT*xi(i)
 else
  hn = hn + 0.5*v_xi(i)*v_xi(i)*Q_i(i) + kT*xi(i)
 end if
end do

write(*,*)'NHC hamiltonian ',hn
if (mod(int((time-dt)/dt),2).eq.0) write(22,*)ratom(1,1),' ',vatom(1,1),' ',hn

end subroutine
