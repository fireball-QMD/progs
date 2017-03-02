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

! Program Description
! ===========================================================================
!       This resets the temperature of the nose-hoover chain.
!
! ===========================================================================
! Code written by:
! J. Keith
! ===========================================================================
!
! Program Declaration
! ===========================================================================
subroutine resetNHC(natoms, T_want, T_wantPrev)
use noseHoover
use constants_fireball
implicit none
! Passed variables
integer, intent(in) :: natoms
real, intent(in) :: T_want, T_wantPrev
! Procedure
kT=kb*T_want
gkT=natoms*3*kT
Q_i = Q_i*T_want/T_wantPrev
v_xi = v_xi*sqrt(T_want/T_wantPrev)
!print*,'T_want/T_wantPrev',T_want/T_wantPrev
!print*,'kT, gkT, Q_i, v_xi',kT,gkT,Q_i,v_xi

end subroutine
