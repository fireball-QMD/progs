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
! fireball-qmd is a free (GPLv3) open project.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


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
