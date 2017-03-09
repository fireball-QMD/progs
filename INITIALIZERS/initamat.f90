! copyright info:
!
!                             @Copyright 2003
!                            Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
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

 
! initamat.f90
! Program Description
! ===========================================================================
! Sets up the matrix amat, which is used by D-orbitals in twister(d)
! If we do not have any D-orbitals, then just return
!
! ===========================================================================
! Code written by:
! Kurt R. Glaesemann
! LLNL
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine initamat(nspecies)
        use constants_fireball
        use interactions
        implicit none

! Argument Declaration and Description
! ===========================================================================
        integer, intent(in) :: nspecies
 
! Local Variable Declaration and Description
! ===========================================================================
        integer in1
        integer issh

! Procedure
! ===========================================================================
! Initialize the a-matrices:  
        haveDorbitals = .false.
        do in1 = 1, nspecies
         do issh = 1, nssh(in1)
          if (lssh(issh,in1) .eq. 2) haveDorbitals = .true.
         end do
         do issh = 1, nsshPP(in1)
          if (lsshPP(issh,in1) .eq. 2) haveDorbitals = .true.
         end do
        end do

        amat(:,:,:) = 0.0d0

        if (.not. haveDorbitals) return

        amat(2,1,1) = 0.5d0
        amat(1,2,1) = 0.5d0
  
        amat(3,2,2) = 0.5d0
        amat(2,3,2) = 0.5d0
   
        amat(1,1,3) = - 1.0d0/(2.0d0*sqrt(3.0d0))
        amat(2,2,3) = - 1.0d0/(2.0d0*sqrt(3.0d0))
        amat(3,3,3) = 2.0d0/(2.0d0*sqrt(3.0d0))
   
        amat(3,1,4) = 0.5d0
        amat(1,3,4) = 0.5d0
 
        amat(1,1,5) = 0.5d0 
        amat(2,2,5) = - 0.5d0
 
        return
        end
