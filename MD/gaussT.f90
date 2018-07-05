! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Motorola, Physical Sciences Research Labs - Alex Demkov
! University of Regensburg - Juergen Fritsch
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola - Jun Wang
! Ohio University - Dave Drabold

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


! gaussT.f90
! Program Description
! ===========================================================================
!       This routine renormalizes the forces for each atom, such that 
! a constant temperature in the system is maintained.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine gaussT (natoms, vatom, xmass, T_want, ftot) 
        use constants_fireball
        use dimensions
        use options, only : verbosity 
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms

        real, intent (in) :: T_want

        real, intent (in), dimension (3, natoms) :: vatom
        real, intent (in), dimension (natoms) :: xmass

! Output
        real, intent (inout), dimension (3, natoms) :: ftot
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom

        real denominator
        real xksi
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! We are doing gaussian constant T dynamics.
!        write (*,*) '  '
!        write (*,*) ' We are doing constant temperature dynamics! '

! Calculate xksi (See Evans et al. P.R.A 28, 1016 (1983).
! xksi = SUM(1:3,1:natoms) ftot(ix,iatom)*v(ix,iatom)/(3/2*n*kB*T)
        denominator = 2.0d0*(3.0d0/2.0d0)*natoms*T_want/kconvert
        xksi = 0.0d0
        do iatom = 1, natoms
         xksi = xksi + ftot(1,iatom)*vatom(1,iatom)                          &
     &               + ftot(2,iatom)*vatom(2,iatom)                          &
     &               + ftot(3,iatom)*vatom(3,iatom)
        end do
        if (denominator .lt. 1.0d-3) then
         write (*,*) ' T_want = ', T_want
         write (*,*) ' This value will not work. Reset! '
         stop
        end if
        xksi = xksi/denominator

! fovermp is a conversion factor; velocity in A/fs.
        do iatom = 1, natoms
         ftot(:,iatom) =                                                     &
     &    ftot(:,iatom) - xksi*xmass(iatom)*vatom(:,iatom)/fovermp
        end do
        if (verbosity .ge. 2) then
          write (*,*) ' xksi = ', xksi, ' RENORMALIZE ftot! '
          do iatom = 1, natoms
           write (*,100) iatom, ftot(:,iatom)
          end do
          write (*,*) '  '
        end if ! verbosity

! Deallocate Arrays
! ===========================================================================
100     format (2x, ' iatom = ', i4, ' ftot   = ', 3e14.6)
 
! Format Statements
! ===========================================================================
 
        return
        end
