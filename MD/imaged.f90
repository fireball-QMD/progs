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

 
! imaged.f90
! Program Description
! ===========================================================================
!       This routine will translate the atom positions during an MD 
! simulation in order to keep atoms in the central cell always.
!
! ===========================================================================
! Code written by:
! Kurt Glaesemann
! Lawrence Livermore National Laboratory
! 7000 East Ave.
! Livermore, CA 94550
! Office Telephone (925) 423-1579
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine imaged (icluster, iimage, itime_step, nstepi)
        use dimensions
        use configuration
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
        integer, intent (in) :: icluster
        integer, intent (in) :: iimage
        integer, intent (in) :: itime_step
        integer, intent (in) :: nstepi
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ibox

        real from_zero
        real shifted_from_zero
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Transform coordinates to box centered on (0,0,0)  (note: this is shifted
! if ishiftO = 1)
        if (iimage .gt. 0 .and. icluster .eq. 0) then 
                
!  otherwise we would always image first step
         if (mod(itime_step - nstepi, iimage) .eq. 0) then
          do iatom = 1, natoms
           from_zero = ratom(1,iatom)**2 + ratom(2,iatom)**2 + ratom(3,iatom)**2

! FIXME  - there is a 'GO TO' here that we need to eliminate
! Only need six boxes, since others are multiples 
           do ibox = 1, 6
514         shifted_from_zero = (ratom(1,iatom) + xl(1,ibox))**2             &
     &                         + (ratom(2,iatom) + xl(2,ibox))**2            &
     &                         + (ratom(3,iatom) + xl(3,ibox))**2
            if (shifted_from_zero .lt.  from_zero) then
             ratom(:,iatom) = ratom(:,iatom) + xl(:,ibox)
             xdot(0,:,iatom) = xdot(0,:,iatom) + xl(:,ibox)
             ximage(:,iatom) = ximage(:,iatom) - xl(:,ibox) 
             from_zero = shifted_from_zero
             go to 514

! Check again, just in case we jumped two boxes 
            end if
           end do
          end do
         end if
        end if
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
