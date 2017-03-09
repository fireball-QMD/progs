! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
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
!       This routine writes out the acceleration to the file - acfile.
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
        subroutine writeout_ac (acfile, itime_step, time, imass, nzx)
        use dimensions
        use configuration
        use constants_fireball
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itime_step

        integer, intent (in), dimension (natoms) :: imass
        integer, intent (in), dimension (nspecies) :: nzx
 
        real, intent (in) :: time

        character (len=30), intent (in) :: acfile
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer in1
        integer iorder
        integer i
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        write (*,*) ' Writing accelerations to the acfile. '
        open (unit = 64, file = acfile, status = 'unknown')
        write (64,100) natoms
        do iorder = 2, 5
         write (64,101) iorder, itime_step, time
         do iatom = 1, natoms
          in1 = imass(iatom)
          if (iorder .le. gear_order) then
            write (64,102) nzx(in1), (xdot(iorder,i,iatom),i=1,3)
          else
            write (64,102) nzx(in1), 0.0d0, 0.0d0, 0.0d0
          end if
         end do
        end do
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (2x, i4)
101     format (2x, i1, 2x, i6, 2x, f9.3)
102     format (2x, i2, 3(2x,f16.9))

 
        close (unit = 64)
        return
      end subroutine writeout_ac
