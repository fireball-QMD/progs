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


! make_mu2shell.f90
! Program Description
! ===========================================================================
!       This routine determines the shell number of the in1'th atomtype and
! stores that information in the array mu2shell.
!
! ===========================================================================
! Code written by:
! Otto F. Sankey (visiting)
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
        subroutine make_mu2shell (nspecies)
        use dimensions
        use interactions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nspecies
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer in1
        integer issh
        integer imumin, imumax
        integer Lvalue
        integer num_Mvalues
 
! Allocate Arrays
! ===========================================================================
        allocate (mu2shell(numorb_max, nspecies))
 
! Procedure
! ===========================================================================
! The extended hubbard subroutines need array mu2shell(imu,in1) which tells
! the shell number of the imu'th orbital of the in1'th atomtype.
        do in1 = 1, nspecies
         imumin = 0
         imumax = 0
         do issh = 1, nssh(in1)
          Lvalue = lssh(issh,in1)
          num_Mvalues = 2*Lvalue + 1
          imumax = imumin + num_Mvalues
          imumin = imumin + 1
          do imu = imumin, imumax
           mu2shell(imu,in1) = issh
          end do
          imumin = imumax
         end do
        end do
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
