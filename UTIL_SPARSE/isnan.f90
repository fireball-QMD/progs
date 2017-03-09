! copyright info:
!
!                             @Copyright 2001
!                            Fireball Committee
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

 
! isnan.f90
! Program Description
! ===========================================================================
!       This is a logical function which returns a .true. if the input is
! NaN case.
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
        logical function isnan (x)
        implicit none
!       include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent (in) :: x
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        logical result
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
!       if (mpi_whatever_real .eq. MPI_REAL) then
         result = (x .ge. z'7fc00000' .and. x .le. z'7fffffff')              &
     &            .or. (x .le. z'ffc00000' .and. x .ge. z'ffffffff')
!       else
!        result =                                                            &
!    &    (x .ge. z'7ff8000000000000' .and. x .le. z'7fffffffffffffff')      &
!    &     .or. (x .le. z'fff8000000000000' .and. x .ge. z'ffffffffffffffff')
!       end if
        isnan = result

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end

