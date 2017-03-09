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


! invert3x3.f90
! Program Description
! ===========================================================================
!       This subroutine will calculate the inverse of a 3x3 matrix. The
! routine was originally written into the Sankey-Niklewski code by Alex A.
! Demkov and was converted to FORTRAN 90.
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
        subroutine invert3x3 (amatrix, ainverse)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input:
        real, intent (in), dimension (3, 3) :: amatrix
 
! Output:
        real, intent (out), dimension (3, 3) :: ainverse
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        real determinant
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        determinant =  &
     &    amatrix(1,1)*(amatrix(2,2)*amatrix(3,3) - amatrix(2,3)*amatrix(3,2))&
     &  + amatrix(1,2)*(amatrix(3,1)*amatrix(2,3) - amatrix(2,1)*amatrix(3,3))&
     &  + amatrix(1,3)*(amatrix(2,1)*amatrix(3,2) - amatrix(3,1)*amatrix(2,2))
 
! Now calculate inverse of inertia tensor
! ainv(i,j) = (-1**(i+j)) * cofactor(j,i) / det(a)
 
        if (abs(determinant) .gt. 1.0d-5) then
         ainverse(1,1) =    &
     &    (amatrix(2,2)*amatrix(3,3) - amatrix(3,2)*amatrix(2,3))/determinant
         ainverse(2,1) =    &
     &    - (amatrix(2,1)*amatrix(3,3) - amatrix(3,1)*amatrix(2,3))/determinant
         ainverse(3,1) =    &
     &    (amatrix(2,1)*amatrix(3,2) - amatrix(3,1)*amatrix(2,2))/determinant
         ainverse(1,2) =    &
     &    - (amatrix(1,2)*amatrix(3,3) - amatrix(3,2)*amatrix(1,3))/determinant
         ainverse(2,2) =    &
     &    (amatrix(1,1)*amatrix(3,3) - amatrix(3,1)*amatrix(1,3))/determinant
         ainverse(3,2) =    &
     &    - (amatrix(1,1)*amatrix(3,2) - amatrix(3,1)*amatrix(1,2))/determinant
         ainverse(1,3) =    &
     &    (amatrix(1,2)*amatrix(2,3) - amatrix(1,3)*amatrix(2,2))/determinant
         ainverse(2,3) =    &
     &    - (amatrix(1,1)*amatrix(2,3) - amatrix(1,3)*amatrix(2,1))/determinant
         ainverse(3,3) =    &
     &    (amatrix(2,2)*amatrix(1,1) - amatrix(1,2)*amatrix(2,1))/determinant
 
        else
         ainverse = 0.0d0
         write (*,*) ' ********* WARNING ********* '
         write (*,*) ' The determinant of the amatrix in invert3x3 is '
         write (*,*) ' equal to zero. Be careful and make sure that '
         write (*,*) ' you really want to continue. '
        end if
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
