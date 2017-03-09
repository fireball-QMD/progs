! copyright info:
!
!                             @Copyright 2002
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

 
! recoverC.f90
! Program Description
! ===========================================================================
!       This subroutine takes a 1D list of integrals and generates a 2x2 
! array with respect to the shells.
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
        subroutine recoverC (n1, n2, hlist, dhlist, hbox, dhbox)
        use dimensions
        use interactions
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: n1, n2

        real, intent(in), dimension (ME2c_max) :: hlist 
        real, intent(in), dimension (ME2c_max) :: dhlist 

! Output
        real, intent(out), dimension (nsh_max, nsh_max) :: hbox
        real, intent(out), dimension (nsh_max, nsh_max) :: dhbox

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer index
        integer indexcoulomb
        integer ii
        integer kk
 
! Procedure
! ===========================================================================
! The algorithm is simply to go along the vector, and fill in the first row 
! of a matrix. Once the nssh(2) elements are absorbed, we increment the row 
! index ii by one, and reset the column index kk back to one.
        ii = 1
        kk = 0

! Loop over all the non-zero integrals for this interaction:
! n1 = number of shells on atom 1, and n2 = number of shells on atom2.
        indexcoulomb = n1*n2
        do index = 1, indexcoulomb
         kk = kk + 1
         hbox(ii,kk) = hlist(index)
         dhbox(ii,kk) = dhlist(index)
         if (mod(index,n2) .eq. 0) then
          ii = ii + 1
          kk = kk - n2
         end if
        end do
 
! Format Statements
! ===========================================================================
 
        return
        end
 
