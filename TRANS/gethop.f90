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


! doscentros.f90
! Program Description
! ===========================================================================
!      This subroutine calculates the (two-center) matrix elements (mu,nu).
! There used to be five different routines that did this for all of the
! two-center interactions - doscentros.f, dosxcatm.f, dosxcontop.f,
! dosenatm.f, and dosenontop.f.  These have now all been reduced to one
! routine in order to make Fireball more lean.
!
!      This routine also calculates the derivative with respect to the
! position of the orbital of the BRA.
!
! ===========================================================================
! Original code written by Jose Ortega.
!
! Code rewritten by:
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
        subroutine gethop (in1, in2, ind1, ind2, distance, eps, sx)

        use dimensions
        use interactions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: in1
        integer, intent (in) :: in2
        integer, intent (in) :: ind1
        integer, intent (in) :: ind2
        real, intent (inout) :: distance
        real, intent (in), dimension (3, 3) :: eps
 
! Output
        real, intent (out), dimension (numorb_max, numorb_max) :: sx
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer inu
        integer index
 
! -slist = output list of matrix elements
        real, dimension (ME2c_max) :: slist

        real, dimension (numorb_max,numorb_max) :: sm

! Procedure
! ===========================================================================
! Initialize sm, scam and sx to zero.

        sm = 0.0d0
        sx = 0.0d0

! This subroutine calls the subroutine intrp1d as needed to find the value of
! the matrix elements for any given atomic separation distance.
! -slist = output list of matrix elements
        do index = 1, index_max2c(in1,in2)
         call interpolate_hop (index, ind1, ind2, distance, slist(index))
        end do
 
! Now recover sm which are two-dimensional arrays from
! slist and dslist which are one-dimensional arrays.
        call recover_2c (in1, in2, slist, sm)
 
! Rotate sm into crystal-coordinates: sm --> sx
        call rotate_fb (in1, in2, eps, sm, sx)
 
! Format Statements
! ===========================================================================
 
        return
      end subroutine gethop
