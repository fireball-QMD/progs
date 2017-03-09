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


! rotated.f90
! Program Description
! ===========================================================================
!       This routine rotates a matrix derivative from molecular coordinates
! to crystal coordinates.
!
! The variable matm is the matrix-box in molecular coordinates
! The variable dmatm is a matrix which is the derivative of the matrix-box
! in molecular coordinates.
! The output is dmatx which is the derivative of the matrix-box in crystal
! coordinates.
!
! In the molecular coordinates we have atoms 1 and 2 (bondcharge) along the
! z-axis; the third atom is in the XZ-plane.
! The labelling of the orbitals is as follows:
!
!   S-shell :                s
!                            1
!
!   P-shell :           py   pz   px
!                       -1   0    1
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                 -2    -1    0    1      2
!
! ===========================================================================
! Original code written by Alex Demkov.
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
        subroutine rotated (in1, in2, eps, deps, matm, dmatm, dmatx)
        use dimensions
        use interactions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: in1, in2
 
        real, intent(in) :: deps (3, 3, 3)
        real, intent(in) :: eps (3, 3)
 
! Note that all dmatm arrays are vectors. That means the corresponding
! scalar derivatives were multiplied by -eta, which is dD/dr1, after the
! interpolation.
        real, intent(in) :: dmatm (3, numorb_max, numorb_max)
        real, intent(in) :: matm (numorb_max, numorb_max)
 
! Output
        real, intent(out) :: dmatx (3, numorb_max, numorb_max)
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer inu
        integer ix
 
        real ddmat (3, 5, 5)
        real dmat (5, 5)
        real dpmat (3, 3, 3)
        real pmat (3, 3)
 
! Procedure
! ===========================================================================
! Set rotational matrices and their respective derivatives.
        call twister (eps, dmat, pmat)
        call twisterd (eps, deps, ddmat, dpmat)
 
! Now we have to do the cases A, B and C, which  are added in makeDmat
! page 1 10/1/98 notes.
! term A: dAleft/dr*Aright*matrix
! term B: Aleft*dAright/dr*matrix
! term C Aleft*Aright*dmatrix/dr 
        call makeDmat (in1, in2, matm, dmatm,   &
     &                 dmat, pmat, ddmat, dpmat, dmatx)
 
! Format Statements
! ===========================================================================
 
        return
        end
 
