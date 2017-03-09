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


! rotate.f90
! Program Description
! ===========================================================================
!       This routine rotates a matrix from molecular to crystal coordinates.
!
! The variable eps is a 3x3 output of the subroutine epsilon.
! The variable dmat is a 5x5 matrix rotating d-orbitals.
! The variable pmat is a 3x3 matrix rotating p-orbitals.
!
! Here is the famous Ortega convention:
! In the molecular coordinates we have atoms 1 and 2 (bondcharge) along the
! z-axis; the third atom is in the XZ-plane.
!
! The labelling of the orbitals is as follows:
!
!   S-shell :                s
!                            1
!
!   P-shell :           py   pz   px
!                       1    2    3
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                  1     2   3     4      5
!
! ===========================================================================
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
        subroutine rotate_fb (in1, in2, eps, mmatrix, xmatrix)
        use dimensions
        use interactions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: in1
        integer, intent(in) :: in2
 
        real, intent(in) :: eps (3, 3)
        real, intent(in) :: mmatrix (numorb_max, numorb_max)
 
! Output
        real, intent(out) :: xmatrix (numorb_max, numorb_max)
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer issh
        integer jssh
        integer k1, k2
        integer n1, l1, m1
        integer n2, l2, m2
 
        real dmat (5, 5)
        real left (5, 5)
        real pmat (3, 3)
        real right (5, 5)
 
! Procedure
! ===========================================================================
! Set up the matrices to rotate p and d orbitals.
        call twister (eps, dmat, pmat)
 
        xmatrix=0.0d0
        n1 = 0
        do issh = 1, nssh(in1)
         l1 = lssh(issh,in1)
         call chooser (l1, dmat, pmat, left)
 
         n2 = 0
         do jssh = 1, nssh(in2)
          l2 = lssh(jssh,in2)
          call chooser (l2, dmat, pmat, right)
 
! Compute
          do m2 = 1, 2*l2 + 1
           do m1 = 1, 2*l1 + 1
            do k2 = 1, 2*l2 + 1
             do k1 = 1, 2*l1 + 1
              xmatrix(n1+k1,n2+k2) = xmatrix(n1+k1,n2+k2)   &
     &          + left(k1,m1)*mmatrix(n1+m1,n2+m2)*right(k2,m2)
             end do
            end do
           end do
          end do
 
          n2 = n2 + 2*l2 + 1
         end do
         n1 = n1 + 2*l1 + 1
        end do
 
! Format Statements
! ===========================================================================
 
        return
        end
 
