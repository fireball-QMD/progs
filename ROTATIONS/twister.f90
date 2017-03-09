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


! twister.f90
! Program Description
! ===========================================================================
!       This routine prepares the D matrices for a given geometry of the
! matrix element. The set of vectors stored in eps(3,3) are needed.
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
!                            0
!
!   P-shell :           py   pz   px
!                       -1   0    +1
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                 -2    -1   0    +1     +2
! ===========================================================================
! Original code written by Alex Demkov.

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
        subroutine twister (eps, dmat, pmat)
        use constants_fireball
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        real, intent(in) :: eps (3, 3)
 
! Output
        real, intent(out) :: dmat (5, 5)
        real, intent(out) :: pmat (3, 3)
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer jx
        integer ix
        real amat_term
        real xlambda11
        real xlambda12
        real xlambda13
        real xlambda32
        real xlambda33

 
! Procedure
! ===========================================================================
! Nothing for S orbitals
! Set the P matrices: eps - x,y,z => pmat - y,z,x
        pmat(1,1) = eps(2,2) 
        pmat(1,2) = eps(2,3) 
        pmat(1,3) = eps(2,1)
  
        pmat(2,1) = eps(3,2) 
        pmat(2,2) = eps(3,3)
        pmat(2,3) = eps(3,1)
         
        pmat(3,1) = eps(1,2)
        pmat(3,2) = eps(1,3)
        pmat(3,3) = eps(1,1)


        if (.not. haveDorbitals) return

! Set the lambda matrices:
        do imu = 1, 5
         xlambda11 = 0
         xlambda12 = 0
         xlambda13 = 0
         xlambda32 = 0
         xlambda33 = 0
         do jx = 1, 3
          do ix = 1, 3
           if (amat(ix,jx,imu) .ne. 0.0d0) then
            amat_term = amat(ix,jx,imu)
            xlambda11 = xlambda11 + amat_term*eps(jx,1)*eps(ix,1)
            xlambda12 = xlambda12 + amat_term*eps(jx,1)*eps(ix,2)
            xlambda13 = xlambda13 + amat_term*eps(jx,1)*eps(ix,3)
            xlambda32 = xlambda32 + amat_term*eps(jx,3)*eps(ix,2)
            xlambda33 = xlambda33 + amat_term*eps(jx,3)*eps(ix,3)
           end if
          end do
         end do
! Set the D matrices:
         dmat(imu,1) = 2.0d0*xlambda12
         dmat(imu,2) = 2.0d0*xlambda32
         dmat(imu,3) = sqrt(3.0d0)*xlambda33
         dmat(imu,4) = 2.0d0*xlambda13
         dmat(imu,5) = 2.0d0*xlambda11 + xlambda33
        end do
 
! Format Statements
! ===========================================================================
 
        return
        end
 
