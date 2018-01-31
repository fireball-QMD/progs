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


! doscentrosDipY.f90
! Program Description
! ===========================================================================
!      This subroutine calculates the (two-center) matrix elements (mu,nu)
!      for the Y-dipole.
!
!      This routine also calculates the derivative with respect to the
! position of the first atom (the atom of the orbital of the BRA).
!
! ===========================================================================
! Original code written by Jose Ortega.
!
! Code rewritten by:
! Jesús I. Mendieta-Moreno (JIMM)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine doscentrosDipY (interaction, isub, iforce, in1, in2, in3,    & 
     &                         distance, eps, deps, sx, spx)
        use dimensions
        use interactions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: interaction
! Here interaction must be = 10 (Y-dipole)
        integer, intent (in) :: isub
 
! The variable interaction is the type of two-center integral that is needed.
! The variable isub is the subtype of interation (not used here)
!               9            0               z-dipole
!               10           0               y-dipole
!               11           0               x-dipole
 
! Derivatives are computed if and only if iforce = 1.
        integer, intent (in) :: iforce

        integer, intent (in) :: in1
        integer, intent (in) :: in2
        integer, intent (in) :: in3

        real, intent (inout) :: distance
        real, intent (in), dimension (3, 3, 3) :: deps
        real, intent (in), dimension (3, 3) :: eps
 
! Output
! The variable sx is < phi(mu,r-r1) ! vna(r-ratm) ! phi(nu,r-r1)>
! which is a representation in molecular coordinates.
! The variable spx is [d/d(r1) atmna(mu,nu)] where atmna(mu,nu) is
! < phi(mu,r-r1) ! vna(r-ratm) ! phi(nu,r-r1)>.
        real, intent (out), dimension (numorb_max, numorb_max) :: sx
        real, intent (out), dimension (3, numorb_max, numorb_max) :: spx
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer inu
        integer index
 
        real, dimension (3) :: eta
 
! -slist = output list of matrix elements
! -dslist = output list of derivatives of matrix elements
! JIMM: the dimensions here are different for X,Y dipoles
        real, dimension (ME2cDipY_max) :: dslist
        real, dimension (ME2cDipY_max) :: slist

        real, dimension (numorb_max,numorb_max) :: sm
        real, dimension (numorb_max,numorb_max) :: spm
        real, dimension (3,numorb_max,numorb_max) :: spmx

        logical switch

! Procedure
! ===========================================================================
! For the atom case, in3 = in1, but for everything else in3 = in2.
! For the ontop case, in2 = in1 (left) or in3 (right).
! Initialize sm, scam and sx to zero.
! < n1 | n2 | n3 >
        sm = 0.0d0
        sx = 0.0d0
        if (iforce .eq. 1) spm = 0.0d0
        if (iforce .eq. 1) spmx = 0.0d0
        if (iforce .eq. 1) spx = 0.0d0


! This subroutine calls the subroutine intrp1d as needed to find the value of
! the matrix elements for any given atomic separation distance.
! -slist = output list of matrix elements
! -dslist = output list of derivatives of matrix elements
        do index = 1, index_max2cDipY(in1,in3)
          call interpolate_1d (interaction, isub, in1, in2, index, iforce,   &
     &                         distance, slist(index), dslist(index))
        end do
 
! Now recover sm ans spm which are two-dimensional arrays from
! slist and dslist which are one-dimensional arrays.
! JIMM: the recover subroutines are specific for X,Y dipoles.
        call recover_2cDipY (in1, in3, slist, sm)
        call recover_2cDipY (in1, in3, dslist, spm)
 
! Rotate sm into crystal-coordinates: sm --> sx
        call rotate_fb (in1, in3, eps, sm, sx)
 
! ****************************************************************************
!
! FORCES
! ****************************************************************************
! spm   is the "scalar" derivative of the matrix.
!
! When we are done, we get:
! spx   is the vector derivative of the matrix.
!
! Only compute derivative if and only if iforce = 1.
        if (iforce .eq. 1) then
 
! As long as epsilon1 is called with sighat in the second "spot" as
! call epsilon1(R1,sighat,spe), then eps(ix,3) = eta(ix).
         eta(:) = eps(:,3)
 
! First do the matrix.
         do inu = 1, num_orb(in3)
          do imu = 1, num_orb(in1)
 
! Note that if we are calculating the on-site matrix elements, then the
! derivatives should be exactly zero.  This is what Otto referred to as the
! ferbie test.  For example, for the on-site overlap, we get an identity
! matrix, thus surely we should never use a "derivative of the unit matrix".
 
! Note the minus sign. d/dr1 = - eta * d/dd.
           if (distance .gt. 1.0d-3) then
            spmx(:,imu,inu) = - eta(:)*spm(imu,inu)
           end if
          end do
         end do
 
         call rotated (in1, in3, eps, deps, sm, spmx, spx)

        end if

! Format Statements
! ===========================================================================
 
        return
      end subroutine doscentrosDipY
