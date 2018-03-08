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


! initboxes.f90
! Program Description
! ===========================================================================
!       This program sets up the neighboring cells based on the lattice
! vectors input. Originally, this routine was called the boxesP.f routine,
! which came from the original self-consistent code as developed by
! Demkov et al.
! To increase the size of the cube you have to uncomment definitions
! of xxl and msub variables for the desired cube size and corresponding part
! of the loop ix,iy
! Usually, for large systems  the cube 5x5x5 is enough, however for bulk
! calculations you'll need to increase it up to 9x9x9 following the size of
! the unit cell.
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
!
! Modified by P. Jelinek (20/10/2006)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine initboxes (itime)

        use configuration
        use options, only : verbosity
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input:
        integer itime

! Output:

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer lbeta
        integer ix, iy
        integer mbeta
        integer midl
!        real, save, dimension (3, 0:124) :: xxl ! cube 5x5x5
!        real, parameter :: mbox = 2
!        real, save, dimension (3, 0:342) :: xxl ! cube 7x7x7
!        real, parameter :: mbox = 3
        real, save, dimension (3, 0:728) :: xxl ! cube 9x9x9
        integer, parameter :: mbox = 4


! Procedure
! ===========================================================================
! Only set-up xxl on the first time step, otherwise just calculate new xl's.
        if (itime .eq. 1) then

! Initialize 27 "neighboring" cells at lattice vector xl(3,0:124).
! The units are angstrom, and NOT units of atomic lattice.
! The central cell is 0,0,0.
! We set up xxl in units of atomic lattice,
! then put xxl ------> xl in angstom.
         xxl (:,:) = 0.0d0

         xxl (1,1) = 1.0d0

         xxl (1,2) = -1.0d0

         xxl (2,3) = 1.0d0

         xxl (2,4) = -1.0d0

         xxl (3,5) = 1.0d0

         xxl (3,6) = -1.0d0

         xxl (1:2,7) = 1.0d0

         xxl (1,8) = -1.0d0
         xxl (2,8) = 1.0d0

         xxl (1:2,9) = -1.0d0

         xxl (1,10) = 1.0d0
         xxl (2,10) = -1.0d0

         xxl (:,11) = 1.0d0

         xxl (1,12) = -1.0d0
         xxl (2:3,12) = 1.0d0

         xxl (1:2,13) = -1.0d0
         xxl (3,13) = 1.0d0

         xxl (:,14) = 1.0d0
         xxl (2,14) = -1.0d0

         xxl (1:2,15) = 1.0d0
         xxl (3,15) = -1.0d0

         xxl (:,16) = -1.0d0
         xxl (2,16) = 1.0d0

         xxl (:,17) = -1.0d0

         xxl (1,18) = 1.0d0
         xxl (2:3,18) = -1.0d0

         xxl (1,19) = 1.0d0
         xxl (3,19) = 1.0d0

         xxl (1,20) = -1.0d0
         xxl (3,20) = 1.0d0

         xxl (2:3,21) = 1.0d0

         xxl (2,22) = -1.0d0
         xxl (3,22) = 1.0d0

         xxl (1,23) = 1.0d0
         xxl (3,23) = -1.0d0

         xxl (1,24) = -1.0d0
         xxl (3,24) = -1.0d0

         xxl (2,25) = 1.0d0
         xxl (3,25) = -1.0d0

         xxl (2:3,26) = -1.0d0

! Set up the rest of xxl. Set up the bottom layer, the three middle
! layers and then the top layer. The inside 3X3X3 has already been set up
! by the data statement. Set up:
! (i) bottom layer at z=-2, and is 5X5.
! (ii) the three middle layers which enclose the inner 3X3.
! (iii) the top layer at z=2, and is 5X5.
! We have 26 lbeta's so far, so let's start at 26 and increase the number.
         lbeta = 26

! -4 layer for the cube 9x9x9
         do ix = -mbox, mbox
          do iy = -mbox, mbox
           lbeta = lbeta + 1
           xxl(1,lbeta) = real(ix)
           xxl(2,lbeta) = real(iy)
           xxl(3,lbeta) = -4.0d0
          end do
         end do

! -3 layer for the cube 7x7x7
         do ix = -mbox, mbox
          do iy = -mbox, mbox
           lbeta = lbeta + 1
           xxl(1,lbeta) = real(ix)
           xxl(2,lbeta) = real(iy)
           xxl(3,lbeta) = -3.0d0
          end do
         end do

! -2 layer for the cube 5x5x5
         do ix = -mbox, mbox
          do iy = -mbox, mbox
           lbeta = lbeta + 1
           xxl(1,lbeta) = real(ix)
           xxl(2,lbeta) = real(iy)
           xxl(3,lbeta) = -2.0d0
          end do
         end do

! Middle layers
         do midl = -1, 1
          do ix = -mbox, mbox
           do iy = -mbox, mbox

! Skip the "core" 3X3X3 cube, and only include the "skin"
            if (iabs(ix) .gt. 1 .or. iabs(iy) .gt. 1) then
             lbeta = lbeta + 1
             xxl(1,lbeta) = real(ix)
             xxl(2,lbeta) = real(iy)
             xxl(3,lbeta) = real(midl)
            end if
           end do
          end do
         end do

! Now the 2.top layer, the cube 5x5x5
         do ix = -mbox, mbox
          do iy = -mbox, mbox
           lbeta = lbeta + 1
           xxl(1,lbeta) = real(ix)
           xxl(2,lbeta) = real(iy)
           xxl(3,lbeta) = 2.0d0
          end do
         end do

! Now the 3.top layer, the cube 7x7x7
         do ix = -mbox, mbox
          do iy = -mbox, mbox
           lbeta = lbeta + 1
           xxl(1,lbeta) = real(ix)
           xxl(2,lbeta) = real(iy)
           xxl(3,lbeta) = 3.0d0
          end do
         end do

! Now the 4.top layer, the cube 9x9x9
         do ix = -mbox, mbox
          do iy = -mbox, mbox
           lbeta = lbeta + 1
           xxl(1,lbeta) = real(ix)
           xxl(2,lbeta) = real(iy)
           xxl(3,lbeta) = 4.0d0
          end do
         end do

        end if ! if (itime)

         mbeta_max = lbeta
         if (verbosity .ge. 3) write(*,*) 'mbeta_max = ', mbeta_max
! allocate xl
         allocate ( xl(3,0:mbeta_max) )

! put xxl into xl with real angstrom units.
         do mbeta = 0,mbeta_max
          xl(:,mbeta) = xxl(1,mbeta)*a1vec(:) + xxl(2,mbeta)*a2vec(:)   &
     &                + xxl(3,mbeta)*a3vec(:)
         end do

! Format Statements
! ===========================================================================

        return
      end subroutine initboxes
