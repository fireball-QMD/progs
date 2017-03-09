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


! doscentrosPP.f90
! Program Description
! ======================================================================
!       This routine computes the non-local pseudopotential matrix element
! ppx(mu,nu)
!-----------------------------------------------------------------
! IMPORTANT : the range of neighbours could/should be different for
! 3-center pseudopotential matrix elements than for neutral-coulomb
! contributions.  This could/should be taken into account in the future
!
! Program Declaration
! ===========================================================================
        subroutine doscentrosPP (interaction, isub, distance, eps, deps,     &
     &                           iforce, in1, in2, sx, spx)
        use dimensions
        use interactions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: iforce
        integer, intent(in) :: in1
        integer, intent(in) :: in2
        integer, intent(in) :: isub
        integer, intent(in) :: interaction
        real, intent (in) :: distance

! Output
        real, intent (in), dimension (3, 3, 3) :: deps
        real, intent (in), dimension (3, 3) :: eps
 
        real, intent(out), dimension (numorb_max, numorb_max) :: sx
        real, intent(out), dimension (3, numorb_max, numorb_max) :: spx
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer inu
        integer index

        real, dimension (ME2cPP_max) :: dpplist
        real, dimension (3) :: eta
        real, dimension (ME2cPP_max) :: pplist
        real, dimension (numorb_max, numorb_max) :: sm
        real, dimension (numorb_max, numorb_max) :: spm
        real, dimension (3, numorb_max, numorb_max) :: spmx

! Procedure
! ===========================================================================
! For the non-local Pseudopotential:  interaction=5, isub=0.
        if (interaction .ne. 5 .or. isub .ne. 0) then
         write (*,*) ' interaction = ', interaction
         write (*,*) ' This routine is only for the pseudopotential '
         write (*,*) ' interactions.  The wrong subroutine is being called. '
         write (*,*) ' We must stop here! '
         stop
        end if
 
        sm = 0.0d0
        sx = 0.0d0
        if (iforce .eq. 1) spx = 0.0d0

        do index = 1, index_maxPP(in1,in2)
         call interpolate_1d (interaction, isub, in1, in2, index, iforce,    &
     &                        distance, pplist(index), dpplist(index))
        end do

! When you change from < 1 | 2 > to < 2 | 1 > sometimes you need to change the 
! sign; it is the same thing with the overlap: < s | pz > = - < pz | s > ....
!
! This is done by defining eps2 with 180 degrees rotation (i.e. instead of 
! r3-r2, we use r2-r3 when defining eps2)
        call recover_PP (in1, in2, pplist, sm)
        call recover_PP (in1, in2, dpplist, spm)
 
! Now, rotate into crystal-coordinates: ppm --> ppx
        call rotatePP (in1, in2, eps, sm, sx)
 
! Now for derivatives. Only compute derivative if iforce=1.
! Here is what you should do
! Multiply spm by (-1.0d0)*eta in order to get it into vector form
!                     ^_______dR12/dr1 gives -eta
! The last column of eps(:,3) is eta
        if (iforce .eq. 1) then
         eta(:) = eps(:,3)
         do imu = 1, num_orb(in1)
          do inu = 1, num_orbpp(in2)
           spmx(:,imu,inu) = - eta(:)*spm(imu,inu)
          end do
         end do
 
         call rotatedPP (in1, in2, eps, deps, sm, spmx, spx)
        end if
 
! Format Statements
! ===========================================================================
 
        return
        end
