! copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Institute of Physics, Czech Republic - Pavel Jelinek

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

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


! main_loop.f90
! Program Description
! ===========================================================================
!       This routine forks to given time loops: MD, CG, DM, NEB etc.
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine main_loop ()

        use options
        use outputs
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================


! Procedure
! ===========================================================================
! Imported rho for grid projection
	if (idensimport .eq. 1 ) then
	  call main_loop_importrho ()
      return
     endif

! Do TDSE
	if (itdse .eq. 1 ) then
	  call main_loop_TDSE ()
      return
     endif

! Do MDET
	if (imdet .eq. 1 ) then
	  call main_loop_MDET ()
      return
     endif

! Do NAC
	if (imdet .eq. 2 ) then
	  call main_loop_NAC ()
      return
     endif

! Dynamical Matrix loop
     if (idynmat .eq. 1) then
      call main_loop_DM ()
      return
      endif

! Nudged Elastic Band loop
     if (ineb .eq. 1) then
      call main_loop_NEB ()
      return
     endif


! Molecular Dynamics loop
     if ( ( iquench .le. 0 ) .and. ( iquench .gt. -4 ) ) then
      call main_loop_MD ()
      return
     endif

! Conjugated Gradient loop
     if (iquench .eq. -4) then
      call main_loop_CG ()
      return
     endif

! bfgs minimization routine OR  steepest descent minimization routine
     if (iquench .eq. -5) then
      call main_loop_MIN ()
      return
     endif

! FIRE minimization loop
     if (iquench .eq. -6 ) then
      call main_loop_FIRE ()
      return
     endif




! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine main_loop

