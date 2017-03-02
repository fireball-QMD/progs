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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

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

