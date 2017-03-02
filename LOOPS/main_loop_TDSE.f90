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

! main_loop_MD.f90
! Program Description
! ===========================================================================
!       This routine performs main TDSE loop
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine main_loop_TDSE ()

        use options
        use configuration
        use options
        use outputs
        use MD
        use neb
        use forces
        use fragments
        use constants_fireball
        use energy
        use barrier
        use optimization
        use interactions
        use charges


        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer itime_step
        logical itrue

! Procedure
! ===========================================================================

! get scf solution of initial state
        call scf_loop (0)

! get energy of ground state
        call getenergy (0)

! get S, S^(-1/2) matrices
        itrue = .true.
        call tddiag_k (itrue)

! possible excitation
        call initpsi ()


! Begin semiclassical molecular dynamics
! ===========================================================================
! ---------------------------------------------------------------------------
!                 M O L E C U L A R   D Y N A M I C S   L O O P
! ---------------------------------------------------------------------------
! ===========================================================================
! Begin molecular dynamics simulation.
        write (*,*) '  '
        write (*,*) ' ======================================================= '
        write (*,*) '             Begin molecular dynamics simulation. '
        write (*,*) ' ======================================================= '
        write (*,*) '  '
! ===========================================================================
!                          Start loop over ionic time steps
! ===========================================================================
        do itime_step = nstepi, nstepf

! ***************************************************************************

! Calculate Eigenvalues & Eigenvectors for given ionic configuration
         call eigenHS (itime_step)

! Set boundary condition for the edge
!         if (itime_step .ne. 1) call set_tdbc ()

! perform electron time evolution loop
         call ete_loop (itime_step)

! calculate density matrix, energy & forces
         call postete (itime_step)

         write (*,*) '  '
         write (*,*) ' ================= Moving Ions ================= '
! Move ions now
         call move_ions (itime_step)

! ===========================================================================
! ---------------------------------------------------------------------------
!                          E N D   T I M E   L O O P
! ---------------------------------------------------------------------------
       end do

! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))

        return
        end subroutine main_loop_TDSE

