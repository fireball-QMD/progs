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

! main_loop_importrho.f90
! Program Description
! ===========================================================================
!       For reading in rho and projecting density to a grid
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine main_loop_importrho ()

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
        integer iatom
        integer imu
        integer in1
        integer iphase
        character (len = 1) reply

        real vscale
        integer, parameter :: phaserho(8) = (/0,45,90,135,180,225,270,315/)
        real etrans

! Procedure
! ===========================================================================
! ===========================================================================
! ---------------------------------------------------------------------------
!                 M O L E C U L A R   D Y N A M I C S   L O O P
! ---------------------------------------------------------------------------
! ===========================================================================
! Begin molecular dynamics simulation.
        write (*,*) '  '
        write (*,*) ' ======================================================= '
        write (*,*) '          Import rho for density grid projection'
        write (*,*) ' ======================================================= '
        write (*,*) '  '
! ===========================================================================
!                          Start loop over time steps
! ===========================================================================
        do itime_step = nstepi, nstepi ! only one time step
         if (iimage .ge. 1)                                                  &
     &    call imaged (icluster, iimage, itime_step, nstepi)

! ***************************************************************************

! perform SCF LOOP
         call scf_loop (itime_step)

! optionally perform post-processing (DOS etc.)
         reply = 'y'
         do while (reply == 'y')
           write(*,*) 'Input energy of transition to import (eV) '
		   read(*,*) etrans
		   do iphase = 1,16
		    call den2mesh_import (icluster,etrans,iphase)
		   end do
           write(*,*)'Read in another energy for grid projection? (y/n) '
           read(*,*) reply
         end do



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
        end subroutine main_loop_importrho

