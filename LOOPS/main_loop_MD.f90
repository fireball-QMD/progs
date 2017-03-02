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
!       This routine performs MD loop
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine main_loop_MD ()

        use options
        use configuration
        use options
!        use outputs
        use MD
!        use neb
        use forces
!        use fragments
        use constants_fireball
        use energy
!        use barrier
        use optimization
!        use interactions
!        use charges
	use scf

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer itime_step

! Procedure
! ===========================================================================

! in case of constrain DFT, first we perform SCF loop for the ground state;
! This will be used later as a reference state
        if(icDFT .eq. 1 ) then
          call scf_loop (0)
          call getenergy (0)
          cDFT_active = .true.
        endif  ! end icDFT 

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
!                          Start loop over time steps
! ===========================================================================

        do itime_step = nstepi, nstepf
         if (iimage .ge. 1)                                                  &
     &    call imaged (icluster, iimage, itime_step, nstepi)

! ***************************************************************************
! junkermeier: this little section will incrementally increase/decrease
!              the temperature over the course of a calculation.

          if (iendtemp .eq. 1 .and. iquench .eq. 0) then
            T_wantPrev = T_want
            T_want = T_initial + T_increment*itime_step
            call resetNHC(natoms,T_want,T_wantPrev)
          end if
! ***************************************************************************
!CHROM
          if (iclassicMD == 0) then
!END CHROM
! perform SCF LOOP
            call scf_loop (itime_step)

! optionally perform post-processing (DOS etc.)
            call postscf ()

! calculate the total energy
            call getenergy (itime_step)
!CHROM
          endif
!END CHROM
! ============================================================================
! ----------------------------------------------------------------------------
!                                F O R C E S
! ----------------------------------------------------------------------------
         if (iforce .eq. 1) then
          
!CHROM zdenka chromcova - classic 
          if( iclassicMD == 0 )then
!END CHROM 
          write (*,*) '  '
          write (*,*) ' ================= Move on to Forces ================= '
!CHROM
	 endif
!END CHROM 
! Assemble forces

          call getforces ()

! Move ions now
          write (*,*) 'move_ions'
          call move_ions (itime_step)

! ===========================================================================
! ---------------------------------------------------------------------------
!                          E N D          F O R C E S
! ---------------------------------------------------------------------------
         end if ! if(force)



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
        end subroutine main_loop_MD

