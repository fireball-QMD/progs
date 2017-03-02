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

! readdata.f90
! Program Description
! ===========================================================================
!       This routine build up density matrix and evaluate a new charges for
! next SCF-step
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine build_rho (itime_step)

        use scf
        use density
        use outputs
        use options
        use energy
        use configuration
        use interactions

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itime_step

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================


! Procedure
! ===========================================================================


! ===========================================================================
!               compute the density matrices
! ===========================================================================
! Compute the density matrices. The results rho and cape are computed.
! The bandstructure energy is also computed.
! First initialize the density matrices to zero.
          rho = 0.0d0
          cape = 0.0d0
! JOM-MDET
          if (tempfe .le. 50.0d0) tempfe = 50.0d0 ! Can't be zero
! JOM : why not ?
         if (imdet .eq. 1) then 
          call mdetdenmat (ifixcharge, iqout, icluster, iwrtefermi,     &
     &                  tempfe, ebs, iwrtpop)
         else
          call denmat (ifixcharge, iqout, icluster, iwrtefermi, tempfe, ebs, &
     &                  iwrtpop)

         end if
! jel-grid
          if (iwrtewf .eq. 1) then
            write (*,*) ' Call ewf2mesh subroutine. '
            call ew2mesh (icluster)
          endif
! end jel-grid


         call writeout_charges (natoms, ifixcharge, iqout, iwrtcharges,     &
     &                           iwrtdensity, basisfile, symbol)

! ===========================================================================
!                  check input and output charges for scf
! ===========================================================================
! call mixer
         call mixer (natoms, itheory, ifixcharge, iwrtcharges)

! Check convergence of charge; sigmatol is in scf.optional
         if (sigma .lt. sigmatol) scf_achieved = .true.

! If ifixcharge = 1 then do not iterate to self=consistancy...
         if (ifixcharge .eq. 1) then
           write (*,*) '  '
           write (*,*) ' !!! !!! !!! BEWARE !!! !!! !!!'
           write (*,*) ' You have chosen ifixcharge = 1 so we do NOT iterate '
           write (*,*) ' to self consistancy!!! (Normally ifixcharge = 0) '
           scf_achieved = .true.
         end if
         if (.not. scf_achieved .and. itheory .ne. 0)                       &
     &     write (*,*) '            BAD NEWS; results are not self-consistent '

! ===========================================================================
!                            write out cdcoeffs
! ===========================================================================
         if (iwrtcdcoefs .gt. 0 .and. scf_achieved)                        &
          call writeout_cd (icluster, iwrtcdcoefs, itime_step)

! jel-TDSE
! JOM-MDET
! save SCF-wavefunctions for TDSE or MDET run
         if (itdse .eq. 1 or imdet .eq. 1) then
          write (*,*) ' save eigenstuff '
!
         else
! Deallocate arrays
          if (iqout .ne. 2 .and. icluster .ne. 1) deallocate (blowim)
          if (iqout .ne. 2) deallocate (blowre)
          deallocate (eigen_k)
! Deallocate arrays
          if (iordern .ne. 1) then
           if (icluster .ne. 1) deallocate (bbnkim)
           deallocate (bbnkre)
          end if
         endif ! if (tdse, mdet)

! Write out the components of the Hamiltonian if iwrtcomponents = 1
         if(iwrtcomponents .eq. 1) call writeout_comph (natoms, itheory, ebs)

! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))

        return

        end subroutine build_rho

