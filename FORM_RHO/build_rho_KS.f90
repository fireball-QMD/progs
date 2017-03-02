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

! build_rho_KS.f90
! Program Description
! ===========================================================================
!       This routine build up density matrix and evaluate a new charges for
! next SCF-step for Kohn-Sham theory
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine build_rho_KS (itime_step)
 
        use scf 
        use density
        use outputs
        use options 
        use energy  
        use configuration   
        use interactions
        use kpoints
        
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itime_step

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
	integer imu    
        integer ikpoint
 
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
        if (tempfe .le. 50.0d0) tempfe = 50.0d0 ! Can't be zero

        call denmat_KS (ifixcharge, iqout, icluster, iwrtefermi, tempfe, ebs)



! ===========================================================================
!                  check input and output charges for scf
! ===========================================================================
! call mixer
        call mixer_KS (natoms, ifixcharge)

! calculate output density from the new density matrix
        write (*,*) ' Assemble density'
        call assemble_KS_den (icluster)

! save density matrix
        call writeout_charges_KS (natoms, ifixcharge, iqout, iwrtcharges,     &
     &                           iwrtdensity, basisfile, symbol)



! jel-grid
          if (iwrtewf .eq. 1) then
            write (*,*) ' Call ewf2mesh subroutine. '
            call ew2mesh (icluster)
          end if
! prokop-excitations
         if (iwrtexcit .eq. 1) then
            call excitations
         end if
! prokop kvaziband
         if (iwrtkvaziband .eq. 1) then 
            call kvaziband()
         end if
! prokop band_project
         if (iwrtkvaziband .eq. 2) then 
            	call band_project ()
         end if
! prokop ARPES_LCAO
         if (iwrtkvaziband .eq. 3) then 
            call ARPES_LCAO ()
         end if
! prokop-eigenvectors
         if (iwrtcdcoefs .eq. -1 ) then
			if (allocated ( bbnkim )) then
            	call writeout_eigenvec ()
			else
				call writeout_eigenvec_G ()
			end if
         end if
! prokop-CoefsLCAO
         if (iwrtcdcoefs .eq. -2 ) then
			if (allocated ( bbnkim )) then
            	call writeCoefsLCAO ()
			else
			!	call writeCoefsLCAO_G ()
			end if
         end if


! be carreful to do that with the rho_out (not new rho_in)
!        write (*,*) ' Assemble Double Counting Correction on the Grid'
!        call assemble_KS_dcc (uxcdcc_ks, uhdcc_ks)
!        call assemble_KS_usr (iforce, uiiuee)

! Check convergence of charge; sigmatol is in scf.optional
         if (sigma .lt. sigmatol) then
           scf_achieved = .true.
! cDFT
           if (icDFT .eq. 1 ) then

              write (*,*) 'write out WF-cDFT'
              open  (322, file = 'el-h.dat', status = 'unknown')

! save e and h wf for next projection
              do ikpoint = 1, nkpoints
               write (*,*) 'EIG-h =',eigen_k(id_hole,ikpoint)
               write (*,*) 'EIG-e =',eigen_k(id_elec,ikpoint)
               wf_hole(:,ikpoint) = blowre(:,id_hole,ikpoint)
               wf_elec(:,ikpoint) = blowre(:,id_elec,ikpoint)
               write (323,'(<norbitals>f10.4)')  (eigen_k(imu,ikpoint), imu=1,norbitals) 
! just for a case also into file
               write (322,*)  'ikpts =',ikpoint
               write (322,*)  ((wf_hole(imu, ikpoint)),imu=1, norbitals)
               write (322,*)  ((wf_elec(imu, ikpoint)),imu=1, norbitals)
              end do ! ikpoint
              close (322)

           endif ! endif icDFT
! end cDFT
         endif ! (sigma)

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

! Deallocate arrays
        if (iqout .ne. 2 .and. icluster .ne. 1) deallocate (blowim)
        if (iqout .ne. 2) deallocate (blowre)
        deallocate (eigen_k)


! Deallocate arrays
        if (iordern .ne. 1) then 
          if (icluster .ne. 1) deallocate (bbnkim)
          deallocate (bbnkre)
        end if
 
! Write out the components of the Hamiltonian if iwrtcomponents = 1
        if(iwrtcomponents .eq. 1) call writeout_comph (natoms, itheory, ebs)
        
! Deallocate Arrays
! ===========================================================================

 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))

        return
        
        end subroutine build_rho_KS
 
