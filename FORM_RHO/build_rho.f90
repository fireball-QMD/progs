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
        use hartree_fock
        use kpoints
        use nonadiabatic
        use dynamo ! e-ph coupling

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
! JOM : why not ?
! JOM-MDET
! jel-nac
!         if (imdet .eq. 1) then
         if (imdet .ne. 0) then
          if (iProjWF .eq. 1) then
              call project_wfmdet ()
           end if
! jel-nac
!          call check_swap (itime_step,Kscf)
          call mdetdenmat (ifixcharge, iqout, icluster, iwrtefermi,     &
     &                  tempfe, ebs, iwrtpop)
         else
!         if (tempfe .le. 50.0d0) tempfe = 50.0d0 ! Can't be zero
          call denmat (ifixcharge, iqout, icluster, iwrtefermi, tempfe, ebs, &
     &                  iwrtpop,bmix,Kscf,igap)

     ! GAP ENRIQUE-FF
           if ((igap .eq. 1).or.(igap .eq. 2)) then
	     call build_ji
           end if
! End GAP ENRIQUE-FF
         end if


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
			    write (*,*) "writing LCAO coefficients for PBC systems"
			    call writeCoefsLCAO ()
			else
			    write (*,*) "writing LCAO coefficients for clusters"
			    call writeCoefsLCAO_G ()
			end if
         end if


         call writeout_charges (natoms, ifixcharge, iqout, iwrtcharges,     &
     &                           iwrtdensity, basisfile, symbol)

! ===========================================================================
!                  check input and output charges for scf
! ===========================================================================
! call mixer
         call mixer (natoms, itheory, ifixcharge, iwrtcharges)

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
               
               !write (323,'(<norbitals>f10.4)')  (eigen_k(imu,ikpoint), imu=1,norbitals) 
               do imu=1,norbitals
                write (323,'(f10.4)', advance='no')  eigen_k(imu,ikpoint)
               end do
! just for a case also into file
               write (322,*)  'ikpts =',ikpoint
               write (322,*)  ((wf_hole(imu, ikpoint)),imu=1, norbitals)
               write (322,*)  ((wf_elec(imu, ikpoint)),imu=1, norbitals)
              end do ! ikpoint
              close (322)

           endif ! endif icDFT
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
         if (iwrtcdcoefs .gt. 0 .and. (scf_achieved .or. (igap .eq. 2)) )   &
          call writeout_cd (icluster, iwrtcdcoefs, itime_step)

! begin VLADA-MDET  
         if ( itime_step .eq. 1 .and. Kscf .eq. 1 ) then
           allocate (bbnkre_o(norbitals,norbitals,nkpoints))
           allocate (blowre_o(norbitals,norbitals,nkpoints))
         end if
         if (scf_achieved) then
           blowre_o(:,:,1)=blowre(:,:,1)
         end if
! end VLADA-MDET
! JOM-MDET
! save SCF-wavefunctions for TDSE or MDET run
         if (itdse .eq. 1 .or. imdet .ne. 0) then
          write (*,*) ' save eigenstuff '
!
         else
! Deallocate arrays
          if (iqout .ne. 2 .and. icluster .ne. 1) deallocate (blowim)
          if (iqout .ne. 2) deallocate (blowre)
! save eigenvalues for e-ph couplig in phimat
          if ((idynmat .eq. 1) .and. (iephc .eq. 1)) then
            write (*,*) ' save eigenstuff for phimat'
! allocate aux array to save eigenvalues from eigen_k
            if(.not. allocated(eigsave)) allocate(eigsave(norbitals, nkpoints))
            eigsave = eigen_k
          endif

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

