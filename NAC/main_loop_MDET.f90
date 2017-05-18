! copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
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


! main_loop_MDET.f90
! Program Description
! ===========================================================================
!       This routine performs MD loop for MDET
!  (Molecular Dynamics with Electronic Transitions)
!
! ===========================================================================

! JOM-info : adapted
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine main_loop_MDET ()

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
        use nonadiabatic
        use density
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer itime_step
        integer iele  
        integer ix    
        integer iatom 
! Procedure
! ===========================================================================
! First initialize:
! We need to define the initial electronic state, as well as the initial
! excited state (for nuclear motion).

! get scf solution of initial state : call scf_loop_na
!       call scf_loop (0)

! get energy of ground state
!       call getenergy (0)

! Initial state for nuclear motion: define foccupy_na and ioccupy_na
! also:
! Initial electronic states: define wfs-coefficients cna for TD-wfs
        write(*,*)'call init_mdet'
        call init_mdet (natoms)
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

! Open the files with the things we want to store (ENRIQUE-JOM):
        open (unit = 210, file = 'dks.dat', status = 'unknown')
        open (unit = 211, file = 'vatom.dat', status = 'unknown')
        open (unit = 212, file = 'energies.dat', status = 'unknown')
        open (unit = 213, file = 'occupancy_MD.dat', status = 'unknown')
        open (unit = 214, file = 'c_na.dat', status = 'unknown')
        open (unit = 215, file = 'dV.dat', status = 'unknown')
	open (unit = 216, file = 'en.dat', status = 'unknown')
!       open (unit = 2121, file = 'random.dat', status = 'old')

        do itime_step = nstepi, nstepf
         itime_step_g = itime_step
 	 write (210,*) " MD step = " , itime_step
	 write (211,*) " MD step = " , itime_step
	 write (212,*) " MD step = " , itime_step
	 write (214,*) " MD step = " , itime_step
	 write (215,*) " MD step = " , itime_step
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

! perform SCF LOOP
! JOM-info : in denmat (now mdetdenmat) we use foccupy_na and ioccupy_na
!            instead of foccupy and ioccupy_k
        write(*,*)'call scf_loop'
        call scf_loop (itime_step)

! DEBUG vlada begin
!         if (iwrtewf .eq. 1) then
!          if (mod(itime_step,10) .eq. 0) then
!            write (*,*) ' Call ewf2mesh subroutine '            
!              call ew2mesh (icluster,itime_step)
!          end if   
!         end if
! DEBUG vlada end

! JOM-info : if (itheory .ne. 0) then we will have two different sets of
! eigenvectors and eigenvalues, one for the nuclear motion and the other
! for the adiabatic KS-states used in the time evolution of electrons.
! the nuclear motion depends on foccupy_na, while the adiabatic
! KS-states will be defined from the Ground State (GS) 
! (or Transition State (TS)) Hamiltonian

! optionally perform post-processing (DOS etc.)
!        call postscf ()

! calculate the total energy
         write(*,*)'call getenergy'
         call getenergy (itime_step)

! ============================================================================
! ----------------------------------------------------------------------------
!                                F O R C E S
! ----------------------------------------------------------------------------
        if (iforce .eq. 1) then
         write (*,*) '  '
         write (*,*) ' =========== Move on to Forces ================= '

! Check if there is not a global sign between old and new wfunction
         write(*,*)'call overlap_sign'
	 call overlap_sign(itime_step)

! Assemble forces and calculate nonadiabatic vector gks
         write(*,*)'call getforces_mdet'
         call getforces_mdet ()

! Calculate numerically the non-adiabatic contribution from the time evolution of the kohn-sham states
! Numerical non-adiabatic contribution is not used  for time evolution or calculating of probabilities.
         write(*,*)'call delta_t_ks'
	 call delta_t_ks (itime_step)

! Time evolution of Kohn-Sham states
         write(*,*)'call evolve_ks_states'
         call evolve_ks_states (itime_step)

! Fewest switches mechanism
         write(*,*)'call fewest_switches'
         call fewest_switches (itime_step)

! Write out occupancy
    ! write(213,'(i4,<nele>f6.1)') itime_step, (2*foccupy_na(map_ks(iele),1),iele=1,nele)
    write (213,'(i4)',advance='no') itime_step
    do iele=1,nele
       write(213,'(f6.1)',advance='no') 2*foccupy_na(map_ks(iele),1) 
    end do
    write (213)

! ----------------------------------------------------------------------------
! JOM-  we need to save eigen_old, psi_old,
! ratom_old, vatom_old...
         write(*,*)'call save_mdstuff'
         call save_mdetstuff
! ----------------------------------------------------------------------------

! Move ions now
! JOM-info : I think it is best to use gear_order = 2 (velocity verlet)
! In this way, the atomic positions are not "corrected" and the
! positions at which the Hamiltonian is solved correspond directly to
! the positions stored in ratom
         write(*,*)'call move_ions "CORRECTOR" part'
         call move_correc (itime_step)

! ----------------------------------------------------------------------------
 
         write(*,*)'call move_ions "PREDICTOR" part'
         call move_predic (itime_step)

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
        end subroutine main_loop_MDET

