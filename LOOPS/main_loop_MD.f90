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
        use MD
        use forces
        use constants_fireball
        use energy
        use optimization
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
         call initcDFT ()
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
         itime_step_g = itime_step
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

