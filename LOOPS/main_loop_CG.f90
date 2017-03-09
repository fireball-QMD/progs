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


! main_loop_CG.f90
! Program Description
! ===========================================================================
!       This routine performs Conjugated Gradient optimization loop
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine main_loop_CG ()
 
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
        
        real vscale
         
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

!CHROM this part is omitted in case of empirical potentials
          if( iclassicMD == 0 )then
! perform SCF LOOP
         	call scf_loop (itime_step)
! optionally perform post-processing (DOS etc.)
         	call postscf () 

! calculate the total energy
         	call getenergy (itime_step)
         end if
!CHROM end
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
         
! For CG optimization         
          write (*,*) '  Call CG subroutine with status = ',istatus
          call cgo (etotper, ftot, iqout, xvfile, iforce, itheory, iwrtxyz)
          if (istatus .eq. 10) then
            exit
          endif
! addtional optimization via quenching (not allowed now)
          if (istatus .eq. 20) then
           exit
          endif

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
        end subroutine main_loop_CG
 
