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


! main_loop_NEB.f90
! Program Description
! ===========================================================================
!       This routine performs Nudged Elastic Band (NEB) loop
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine main_loop_NEB ()
 
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

! perform SCF LOOP
         call scf_loop (itime_step)

! optionally perform post-processing (DOS etc.)
         call postscf () 

! calculate the total energy
         call getenergy (itime_step)


!         write (*,*) '  '
!         write (*,*) ' ================= Move on to Forces ================= '
 
! Assemble forces 
         call getforces ()
         
! doing NEB
         call move_neb (itime_step, etot, ftot)

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
600     format (2x, ' average temperature       = ', f8.2)
601     format (2x, '                    nuclear kinetic energy = ', f18.8)
602     format (2x, ' Grand Total = nuclear kinetic + potential = ', f18.8)
603     format (2x, '                     total energy per atom = ', f18.8)
604     format (2x, '               grand total energy per atom = ', f18.8)
605     format (2x, '                        deltaE/atom  (meV) = ', f18.8)
700     format (2x, i2, 3(2x,f12.6))
701     format (2x, a2, 3(2x,f12.6))
1889    format (2x, ' deltae, etotnew, etotold = ', 3f15.8)
1888    format (2x, ' ACVE ', 3f15.8)
2000    format (3x, ' Checking MD convergence : ')
2001    format (2x,' ++++ iter = ',i8,' Etot= ',f16.8,' Fi_max= ',f14.6)
2010    format (2x,'  +++ Etot  RES =', f16.8,'  TOL = ',f16.8,'   CONVERGED ')
2011    format (2x,'  +++ Etot  RES =', f16.8,'  TOL = ',f16.8,'   NOT CONVERGED ')
2020    format (2x,'  +++ Fmax  RES =', f16.8,'  TOL = ',f16.8,'   CONVERGED ')
2021    format (2x,'  +++ Fmax  RES =', f16.8,'  TOL = ',f16.8,'   NOT CONVERGED ')
        return
        end subroutine main_loop_NEB
 
