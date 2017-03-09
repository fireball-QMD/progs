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
!       This routine performs static calculations of NAC d_ij= <Psi_i|d/dR Psi_j>
!
! ===========================================================================

! JOM-info : adapted
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine main_loop_NAC ()

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

! Procedure
! ===========================================================================
! First initialize:
! We need to define the initial electronic state, as well as the initial
! excited state (for nuclear motion).

! Initial state for nuclear motion: define foccupy_na and ioccupy_na
! also:
! Initial electronic states: define wfs-coefficients cna for TD-wfs
        write(*,*)'call init_mdet'
        call init_mdet (natoms)

! Open the files with the things we want to store (ENRIQUE-JOM):
        open (unit = 210, file = 'gks.dat', status = 'unknown')
        open (unit = 211, file = 'vatom.dat', status = 'unknown')
        open (unit = 212, file = 'energies.dat', status = 'unknown')

! Do MC-loop 
        do itime_step = nstepi, nstepf
! get scf solution of initial state : call scf_loop_na
          call scf_loop (itime_step)

! get energy of ground state
          call getenergy (itime_step)

! Check if there is not a global sign between old and new wfunction

! Assemble forces & NAC
          call getnac ()

! enddo MC Loop
          call MCsolar(itime_step)
          
        enddo
        
! close files
        close (210)
        close (211)
        close (212)
        
! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))

        return
        end subroutine main_loop_NAC

