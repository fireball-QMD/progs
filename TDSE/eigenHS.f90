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


! eigenHS.f90
! Program Description
! ===========================================================================
!       This routine gets noSCF eigenvalues & eigenvectors of H,S
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine eigenHS (itime_step)

        use options
        use tdse

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itime_step

! Local Parameters and Data Declaration
! ===========================================================================


! Local Variable Declaration and Description
! ===========================================================================
        integer iclock
        logical isH

! Procedure
! ===========================================================================
! do we dump psi in this time step?
        iclock = mod(itime_step,np2es)
        write (*,*) 'iclk',itime_step,np2es,iclock
        if ((iclock .eq. 0) .or. (itime_step .eq. 1)) then
          isH = .true.
          write (*,*) ' In this time step we project actual psi onto eigenstate;'
          write (*,*) ' full diagonalization will be performed.'
        else
          isH = .false.
        endif

! assemble H,S
        if (isH) then

! assemble Hamiltonian
         call assemble_h ()
! diagonalize H,S; get S^(1/2),S^(-1/2)
         write (*,*)  ' Call tddiag_k'
         call tddiag_k (isH)

! project psi onto set of eigenstates
         write (*,*) 'call psi2es'
         call get_psi2es (itime_step)

! get charges and rho (??)
!         call build_rho (0)

        else

! assemble S
         call assemble_h ()
! diagonalize S; get S^(1/2),S^(-1/2) and orthogonalize H
         write (*,*)  ' Call tddiag_k'
         call tddiag_k (isH)

        endif

! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine eigenHS

