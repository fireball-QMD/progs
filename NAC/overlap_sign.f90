! copyright info:
!
!                             @Copyright 2009
!                FAST (Fireball Atomic Simulation Techniques)
! West Virginia University - James P. Lewis, Chair
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


! overlap_sign.f90
! Program Description
! ===========================================================================
! Calculate non-adiabatic coupling(d_{jk}  contribution V.d_{jk}
! using Kohn-Sham states at different
! time steps, and compares with the equivalent contribution obtained
! directly using the non-adiabtic couplings calculated in
! nacouplings.f90
!
! ===========================================================================
! Code written by Enrique Abad Gonzalez
! Departmento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine overlap_sign (itime_step)

        use configuration
        use nonadiabatic
        use density
        use interactions
        use kpoints
        use MD

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer itime_step


! Local Parameters and Data Declaration
! ===========================================================================
       real, parameter :: hbar = 0.65822d0


! Local Variable Declaration and Description
! ===========================================================================
        integer ikpoint, imu, iorbital


! ===========================================================================
! Calculate non-adiabatic couplings using Kohn-Sham states at different
! time steps
! (1) Define overlap matrix between orbitals mu at time "t" and orbitals
! nu at time "t+dt"
    call overlap_numeric(itime_step)

! Check if this overlap matrix is not unit matrix and change in consequently
    do ikpoint = 1, nkpoints
      do imu = 1, nele
            if ( sumb(imu,imu) .lt. -0.1 ) then
              write (*,*) 'The arbitrary sign of wavefunctions of this and previous time step are different'
              write (*,*) 'We will change the sign in order to have the same sign in all time steps!!'
              do iorbital = 1, norbitals
                bbnkre(iorbital,map_ks(imu),ikpoint) = - bbnkre(iorbital,map_ks(imu),ikpoint)
              end do
            end if
      end do
    end do ! end do ikpoints



! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format ('c_na',2i4,f8.4,2(f7.3))
300     format ('NAC-SUMS',2i4,2f8.4)
301     format ('S(t,tprime)',2i4,1f8.4)
400     format ('S',4f7.3)


        return
        end subroutine overlap_sign

