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


! postete.f90
! Program Description
! ===========================================================================
!      This subroutine assemble density matrix, energy & forces after
! time-dependent evoution of wavefunction with fixed ionic positions
! ===========================================================================
! Code rewritten by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine postete (itime_step)
        use energy
        use tdse


        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
       integer, intent (in) :: itime_step

! Local Parameters and Data Declaration
! ===========================================================================

       integer  ielec
       integer  imu
       integer  inu
       integer  ikpoint

! Local Variable Declaration and Description
! ===========================================================================


! Procedure
! ===========================================================================

! assemble the density matrix & band structure energy
        call tddenmat ()

! calculate energy
        write (*,*) ' Assemble total energy'
        call getenergy (itime_step)

! calc forces
        call getforces (itime_step)

! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine postete

