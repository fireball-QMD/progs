! copyright info:
!
!                             @Copyright 2006
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

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

 
! assemble_2c_ordern_init.f90
! Program Description
! ===========================================================================
!       This routine intitializes some things by communicating to all 
! processors information read from the master processor.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_2c_ordern_init (natoms, nspecies, ivdw)
        use dimensions
        use forces
        use interactions
        use neighbor_map
        use integrals
        use constants_fireball
        use mpi_declarations
        implicit none
 
        include 'mpif.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms
        integer, intent (in) :: nspecies
        integer, intent (in) :: ivdw
 
! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ierror
        integer packsize
        integer packpos
 
! Broadcast buffers & sizes
        integer kbsize
        integer krbsize
        integer*1, dimension (:), allocatable :: kbcb
        integer*1, dimension (:), allocatable :: krbcb

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        call bcast (neigh_self, natoms, MPI_INTEGER)
        if (ivdw .eq. 1) call bcast (neigh_vdw_self, natoms, MPI_INTEGER)

        call bcast (neigh_back, natoms*neigh_max, MPI_INTEGER)
        call bcast (sVNL, numorb_max*numorb_max*natoms*neigh_max,            &
     &              mpi_whatever_real)
        call bcast (spVNL, 3*numorb_max*numorb_max*natoms*neigh_max,         &
     &              mpi_whatever_real)
        call bcast (vxc_1c, numorb_max*numorb_max*natoms*neigh_max,          &
     &              mpi_whatever_real)
        call bcast (xcnu1c, nsh_max*nsh_max*nspecies, mpi_whatever_real)
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
