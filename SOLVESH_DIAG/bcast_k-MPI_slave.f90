! copyright info:
!
!                             @Copyright 2005
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek
! Brigham Young University - Hao Wang

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
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


! diag_k.f90
! Program Description
! ===========================================================================
!       This is  a subroutine of k-loop
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================

! Broadcast data only once
        subroutine Bcast_data_0_slave (natoms, iqout, Kscf)
          use mpi_declarations
          use kpoints
          use neighbor_map
          use interactions
          use configuration
          use density
          implicit none

          include 'mpif.h'

          integer, intent (out) :: natoms
          integer, intent (out) :: iqout
          integer, intent (out) :: Kscf

          integer ierror
          integer int_mybuffer (7)
          integer nspecies

! Receive nspecies from master
          call MPI_BCAST (nspecies, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

! Allocate data
          allocate (num_orb(nspecies))

! Receive num_orb
          call MPI_BCAST (num_orb, nspecies, MPI_INTEGER, 0, MPI_COMM_WORLD, &
                          ierror)

! Receive data from master
          call MPI_BCAST (int_mybuffer, 7, MPI_INTEGER, 0, MPI_COMM_WORLD, &
                          ierror)
          natoms = int_mybuffer (1) 
          norbitals = int_mybuffer (2)   
          iqout = int_mybuffer (3)  
          nkpoints = int_mybuffer (4)  
          Kscf = int_mybuffer (5)  
          mbeta_max = int_mybuffer (6)  
          numorb_max = int_mybuffer (7)  

! Allocate  arrays 
          allocate (special_k(3,nkpoints))
          allocate (imass(natoms))
          allocate (degelec(natoms))
          allocate (xl(3,0:mbeta_max))
          allocate (ratom(3,natoms))

! Allocate the arrays for wavefunction coefficients
          if (iqout .ne. 2) then
           allocate (blowre (norbitals, norbitals, nkpoints))
           allocate (blowim (norbitals, norbitals, nkpoints))
          endif
          allocate (bbnkre (norbitals, norbitals, nkpoints))
          allocate (bbnkim (norbitals, norbitals, nkpoints))
          allocate (eigen_k (norbitals, nkpoints))
! Allocate neighbors
          allocate (neighn (natoms))
          allocate (neighPPn (natoms))

! Receive data from master
          call MPI_BCAST (imass, natoms, MPI_INTEGER, 0, MPI_COMM_WORLD,   &
                          ierror)
          call MPI_BCAST (degelec, natoms, MPI_INTEGER, 0, MPI_COMM_WORLD, &
                          ierror)
          call MPI_BCAST (xl, 3*(1+mbeta_max), mpi_whatever_real, 0, MPI_COMM_WORLD, &
                          ierror)
          call MPI_BCAST (special_k, 3*nkpoints, mpi_whatever_real, 0,  &
                          MPI_COMM_WORLD, ierror)


          return
        end subroutine Bcast_data_0_slave
!==========================================================================
!==========================================================================

        subroutine Bcast_data_1_slave (natoms)
          use mpi_declarations
          use neighbor_map
          use interactions
          use configuration
          implicit none

          include 'mpif.h'

          integer, intent (in) :: natoms

          integer ierror
          integer int_mybuffer (2)
          integer ndim

! Receive nspecies from master (called only 1. step in SCF)
          call MPI_BCAST (int_mybuffer, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, &
                          ierror)
          neigh_max = int_mybuffer(1)
          neighPP_max = int_mybuffer(2)

          if (allocated (s_mat)) then 
            deallocate (s_mat)
            deallocate (h_mat)
            deallocate (vnl)
            deallocate (neigh_b)
            deallocate (neigh_j)
            deallocate (neighPP_b)
            deallocate (neighPP_j)
          endif
            
          allocate (s_mat (numorb_max, numorb_max, neigh_max, natoms))
          allocate (h_mat (numorb_max, numorb_max, neigh_max, natoms))
          allocate (vnl (numorb_max, numorb_max, neighPP_max**2, natoms))
          allocate (neigh_b (neigh_max, natoms))
          allocate (neigh_j (neigh_max, natoms))
          allocate (neighPP_b (neighPP_max**2, natoms))
          allocate (neighPP_j (neighPP_max**2, natoms))


! Receive data from master
! ratom
          ndim = 3*natoms
          call MPI_BCAST (ratom, ndim, mpi_whatever_real, 0, MPI_COMM_WORLD, &
                          ierror)
! neigh
          ndim = natoms
          call MPI_BCAST (neighn, ndim, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
          ndim = neigh_max*natoms
          call MPI_BCAST (neigh_j, ndim, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
          call MPI_BCAST (neigh_b, ndim, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
! neighPP
          ndim = natoms
          call MPI_BCAST (neighPPn, ndim, MPI_INTEGER, 0, MPI_COMM_WORLD,  &
                         ierror)
          ndim = neighPP_max*neighPP_max*natoms
          call MPI_BCAST (neighPP_j, ndim, MPI_INTEGER, 0, MPI_COMM_WORLD, &
                          ierror)
          call MPI_BCAST (neighPP_b, ndim, MPI_INTEGER, 0, MPI_COMM_WORLD, &
                          ierror)

          ndim = numorb_max*numorb_max*neigh_max*natoms
          call MPI_BCAST (s_mat, ndim, mpi_whatever_real, 0, MPI_COMM_WORLD, &
                         ierror)
          call MPI_BCAST (h_mat, ndim, mpi_whatever_real, 0, MPI_COMM_WORLD, &
                         ierror)
          ndim = numorb_max*numorb_max*neighPP_max*neighPP_max*natoms
          call MPI_BCAST (vnl, ndim, mpi_whatever_real, 0, MPI_COMM_WORLD, &
                         ierror)


          return
        end subroutine Bcast_data_1_slave

       subroutine Bcast_data_2_slave (natoms)
          use mpi_declarations
          use neighbor_map
          use interactions
          implicit none

          include 'mpif.h'

          integer, intent (in) :: natoms

          integer ierror
          integer ndim

          ndim = numorb_max*numorb_max*neigh_max*natoms
          call MPI_BCAST (h_mat, ndim, mpi_whatever_real, 0, MPI_COMM_WORLD, &
                         ierror)
          ndim = numorb_max*numorb_max*neighPP_max*neighPP_max*natoms
          call MPI_BCAST (vnl, ndim, mpi_whatever_real, 0, MPI_COMM_WORLD, &
                         ierror)

          return
        end subroutine Bcast_data_2_slave



