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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

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
        subroutine  Bcast_data_0 (natoms, iqout, Kscf)
          use mpi_declarations
          use interactions
          use configuration
          use kpoints
          implicit none

          include 'mpif.h'

          integer, intent (in) :: natoms
          integer, intent (in) :: iqout
          integer, intent (in) :: Kscf

          integer ierror
          integer mybuffer (7)


          mybuffer (1) = natoms
          mybuffer (2) = norbitals  
          mybuffer (3) = iqout 
          mybuffer (4) = nkpoints 
          mybuffer (5) = Kscf
          mybuffer (6) = mbeta_max
          mybuffer (7) = numorb_max
! integer variables
          call MPI_BCAST (mybuffer, 7, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
! imass
          call MPI_BCAST (imass, natoms, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
! degelec
          call MPI_BCAST (degelec, natoms, MPI_INTEGER, 0, MPI_COMM_WORLD, &
                          ierror)
! xl
         call MPI_BCAST (xl, 3*mbeta_max, mpi_whatever_real, 0, MPI_COMM_WORLD, &
                          ierror)
! special_k
          call MPI_BCAST (special_k, 3*nkpoints, mpi_whatever_real, 0,  &
                          MPI_COMM_WORLD, ierror)

          return
        end subroutine Bcast_data_0
!=============================================================================
!=============================================================================

! Broadcast data only once
        subroutine Bcast_data_1 (natoms)
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

          int_mybuffer(1) = neigh_max
          int_mybuffer(2) = neighPP_max

! Send nspecies from master (called only 1. step in SCF)
          call MPI_BCAST (int_mybuffer, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, &
                          ierror)
!ratom
          ndim = 3*natoms
          call MPI_BCAST (ratom, ndim, mpi_whatever_real, 0, MPI_COMM_WORLD, &
                          ierror)
! neigh
          ndim = natoms
          call MPI_BCAST (neighn, ndim, MPI_INTEGER, 0, MPI_COMM_WORLD, &
                          ierror)
          ndim = neigh_max*natoms
          call MPI_BCAST (neigh_j, ndim, MPI_INTEGER, 0, MPI_COMM_WORLD, &
                          ierror)
          call MPI_BCAST (neigh_b, ndim, MPI_INTEGER, 0, MPI_COMM_WORLD, &
                          ierror)
! neighPP
          ndim = natoms
          call MPI_BCAST (neighPPn, ndim, MPI_INTEGER, 0, MPI_COMM_WORLD, &
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
        end subroutine Bcast_data_1

!============================================================================
!=============================================================================

      subroutine Bcast_data_2 (natoms)
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
        end subroutine Bcast_data_2

!============================================================================
!=============================================================================

       subroutine Gather_k (natoms, iqout, norbitals)
          use mpi_declarations
          use density
          use kpoints
          implicit none

          include 'mpif.h'
! Argument Declaration and Description
! ===========================================================================
! Input
          integer, intent (in) :: natoms
          integer, intent (in) :: iqout
          integer, intent (in) :: norbitals 

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
          integer ierror
          real, dimension (:, :), allocatable :: r2_tmp
          real, dimension (:, :, :), allocatable :: r3_tmp

! Allocate Arrays
! ===========================================================================
          allocate (r2_tmp (norbitals, nkpoints))
          allocate (r3_tmp (norbitals, norbitals, nkpoints))

! eigenvalue
          call MPI_ALLREDUCE (eigen_k, r2_tmp, norbitals*nkpoints,            &
     &                      mpi_whatever_real, MPI_SUM, MPI_COMM_WORLD, ierror)
          eigen_k = r2_tmp
! eigenvectors AO
          call MPI_ALLREDUCE (bbnkre, r3_tmp, norbitals*norbitals*nkpoints,  &
     &                      mpi_whatever_real, MPI_SUM, MPI_COMM_WORLD, ierror)
          bbnkre = r3_tmp
          call MPI_ALLREDUCE (bbnkim, r3_tmp, norbitals*norbitals*nkpoints,  &
     &                      mpi_whatever_real, MPI_SUM, MPI_COMM_WORLD, ierror)
          bbnkim = r3_tmp
! eigenvectors MO (Lowdin)
          if (iqout .ne. 2) then
           call MPI_ALLREDUCE (blowre, r3_tmp, norbitals*norbitals*nkpoints,  &
     &                       mpi_whatever_real, MPI_SUM, MPI_COMM_WORLD, ierror)
           blowre = r3_tmp
           call MPI_ALLREDUCE (blowim, r3_tmp, norbitals*norbitals*nkpoints,  &
     &                       mpi_whatever_real, MPI_SUM, MPI_COMM_WORLD, ierror)
           blowim = r3_tmp
          endif

! Deallocate Arrays
! ===========================================================================
          deallocate (r2_tmp, r3_tmp)

! Format Statements
! ===========================================================================

        return
        end

!============================================================================
!=============================================================================
        subroutine Bcast_np (nspecies)
          use interactions
          implicit none
          include 'mpif.h'
          integer, intent (in) :: nspecies

          integer ierror 
! send nspecies
          call MPI_BCAST (nspecies, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
! send num_orb(:)
          call MPI_BCAST (num_orb, nspecies, MPI_INTEGER, 0, MPI_COMM_WORLD, &
                          ierror)

          return          
        end subroutine bcast_np
