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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
! assemble_2c_ordern_final.f90
! Program Description
! ===========================================================================
!       Take contributions from each processor and perform an all_reduce 
! command to get all contributions combined to each processor.
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
        subroutine assemble_2c_ordern_final (natoms)
        use forces
        use interactions
        use neighbor_map
        use mpi_declarations
        implicit none
 
        include 'mpif.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ierror

        real, dimension (:, :, :, :), allocatable :: r4_tmp
        real, dimension (:, :, :, :, :), allocatable :: r5_tmp

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Allocate Arrays
! ===========================================================================
        allocate (r4_tmp (numorb_max, numorb_max, neigh_max, natoms))
        allocate (r5_tmp (3, numorb_max, numorb_max, neigh_max, natoms))

! Procedure
! ===========================================================================
        call MPI_ALLREDUCE (s_mat, r4_tmp,                                   &
     &                      numorb_max*numorb_max*natoms*neigh_max, &
     &                      mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        s_mat = r4_tmp
        call MPI_ALLREDUCE (t_mat, r4_tmp,                                   &
     &                      numorb_max*numorb_max*natoms*neigh_max,          &
     &                      mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        t_mat = r4_tmp
        call MPI_ALLREDUCE (sp_mat, r5_tmp,                                  &
     &                      3*numorb_max*numorb_max*natoms*neigh_max,        &
     &                      mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        sp_mat = r5_tmp
        call MPI_ALLREDUCE (tp_mat, r5_tmp,                                  &
     &                      3*numorb_max*numorb_max*natoms*neigh_max,        &
     &                      mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        tp_mat = r5_tmp

! Deallocate Arrays
! ===========================================================================
        deallocate (r4_tmp, r5_tmp)
 
! Format Statements
! ===========================================================================
 
        return
        end
