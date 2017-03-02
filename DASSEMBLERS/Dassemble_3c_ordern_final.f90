! copyright info:
!
!                             @Copyright 2002
!                            Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega

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
 
! Dassemble_3c_ordern_final.f90
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
        subroutine Dassemble_3c_ordern_final (natoms)
        use forces
        use mpi_declarations
        use neighbor_map
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

        real, dimension (:, :), allocatable :: r2_tmp

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Allocate Arrays
! ===========================================================================
        allocate (r2_tmp (3, natoms))

! Procedure
! ===========================================================================
        call MPI_ALLREDUCE (f3naa, r2_tmp, 3*natoms,                         &
     &                      mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        f3naa = r2_tmp
        call MPI_ALLREDUCE (f3nab, r2_tmp, 3*natoms,                         &
     &                      mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        f3nab = r2_tmp
        call MPI_ALLREDUCE (f3nac, r2_tmp, 3*natoms,                         &
     &                      mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        f3nac = r2_tmp
        call MPI_ALLREDUCE (f3nla, r2_tmp, 3*natoms,                         &
     &                      mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        f3nla = r2_tmp
        call MPI_ALLREDUCE (f3nlb, r2_tmp, 3*natoms,                         &
     &                      mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        f3nlb = r2_tmp
        call MPI_ALLREDUCE (f3nlc, r2_tmp, 3*natoms,                         &
     &                      mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        f3nlc = r2_tmp
        call MPI_ALLREDUCE (f3xca, r2_tmp, 3*natoms,                         &
     &                      mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        f3xca = r2_tmp
        call MPI_ALLREDUCE (f3xcb, r2_tmp, 3*natoms,                         &
     &                      mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        f3xcb = r2_tmp
        call MPI_ALLREDUCE (f3xcc, r2_tmp, 3*natoms,                         &
     &                      mpi_whatever_real, MPI_SUM, MPI_BTN_WORLD, ierror)
        f3xcc = r2_tmp

! Deallocate Arrays
! ===========================================================================
        deallocate (r2_tmp)
 
! Format Statements
! ===========================================================================
 
        return
        end
