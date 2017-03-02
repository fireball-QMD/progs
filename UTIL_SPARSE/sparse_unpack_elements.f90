! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
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

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! sparse_unpack_elements.f90
! Program Description
! ===========================================================================
!
!
! ===========================================================================
! Original code written by Spencer Shellman.

! Code rewritten by:
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
        subroutine sparse_unpack_elements (packarray, packsize, nAmax,       &
     &                                     nprowsmax2, numA, elementsA, packpos)
        use mpi_declarations
        implicit none
 
        include 'mpif.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nAmax
        integer, intent (in) :: nprowsmax2
        integer, intent (in) :: packsize

        integer, intent (in), dimension (nprowsmax2) :: numA

! input byte array; should have size at least
! kind(integer)*(3+isendrows*(1+imax)) + kind(real)*2*isendrows*imax
! where 'imax' = the maximum nonzero element count in a column (this value can
! be obtained by calling sparse_getpacksize)
        integer*1, intent (in), dimension (packsize) :: packarray

! Output
! packpos = on input, the position within the array at which to begin packing;
! on output, number of bytes in the packed array plus initial packsize value.
        integer, intent (inout) :: packpos

        real, intent (out), dimension (nAmax, nprowsmax2) :: elementsA
!$ volatile numA,elementsA
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ierror
        integer imu
        integer index
        integer isendrows
        integer itotal
        integer maxcol

        real, dimension (:), allocatable :: hold_elements

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Initialize the temporary array to zero:
        elementsA = 0

        call MPI_UNPACK (packarray, packsize, packpos, maxcol, 1,            &
     &                   MPI_INTEGER, MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (packarray, packsize, packpos, isendrows, 1,         &
     &                   MPI_INTEGER, MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (packarray, packsize, packpos, itotal, 1,            &
     &                   MPI_INTEGER, MPI_BTN_WORLD, ierror)

        allocate (hold_elements (itotal))

        call MPI_UNPACK (packarray, packsize, packpos, hold_elements,        &
     &                   itotal, mpi_whatever_real, MPI_BTN_WORLD, ierror)
        itotal = 0
        do imu = 1, isendrows
         do index = 1, numA(imu)
          itotal = itotal + 1
          elementsA(index,imu) = hold_elements(itotal)
         end do
        end do

! Deallocate Arrays
! ===========================================================================
        deallocate (hold_elements)
 
! Format Statements
! ===========================================================================
 
        return
        end

