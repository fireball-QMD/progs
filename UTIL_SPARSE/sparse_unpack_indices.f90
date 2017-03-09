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


! sparse_unpack_indices.f90
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
        subroutine sparse_unpack_indices (packarray, packsize, nAmax,        &
     &                                    nprowsmax, numA, listA, packpos)  
        implicit none
 
        include 'mpif.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nAmax
        integer, intent (in) :: nprowsmax
        integer, intent (in) :: packsize

! input byte array; should have size at least
! kind(integer)*(3+isendrows*(1+imax)) + kind(real)*2*isendrows*imax
! where 'imax' = the maximum nonzero element count in a column (this value can
! be obtained by calling sparse_getpacksize)
        integer*1, intent (in), dimension (packsize) :: packarray

! Output
! packsize = on input, the position within the array at which to begin packing;
! on output, number of bytes in the packed array plus initial packsize value.
        integer, intent (inout) :: packpos

        integer, intent (out), dimension (nprowsmax) :: numA
        integer, intent (out), dimension (nAmax, nprowsmax) :: listA
!FIXME: It is necessary to make these parameters volatile when using
!OpenMP.  I hope when module-private data is added to the standard
!this will not be necessary.
!$ volatile numA,listA
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

        integer, dimension (:), allocatable :: hold_indices

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Initialize the temporary array to zero:
        numA = 0
        listA = 0
       
        call MPI_UNPACK (packarray, packsize, packpos, maxcol, 1,            &
     &                   MPI_INTEGER, MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (packarray, packsize, packpos, isendrows, 1,         &
     &                   MPI_INTEGER, MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (packarray, packsize, packpos, itotal, 1,            &
     &                   MPI_INTEGER, MPI_BTN_WORLD, ierror)

        allocate (hold_indices (itotal))

        call MPI_UNPACK (packarray, packsize, packpos, numA, isendrows,      &
     &                   MPI_INTEGER, MPI_BTN_WORLD, ierror)
        call MPI_UNPACK (packarray, packsize, packpos, hold_indices, itotal, &
     &                   MPI_INTEGER, MPI_BTN_WORLD, ierror)

        itotal = 0
        do imu = 1, isendrows
         do index = 1, numA(imu)
          itotal = itotal + 1
          listA(index,imu) = hold_indices(itotal)
         end do
        end do

! Deallocate Arrays
! ===========================================================================
        deallocate (hold_indices)
 
! Format Statements
! ===========================================================================
 
        return
        end

