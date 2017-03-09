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


! sparse_pack.f90
! Program Description
! ===========================================================================
!       This routine packs a sparse matrix into a byte array so it can be
! conveniently communicated.
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
        subroutine sparse_pack_elements (packsize, isendrows, nAmax,         &
     &                                   nprowsmax2, numA, elementsA,        &
     &                                   packarray, packpos)
        use mpi_declarations
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: isendrows
        integer, intent (in) :: nAmax
        integer, intent (in) :: nprowsmax2
        integer, intent (in) :: packsize

        integer, intent (in), dimension (nprowsmax2) :: numA

        real, intent (in), dimension (nAmax, nprowsmax2) :: elementsA

! Output
! output byte array; should have size at least 
! kind(integer)*(3+isendrows*(1+imax)) + kind(real)*2*isendrows*imax
! where 'imax' = the maximum nonzero element count in a column (this value can 
! be obtained by calling sparse_getpacksize)
        integer*1, intent (out), dimension (packsize) :: packarray

! packsize = on input, the position within the array at which to begin packing;
! on output, number of bytes in the packed array plus initial packsize value.
        integer, intent (inout) :: packpos
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ierror
        integer imu
        integer index
        integer itotal
        integer maxcol

        real, dimension (:), allocatable :: hold_elements 
 
! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Allocate Arrays
! ===========================================================================
        allocate (hold_elements (nAmax*isendrows))
 
! Procedure
! ===========================================================================
! Initialize the hold_indices and hold_elements to zero:
        hold_elements = 0

        itotal = 0
        maxcol = 0
        do imu = 1, isendrows
         if (numA(imu) .gt. maxcol) maxcol = numA(imu)
         do index = 1, numA(imu)
          itotal = itotal + 1
          hold_elements (itotal) = elementsA(index,imu)
         end do
        end do

! FIXME Would it be better to call MPI_PACK once per column and avoid the need
! for the hold arrays?
        call MPI_PACK (maxcol, 1, MPI_INTEGER, packarray, packsize, packpos, &
     &                 MPI_BTN_WORLD, ierror)
        call MPI_PACK (isendrows, 1, MPI_INTEGER, packarray, packsize,       &
     &                 packpos, MPI_BTN_WORLD, ierror)
        call MPI_PACK (itotal, 1, MPI_INTEGER, packarray, packsize, packpos, &
     &                 MPI_BTN_WORLD, ierror)
        call MPI_PACK (hold_elements, itotal, mpi_whatever_real, packarray,  &
     &                 packsize, packpos, MPI_BTN_WORLD, ierror)

! Deallocate Arrays
! ===========================================================================
        deallocate (hold_elements)
 
! Format Statements
! ===========================================================================
 
        return
        end

