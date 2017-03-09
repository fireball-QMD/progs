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


! sparse_norm2.f90
! Program Description
! ===========================================================================
!       This routine calculates the norm product of a vector.
!
! ===========================================================================
! Code written by:
! Spencer Shellman
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine sparse_norm2 (nprows, nAmax, ncolsAmax, numA, listA,      &
     &                           elementsA, istartrowA, reduce, prod)
        use interactions
        use mpi_declarations
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
          integer, intent (in) :: nprows
          integer, intent (in) :: istartrowA
          integer, intent (in) :: nAmax
          integer, intent (in) :: ncolsAmax
 
          integer, intent (in), dimension (ncolsAmax) :: numA
          integer, intent (in), dimension (nAmax, ncolsAmax) :: listA

          real, intent (in), dimension (nAmax, ncolsAmax) :: elementsA

          logical, intent (in) :: reduce

!$ volatile numA,listA,elementsA

! Output
          double precision, intent (out) :: prod

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
          integer irow, iA, pA
          integer ierror
          double precision prodt

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
          prodt = 0.0d0
!$omp parallel do private(iA,pA) reduction(+:prodt)
          do irow = 1, nprows
           iA = 1
           do while (iA .le. numA(irow) .and. listA(iA,irow) .lt. istartrowA)
            iA = iA + 1
           end do
           do while (iA .le. numA(irow))
            pA = listA(iA,irow) - istartrowA + 1
            prodt = prodt + elementsA(iA,irow)*elementsA(iA,irow)
            iA = iA + 1
           end do
          end do
          if (reduce) then
           call MPI_ALLREDUCE (prodt, prod, 1, mpi_whatever_double, MPI_SUM,&
     &                         MPI_BTN_WORLD, ierror)
          else
           prod = prodt
          end if

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
          return
          end subroutine sparse_norm2
