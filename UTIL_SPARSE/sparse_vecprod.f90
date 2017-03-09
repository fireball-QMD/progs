! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University of Utah - James P. Lewis, Chair
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

 
! sparse_vecprod.f90
! Program Description
! ===========================================================================
!       This routine calculates the vector product of parts of matrix A and
! matrix B, treating them as 1-dimensional vectors.  The segments multiplied 
! have dimensions rows x 'nprows'.  The segment in A starts at row index 
! 'startrowA'.  The segment in A starts at row index 'startrowB'.  If 'reduce' 
! is true then the vector products are summed over all processors.
!
! ===========================================================================
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
        subroutine sparse_vecprod (nprows, nAmax, nBmax, ncolsAmax,          &
     &                             ncolsBmax, numA, listA, elementsA,        &
     &                             istartrowA, numB, listB, elementsB,       &
     &                             istartrowB, reduce, prod)
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
        integer, intent (in) :: istartrowB
        integer, intent (in) :: nBmax
        integer, intent (in) :: ncolsBmax
 
        integer, intent (in), dimension (ncolsAmax) :: numA
        integer, intent (in), dimension (ncolsBmax) :: numB
        integer, intent (in), dimension (nAmax, ncolsAmax) :: listA
        integer, intent (in), dimension (nBmax, ncolsBmax) :: listB

        real, intent (in), dimension (nAmax, ncolsAmax) :: elementsA
        real, intent (in), dimension (nBmax, ncolsBmax) :: elementsB

        logical, intent (in) :: reduce

!$ volatile numA,listA,elementsA,numB,listB,elementsB

! Output
        double precision, intent (out) :: prod

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer irow, iA, iB, pA, pB
        integer ierror
        double precision prodt1, prodt2, prodt3

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        prodt1 = 0.0d0
        prodt2 = 0.0d0
!$omp parallel do private(iA,iB,pA,pB) reduction(+:prodt1,prodt2)
        do irow = 1, nprows
         iA = 1
         iB = 1
         do while (iA .le. numA(irow) .and. listA(iA,irow) .lt. istartrowA)
          iA = iA + 1
         end do
         do while (iB .le. numB(irow) .and. listB(iB,irow) .lt. istartrowB)
          iB = iB + 1
         end do
         do while (iA .le. numA(irow) .and. iB .le. numB(irow))
          pA = listA(iA,irow) - istartrowA + 1
          pB = listB(iB,irow) - istartrowB + 1
          if (pA .lt. pB) then
           iA = iA + 1
          else if (pA .gt. pB) then
           iB = iB + 1
          else
           if (elementsA(iA,irow) * elementsB(iB,irow) .ge. 0) then
            prodt1 = prodt1 + elementsA(iA,irow)*elementsB(iB,irow)
           else
            prodt2 = prodt2 + elementsA(iA,irow)*elementsB(iB,irow)
           end if
           iA = iA + 1
           iB = iB + 1
          end if
         end do
        end do
        prodt3 = prodt1 + prodt2
        if (reduce) then
         call MPI_ALLREDUCE (prodt3, prod, 1, mpi_whatever_double, MPI_SUM,&
     &                       MPI_BTN_WORLD, ierror)
        else
         prod = prodt3
        end if

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
