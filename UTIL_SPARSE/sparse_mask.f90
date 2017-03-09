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


! sparse_mask.f90
! Program Description
! ===========================================================================
!       This routine masks matrix A by matrix B; that is, the resulting
! matrix C is zero wherever B is zero, and equal to A wherever B is nonzero.
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
! Office telephone 801-257-9796
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine sparse_mask (nprows, nAmax, nBmax, nCmax, ncolsAmax,       &
     &                         ncolsBmax, ncolsCmax, numA, listA, elementsA, &
     &                         istartrowA, numB, listB,     &
     &                         istartrowB, numC, listC, elementsC)
        use interactions
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: istartrowA
        integer, intent (in) :: istartrowB
        integer, intent (in) :: nprows
        integer, intent (in) :: nAmax
        integer, intent (in) :: nBmax
        integer, intent (in) :: nCmax
        integer, intent (in) :: ncolsAmax
        integer, intent (in) :: ncolsBmax
        integer, intent (in) :: ncolsCmax

!$ volatile numA,listA,elementsA,numB,listB,numC,listC,elementsC

        integer, intent (in), dimension (ncolsAmax) :: numA
        integer, intent (in), dimension (ncolsBmax) :: numB
        integer, intent (in), dimension (nAmax, ncolsAmax) :: listA    
        integer, intent (in), dimension (nBmax, ncolsBmax) :: listB

        real, intent (in), dimension (nAmax, ncolsAmax) :: elementsA

! Output
        integer, intent (out), dimension (ncolsCmax) :: numC
        integer, intent (out), dimension (nCmax, ncolsCmax) :: listC

        real, intent (out), dimension (nCmax, ncolsCmax) :: elementsC

! Argument Declaration and Description
! ===========================================================================
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iA, iB
        integer imu
        integer indexA, indexB, indexC
        integer inu
        integer, dimension (:), allocatable :: indexloc

! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Initialize the resultant matrix to zero.
        numC = 0
        listC = 0
        elementsC = 0.0d0

! Set matrix C to matrix A masked by the index set B
! indexloc contains a nonzero index everywhere that the column in B is nonzero.
!$omp parallel do private(indexA,indexB,indexC,indexloc)
        do iA = 1, nprows
         allocate(indexloc(norbitals))
         indexloc = 0
         do inu = 1, numB(iA)
          indexB = listB(inu,iA)
          indexloc(indexB) = inu
         end do
         do imu = 1, numA(iA)
          indexA = listA(imu,iA) 
          if (indexloc(indexA).gt.0) then
           numC(iA) = numC(iA) + 1
           indexC = numC(iA)
           listC(indexC,iA) = indexA 
           elementsC(indexC,iA) = elementsA(imu,iA)
           indexloc(indexA) = 0
          end if
         end do
         deallocate(indexloc)
        end do

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
