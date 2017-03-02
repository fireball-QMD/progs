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

! sparse_add.f90
! Program Description
! ===========================================================================
!       This routines adds segments of two matrices, represented in a compact 
! form, together.  A multiplication factor for the second matrix may be 
! included, such that C = multA*A + multB*B.
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
        subroutine sparse_add (nprows, nAmax, nBmax, nCmax, ncolsAmax,       &
     &                         ncolsBmax, ncolsCmax, numA, listA, elementsA, &
     &                         multA, istartrowA, numB, listB, elementsB,    &
     &                         multB, istartrowB, numC, listC, elementsC)
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

!$ volatile numA,listA,elementsA,numB,listB,elementsB,numC,listC,elementsC

        integer, intent (in), dimension (ncolsAmax) :: numA
        integer, intent (in), dimension (ncolsBmax) :: numB
        integer, intent (in), dimension (nAmax, ncolsAmax) :: listA    
        integer, intent (in), dimension (nBmax, ncolsBmax) :: listB

        real, intent (in) :: multA
        real, intent (in) :: multB
        real, intent (in), dimension (nAmax, ncolsAmax) :: elementsA
        real, intent (in), dimension (nBmax, ncolsBmax) :: elementsB

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

        logical match

! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Initialize the resultant matrix to zero.
        numC = 0
        listC = 0
        elementsC = 0.0d0

! Add matrices C = multA*A + multB*B
! indexloc contains a nonzero index everywhere that the column in B is nonzero.
!$omp parallel do private(indexA,indexB,indexC,match,indexloc)
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
           elementsC(indexC,iA) =                                           &
     &      multA*elementsA(imu,iA) + multB*elementsB(indexloc(indexA),iA)
           indexloc(indexA) = 0
          else

! Get all elements of matrix A where matrixB is zero.
           numC(iA) = numC(iA) + 1
           indexC = numC(iA)
           listC(indexC,iA) = indexA 
           elementsC(indexC,iA) = multA*elementsA(imu,iA)
          end if
         end do
! Get all elements of matrixB where matrixA is zero. 
         if (nCmax .gt. nBmax) then
          do imu = 1, norbitals
           if (indexloc(imu).gt.0) then
            numC(iA) = numC(iA) + 1
            indexC = numC(iA)
            listC(indexC,iA) = imu
            elementsC(indexC,iA) = multB*elementsB(indexloc(imu),iA)
           end if
          end do
         end if
         deallocate(indexloc)
        end do

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
