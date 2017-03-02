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

! sparse_mult.f90
! Program Description
! ===========================================================================
! Multiplies matrix A^T by matrix B^T to get matrix C^T.
! Matrix A has the following format:
!  numA contains npcols integers indicating the # of nonzero elements in 
! each column of A.
!  listA contains for each column of A a list of the positions of the nonzero 
! elements.
!  elementsA contains for each column of A a list of the elements in those 
! positions.

! Matrices B and C have similar structure.
! A and B are transposed before they are multiplied, and the output matrix C 
! is the transpose of A^T*B^T.

!  npcols = number of columns to operate on in A
!  nAmax, nBmax, nCmax = row dimensions of A, B, C
!  ncolsAmax, ncolsBmax, ncolsCmax = column dimensions of A, B, C; ncolsCmax 
! must be >= npcols 
!  ncolsB = number of columns to operate on in B; should be comparable to 
! nAmax, since they define the common dimension of A and B.
!  nCmax must be >= nBmax
!  append = if true, the product is added to the matrix already contained 
! in C; if false, C is set to zero first
!  need_list = if true, numC and listC are set up and returned; if false, the 
! values passed in numC and listC are used and they are not changed.
!  istartrowA = currently ignored, should be removed
!  istartrowB = the row within A at which the multiplication should begin; 
! allows extracting and muliplying part of a matrix (the name should be changed)
 
! ===========================================================================
! Original code written by Spencer Shellman.
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
        subroutine sparse_mult (npcols, nAmax, nBmax, nCmax, ncolsAmax,      &
     &                          ncolsB, ncolsBmax, ncolsCmax, numA, listA,   &
     &                          elementsA, istartrowA, numB, listB,          &
     &                          elementsB, istartrowB, append, need_list,    &
     &                          numC, listC, elementsC)         
        use interactions
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: istartrowA
        integer, intent (in) :: istartrowB
        integer, intent (in) :: npcols
        integer, intent (in) :: nAmax
        integer, intent (in) :: nBmax
        integer, intent (in) :: nCmax
        integer, intent (in) :: ncolsAmax
        integer, intent (in) :: ncolsB
        integer, intent (in) :: ncolsBmax
        integer, intent (in) :: ncolsCmax

!$ volatile numA, listA, elementsA, numB, listB, elementsB, numC, listC
!$ volatile elementsC, istartrowA, istartrowB, npcols, nAmax, nBmax, nCmax
!$ volatile ncolsAmax, ncolsB, ncolsBmax, ncolsCmax, append, need_list
 
        integer, intent (in), dimension (ncolsAmax) :: numA
        integer, intent (in), dimension (ncolsBmax) :: numB
        integer, intent (in), dimension (nAmax, ncolsAmax) :: listA
        integer, intent (in), dimension (nBmax, ncolsBmax) :: listB

        real, intent (in), dimension (nAmax, ncolsAmax) :: elementsA
        real, intent (in), dimension (nBmax, ncolsBmax) :: elementsB

        logical, intent (in) :: append
        logical, intent (in) :: need_list

! Output
        integer, intent (inout), dimension (ncolsCmax) :: numC
        integer, intent (inout), dimension (nCmax, ncolsCmax) :: listC

        real, intent (inout), dimension (nCmax, ncolsCmax) :: elementsC

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iA, iB, iC
        integer imu
        integer index1, index2
        integer inu
        integer mmu

        integer, dimension (norbitals) :: indexloc
        integer, dimension (norbitals) :: indexv

! separate buffers for adding positive and negative numbers, for stability
        double precision, dimension (norbitals) :: innerprod1, innerprod2
        double precision r8a, r8b

! Allocate Arrays
! ===========================================================================
! Procedure
! ===========================================================================
! Initialize the resultant matrix to zero.
        if (.not. append) then 
         if (need_list) then 
          numC = 0
          listC = 0
         end if
         elementsC = 0.0d0
        end if
! First, set up indices for C sparse matrix.
        if (need_list) then
! The parameter indexloc contains a 1 for each index at which a nonzero 
! element is contained in the column indexv holds the locations of nonzero 
! elements in the column
         indexloc = 0
         indexv = 0
!$omp parallel do private(index1,index2,iB) firstprivate(indexloc,indexv)
         do iA = 1, npcols
          index2 = numC(iA)
          if (append) then               ! If we are appending, then do not 
           do imu = 1, index2            ! count this index in listC again.
            index1 = listC(imu,iA)  
            indexloc(index1) = 1
            indexv(imu) = index1
           end do
          end if
          do imu = 1, numA(iA)
           iB = listA(imu,iA)
           if (iB .ge. istartrowB .and. iB .le. ncolsB + istartrowB - 1) then
            iB = iB - istartrowB + 1
            do inu = 1, numB(iB)
             index1 = listB(inu,iB)
             if (indexloc(index1) .eq. 0) then
              indexloc(index1) = 1
              index2 = index2 + 1
              indexv(index2) = index1
             end if
            end do
           end if
          end do
          numC(iA) = index2 
          do imu = 1, index2
           index1 = indexv(imu)
           indexv(imu) = 0
           indexloc(index1) = 0
           listC(imu,iA) = index1
          end do
         end do
        end if

! Second, multiply matrices C = A*B
! innerprod accumulates each column of C
        innerprod1 = 0.0d0
        innerprod2 = 0.0d0
!$omp parallel do private(iB,iC,index1,r8a,r8b) firstprivate(innerprod1,innerprod2)
        do iA = 1, npcols
         do imu = 1, numA(iA)
          iB = listA(imu,iA)
          if (iB .ge. istartrowB .and. iB .le. ncolsB + istartrowB - 1) then
           iB = iB - istartrowB + 1
           do inu = 1, numB(iB)
            index1 = listB(inu,iB)
            r8a = elementsA(imu,iA)
            r8b = elementsB(inu,iB)
            if (r8a*r8b.ge.0) then
               innerprod1(index1) = innerprod1(index1) + r8a*r8b
            else
               innerprod2(index1) = innerprod2(index1) + r8a*r8b
            end if
           end do
          end if
         end do
         do imu = 1, numC(iA)
          iC = listC(imu,iA)
          if (append) then
           elementsC(imu,iA) =                                               &
     &      elementsC(imu,iA) + innerprod1(iC) + innerprod2(iC)
          else 
           elementsC(imu,iA) = innerprod1(iC) + innerprod2(iC)
          end if
          innerprod1(iC) = 0.0d0
          innerprod2(iC) = 0.0d0
         end do
        end do

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end

