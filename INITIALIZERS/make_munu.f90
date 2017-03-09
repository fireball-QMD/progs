! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young of Utah - James P. Lewis, Chair
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


! make_munu.f90
! Program Description
! ===========================================================================
!       This subroutine calculates the following information (for all pairs
! of atoms (in1,in2)) :
!
! num_orb (in1) : number of orbitals in atom-type in1
! mu (index,in1,in2) : the mu-position for each matrix-element (index) between
!                      atom-type in1 and atom-type in2
! nu (index,in1,in2) : the nu-position for each matrix-element (index) between
!                      atom-type in1 and atom-type in2
 
! (on the BOX ( num_orb(in1) x num_orb(in2)))
!
! Atoms 1 and 2 (bondcharge) are along the z-axis; the third atom is in the
! xz-plane. The labelling of the orbitals is as follows:
!
!   S-shell :                s
!                            0
!
!   P-shell :           py   pz   px
!                       -1   0    +1
!
!   D-shell :     xy    yz   z^2  xz   x^2-y^2
!                 -2    -1   0    +1     +2
!
! ***************************************************************************
! JOM: IMPORTANT: in this subroutine, the order of the different matrix
! elements (index) is the SAME as the one given in CREATOR
! (i.e. in MK_INDEX.f)
! ===========================================================================
! Original code written by Jose Ortega.
 
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
        subroutine make_munu (nspecies)
        use dimensions
        use interactions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nspecies
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer index
        integer in1, in2
        integer issh1, issh2
        integer l1, l2
        integer n1, n2
 
! Allocate Arrays
! ===========================================================================
        allocate (num_orb (nspecies))
        allocate (index_max2c (nspecies, nspecies))
        allocate (index_max3c (nspecies, nspecies))
 
! Procedure
! ===========================================================================
! First, calculate the number of orbitals in each atom-type
        do in1 = 1, nspecies
         num_orb(in1) = 0
         do issh1 = 1 , nssh(in1)
          l1 = lssh(issh1,in1)
          num_orb(in1) = num_orb(in1) + 2*l1 + 1
         end do
        end do
 
! Now we calculate ME3c_max and ME2c_max
! ME2c_max is the maximum number of two-center matrix elements: 
! Examples: s ==> 1, sp^3 ==>  6, ss*p^3p*^3 ==> 24, sp^3d^5 ==> 19
! ME3c_max is the maximum number of three-center matrix elements:
! Examples: s ==> 1, sp^3 ==> 10, ss*p^3p*^3 ==> 40, sp^3d^5 ==> 45
        ME2c_max = 0
        ME3c_max = 0

        do in1 = 1, nspecies
         do in2 = 1, nspecies
          index = 0
          do issh1 = 1 , nssh(in1)
           l1 = lssh(issh1,in1)
           do issh2 = 1, nssh(in2)
            l2 = lssh(issh2,in2)
            do imu = -min(l1,l2), min(l1,l2)
             index = index + 1
            end do
           end do
          end do
          if (index .gt. ME2c_max) ME2c_max = index
          do issh1 = 1, nssh(in1)
           l1 = lssh(issh1,in1)
           do issh2 = 1, nssh(in2)
            l2 = lssh(issh2,in2)
 
            if (l1 .eq. 0 .and. l2 .ne. 0) then
             index = index + 1
            end if
 
            if (l1 .eq. 1) then
             index = index + 1
             if (l2 .ne. 0) then
              index = index + 1
             end if
             if (l2 .eq. 2) then
              index = index + 2
             end if
            end if
 
            if (l1 .eq. 2) then
             index = index + 1
             if (l2 .ne. 0) then
              index = index + 3
             end if
             if (l2 .eq. 2) then
              index = index + 2
             end if
            end if

           end do
          end do
 
          do issh1 = 1, nssh(in1)
           l1 = lssh(issh1,in1)
           do issh2 = 1, nssh(in2)
            l2 = lssh(issh2,in2)
 
            if (l1 .eq. 2) then
             index = index + 1
            end if
 
            if (l2 .eq. 2) then
             index = index + 1
            end if

           end do
          end do
          if (index .gt. ME3c_max) ME3c_max = index
 
         end do
        end do
        
! Allocate arrays
        allocate (mu (ME3c_max, nspecies, nspecies))
        allocate (mvalue (ME3c_max, nspecies, nspecies))
        allocate (nu (ME3c_max, nspecies, nspecies))

! Now, calculate the mu-nu-map of the matrix elements on the box
! (num_orb(in1) x num_orb(in2)).  For each matrix-element (index) we calculate
! mu(index) and nu(index).
 
! First, the non-zero two-center (and also three-center) interactions
        do in1 = 1, nspecies
         do in2 = 1, nspecies
          index = 0
          n1 = 0
          do issh1 = 1 , nssh(in1)
           l1 = lssh(issh1,in1)
           n1 = n1 + l1 + 1
           n2 = 0
           do issh2 = 1, nssh(in2)
            l2 = lssh(issh2,in2)
            n2 = n2 + l2 + 1
            do imu = -min(l1,l2), min(l1,l2)
             index = index + 1
             mu(index,in1,in2) = n1 + imu
             nu(index,in1,in2) = n2 + imu
             mvalue(index,in1,in2) = 0
            end do
            n2 = n2 + l2
           end do
           n1 = n1 + l1
          end do
          index_max2c(in1,in2) = index
 
! For three-center interactions there are some extra non-zero interactions; in
! this case, we find out (from symmetry considerations) that non-negative
! values of m_i mix with non-negative values of m_j.  Also, negative values of
! m_i mix with the negative values of m_j.
!
! Case 1, the interactions with M1 = M2 +- 1      (M=1 case)
          n1 = 0
          do issh1 = 1, nssh(in1)
           l1 = lssh(issh1,in1)
           n1 = n1 + l1 + 1
           n2 = 0
           do issh2 = 1, nssh(in2)
            l2 = lssh(issh2,in2)
            n2 = n2 + l2 + 1
 
            if (l1 .eq. 0 .and. l2 .ne. 0) then
             index = index + 1
             mu(index,in1,in2) = n1
             nu(index,in1,in2) = n2 + 1
             mvalue(index,in1,in2) = 1
            end if
 
            if (l1 .eq. 1) then
             index = index + 1
             mu(index,in1,in2) = n1 + 1
             nu(index,in1,in2) = n2
             mvalue(index,in1,in2) = 1
 
             if (l2 .ne. 0) then
              index = index + 1
              mu(index,in1,in2) = n1
              nu(index,in1,in2) = n2 + 1
              mvalue(index,in1,in2) = 1
             end if
 
             if (l2 .eq. 2) then
              index = index + 1
              mu(index,in1,in2) = n1 + 1
              nu(index,in1,in2) = n2 + 2
              mvalue(index,in1,in2) = 1
 
              index = index + 1
              mu(index,in1,in2) = n1 - 1
              nu(index,in1,in2) = n2 - 2
              mvalue(index,in1,in2) = 1
             end if
            end if
 
            if (l1 .eq. 2) then
             index = index + 1
             mu(index,in1,in2) = n1 + 1
             nu(index,in1,in2) = n2
             mvalue(index,in1,in2) = 1
 
             if (l2 .ne. 0) then
              index = index + 1
              mu(index,in1,in2) = n1
              nu(index,in1,in2) = n2 + 1
              mvalue(index,in1,in2) = 1
 
              index = index + 1
              mu(index,in1,in2) = n1 - 2
              nu(index,in1,in2) = n2 - 1
              mvalue(index,in1,in2) = 1
 
              index = index + 1
              mu(index,in1,in2) = n1 + 2
              nu(index,in1,in2) = n2 + 1
              mvalue(index,in1,in2) = 1
             end if
 
             if (l2 .eq. 2) then
              index = index + 1
              mu(index,in1,in2) = n1 + 1
              nu(index,in1,in2) = n2 + 2
              mvalue(index,in1,in2) = 1
 
              index = index + 1
              mu(index,in1,in2) = n1 - 1
              nu(index,in1,in2) = n2 - 2
              mvalue(index,in1,in2) = 1
             end if
            end if
            n2 = n2 + l2
           end do
           n1 = n1 + l1
          end do
 
! Case 2, the interactions with M1 = M2 +- 2      (M=2 case)
          n1 = 0
          do issh1 = 1, nssh(in1)
           l1 = lssh(issh1,in1)
           n1 = n1 + l1 + 1
           n2 = 0
           do issh2 = 1, nssh(in2)
            l2 = lssh(issh2,in2)
            n2 = n2 + l2 + 1
 
            if (l1 .eq. 2) then
             index = index + 1
             mu(index,in1,in2) = n1 + 2
             nu(index,in1,in2) = n2
             mvalue(index,in1,in2) = 2
            end if
 
            if (l2 .eq. 2) then
             index = index + 1
             mu(index,in1,in2) = n1
             nu(index,in1,in2) = n2 + 2
             mvalue(index,in1,in2) = 2
            end if
            n2 = n2 + l2
           end do
           n1 = n1 + l1
          end do
          index_max3c(in1,in2) = index
 
! End loops over the species
         end do
        end do
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
