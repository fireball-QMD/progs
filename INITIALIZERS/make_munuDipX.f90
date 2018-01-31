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
! along with this program.  if not, see <http://www.gnu.org/licenses/>.


! make_munuDipX.f90
! Program Description
! ===========================================================================
!       This subroutine calculates the following information (for all pairs
! of atoms (in1,in2)) :  (for the case of X dipole)
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
! (i.e. in mk_indexDipX.f)
! ===========================================================================
! Original code written by Jose Ortega.
 
! Code rewritten by:
! Jesús I. Mendieta-Moreno (JIMM)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine make_munuDipX (nspecies)
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
        allocate (index_max2cDipX (nspecies, nspecies))
 
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
 

! MEDipX_max is the maximum number of two-center matrix elements for X dipole: 

        ME2cDipX_max = 0

        do in1 = 1, nspecies
         do in2 = 1, nspecies
          index = 0
          do issh1 = 1, nssh(in1)
           l1 = lssh(issh1,in1)
           do issh2 = 1, nssh(in2)
            l2 = lssh(issh2,in2)
!         
           if ( l1 .eq. 0 .and. l2 .ne. 0 ) then
             index = index + 1
           end if
!         
           if ( l2 .eq. 0 .and. l1 .ne. 0 ) then
             index = index + 1
           end if
!         
           if ( l1 .eq. 1 .and. l2 .ne. 0 ) then
             index = index + 1
           end if
!         
           if ( l2 .eq. 1 .and. l1 .ne. 0 ) then
             index = index + 1
           end if
          
          
           if ( l1 .eq. 2 .and. l2 .ne. 0 ) then
             index = index + 1
             index = index + 1
           end if
!         
           if ( l2 .eq. 2 .and. l1 .ne. 0 ) then
             index = index + 1
             index = index + 1
           end if
          
           if ( l1 .eq. 2 .and. l2 .ne. 0 ) then
             index = index + 1
           end if
          
           if ( l2 .eq. 2 .and. l1 .ne. 0 ) then
             index = index + 1
           end if

           end do
          end do
          if (index .gt. ME2cDipX_max) ME2cDipX_max = index
          
         end do
       end do

 
! Allocate arrays
        allocate (muDipX (ME2cDipX_max, nspecies, nspecies))
        allocate (nuDipX (ME2cDipX_max, nspecies, nspecies))

        if (ME2cDipX_max .gt. ME2c_max) ME2c_max = ME2cDipX_max

! Now, calculate the mu-nu-map of the matrix elements on the box
! (num_orb(in1) x num_orb(in2)).  For each matrix-element (index) we calculate
! muDipX(index) and nuDipX(index).

        do in1 = 1, nspecies
         do in2 = 1, nspecies
          index = 0 
          n1 = 0
          do issh1 = 1, nssh(in1)
           l1 = lssh(issh1,in1)
           n1 = n1 + l1 + 1
           n2 = 0
           do issh2 = 1, nssh(in2)
            l2 = lssh(issh2,in2)
            n2 = n2 + l2 + 1

            if ( l1 .eq. 0 .and. l2 .ne. 0 ) then
              index = index + 1
              muDipX(index,in1,in2) = n1
              nuDipX(index,in1,in2) = n2 + 1
            end if
!
            if ( l2 .eq. 0 .and. l1 .ne. 0 ) then
              index = index + 1
              muDipX(index,in1,in2) = n1 + 1
              nuDipX(index,in1,in2) = n2 
            end if
!
            if ( l1 .eq. 1 .and. l2 .ne. 0 ) then
              index = index + 1
              muDipX(index,in1,in2) = n1
              nuDipX(index,in1,in2) = n2 + 1
            end if
!
            if ( l2 .eq. 1 .and. l1 .ne. 0 ) then
              index = index + 1
              muDipX(index,in1,in2) = n1 + 1
              nuDipX(index,in1,in2) = n2 
            end if


            if ( l1 .eq. 2 .and. l2 .ne. 0 ) then
              index = index + 1
              muDipX(index,in1,in2) = n1
              nuDipX(index,in1,in2) = n2 + 1

              index = index + 1
              muDipX(index,in1,in2) = n1 + 2
              nuDipX(index,in1,in2) = n2 + 1
            end if
!
            if ( l2 .eq. 2 .and. l1 .ne. 0 ) then
              index = index + 1
              muDipX(index,in1,in2) = n1 + 1
              nuDipX(index,in1,in2) = n2 

              index = index + 1
              muDipX(index,in1,in2) = n1 + 1
              nuDipX(index,in1,in2) = n2 + 2
            end if

            if ( l1 .eq. 2 .and. l2 .ne. 0 ) then
              index = index + 1
              muDipX(index,in1,in2) = n1 - 2
              nuDipX(index,in1,in2) = n2 - 1
            end if

            if ( l2 .eq. 2 .and. l1 .ne. 0 ) then
              index = index + 1
              muDipX(index,in1,in2) = n1 - 1
              nuDipX(index,in1,in2) = n2 - 2
            end if

            n2 = n2 + l2
           end do
           n1 = n1 + l1
          end do
          index_max2cDipX(in1,in2) = index

! end loops over the species
         end do
        end do
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
