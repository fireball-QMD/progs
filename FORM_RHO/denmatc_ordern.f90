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


! denmatb_ordern.f90
! Program Description
! ===========================================================================
!       This routine calculates the charges from the given density matrices.
!
! ===========================================================================
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
        subroutine denmatc_ordern (natoms, ifixcharge, iqout)
        use charges
        use density
        use interactions
        use neighbor_map
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ifixcharge
        integer, intent (in) :: iqout
        integer, intent (in) :: natoms

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer imu
        integer in1, in2
        integer ineigh
        integer inu
        integer issh
        integer jatom
        integer jneigh
        integer mqn

        real, dimension (numorb_max, natoms) :: QMulliken

! Procedure
! ===========================================================================
! ****************************************************************************
!
!  C O M P U T E    M U L L I K E N    C H A R G E S
! ****************************************************************************
        if (iqout .eq. 2) then
 
! Compute Mulliken charges.
         if (ifixcharge .eq. 1) then
          do iatom = 1, natoms
           in1 = imass(iatom)
           QMulliken_TOT(iatom) = 0.0d0
           do issh = 1, nssh(in1)
            Qout(issh,iatom) = Qin(issh,iatom) 
            QMulliken_TOT(iatom) = QMulliken_TOT(iatom) + Qin(issh,iatom) 
           end do
          end do
         else
          do iatom = 1, natoms
           in1 = imass(iatom)
           do issh = 1, nssh(in1)
            Qout(issh,iatom) = 0.0d0
           end do
 
! Initialize
           QMulliken_TOT(iatom) = 0.0d0
           QMulliken(:,iatom) = 0.0d0
 
! Loop over neighbors
           do ineigh = 1, neighn(iatom)
            jatom = neigh_j(ineigh,iatom)
            in2 = imass(jatom)
            jneigh = neigh_back(iatom,ineigh)
 
            do imu = 1, num_orb(in1)
             do inu = 1, num_orb(in2)
              QMulliken(imu,iatom) = QMulliken(imu,iatom)                    &
     &         + 0.5d0*(rho(imu,inu,ineigh,iatom)*s_mat(imu,inu,ineigh,iatom)&
     &                + rho(inu,imu,jneigh,jatom)*s_mat(inu,imu,jneigh,jatom))
             end do
            end do
 
! End loop over neighbors
           end do
 
! Finally the imu loop.
           imu = 0
           do issh = 1, nssh(in1)
            do mqn = 1, 2*lssh(issh,in1) + 1
             imu = imu + 1
             Qout(issh,iatom) = Qout(issh,iatom) + QMulliken(imu,iatom)
             QMulliken_TOT(iatom) = QMulliken_TOT(iatom) + QMulliken(imu,iatom)
            end do
           end do
 
! End loop over atoms
          end do
         end if
        end if
 
 
! Deallocate Arrays
! ===========================================================================

 
! Format Statements
! ===========================================================================
 
        return
        end
