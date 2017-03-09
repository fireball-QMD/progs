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


! writeout_comph.f90
! Program Description
! ===========================================================================
!       This routine is a utility to write out the components of the 
! Hamiltonian and compares the pieces with the total band-structure energy.
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
        subroutine writeout_comph (natoms, itheory, ebs)
        use density
        use dimensions
        use interactions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input 
        integer, intent (in) :: itheory
        integer, intent (in) :: natoms

        real, intent (in) :: ebs

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer imu
        integer in1, in2
        integer ineigh
        integer inu
        integer jatom

        real comp_EBS
        real comp_EWDLR
        real comp_EWDSR
        real comp_KE
        real comp_VCA
        real comp_VNA
        real comp_VNL
        real comp_VNAa
        real comp_VNAo
        real comp_VXC
        real comp_VXC_CA
        real comp_VXC_1C

        real tmp

! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! We will compute the individual components of the total energy. Compute
! rho(mu,nu,i,m)*matrix_element_of_h(imu,inu,iatom,ineigh) to obtain this.
! Changes of this can then be directly compared with various force components.
! Calculate and write out the components only if we are interested in pieces.
        write (*,*) '  '
        write (*,*) ' ----------------------------------------------------- '
        write (*,*) ' Compute the various components of the total '
        write (*,*) ' band-structure energy. '
        write (*,*) ' ----------------------------------------------------- '

! We use comp to mean component. So compKE means Kinetic Energy component
! to the total energy.  First initialize everything to zero.
        comp_KE = 0.0d0
        comp_VNA = 0.0d0
        comp_VXC = 0.0d0
        comp_VNAa = 0.0d0
        comp_VNAo = 0.0d0
        comp_VNL = 0.0d0
        comp_VCA = 0.0d0
        comp_VXC_CA = 0.0d0
        comp_EWDLR = 0.0d0
        comp_EWDSR = 0.0d0
        comp_EBS = 0.0d0
        comp_VXC_1C = 0.0d0
!oslxc-forces-test

        do iatom = 1, natoms
          in1 = imass(iatom)
          do ineigh = 1, neighn(iatom)
           jatom = neigh_j(ineigh,iatom)
           in2 = imass(jatom)
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             comp_KE =                                                       &
     &        comp_KE + rho(imu,inu,ineigh,iatom)*t_mat(imu,inu,ineigh,iatom)
             if (ineigh .eq. neigh_self(iatom)) then
              comp_VNAa =                                                    &
     &         comp_VNAa + rho(imu,inu,ineigh,iatom)*vna(imu,inu,ineigh,iatom)
             else
              comp_VNAo =                                                    &
     &         comp_VNAo + rho(imu,inu,ineigh,iatom)*vna(imu,inu,ineigh,iatom)
             end if
             comp_VNA =                                                      &
     &        comp_VNA + rho(imu,inu,ineigh,iatom)*vna(imu,inu,ineigh,iatom)
             comp_VXC =                                                      &
     &        comp_VXC + rho(imu,inu,ineigh,iatom)*vxc(imu,inu,ineigh,iatom)
             comp_EBS =                                                      &
     &        comp_EBS + rho(imu,inu,ineigh,iatom)*h_mat(imu,inu,ineigh,iatom)
             comp_VXC_1C =                                                   &
     &        comp_VXC_1C + rho(imu,inu,ineigh,iatom)*vxc_1c(imu,inu,ineigh,iatom)
             if (itheory .eq. 1 .or. itheory .eq. 2) then
              comp_VXC_CA = comp_VXC_CA                                      &
     &          + rho(imu,inu,ineigh,iatom)*vxc_ca(imu,inu,ineigh,iatom)
              comp_VCA =                                                     &
     &         comp_VCA + rho(imu,inu,ineigh,iatom)*vca(imu,inu,ineigh,iatom)
              comp_EWDLR = comp_EWDLR                                        &
     &         + rho(imu,inu,ineigh,iatom)*ewaldlr(imu,inu,ineigh,iatom)
              comp_EWDSR = comp_EWDSR                                        &
     &        - rho(imu,inu,ineigh,iatom)*ewaldsr(imu,inu,ineigh,iatom)
             end if
            end do
           end do
          end do
        end do

 
! HAO Sep. 10, 2003. ADD FOR VNL
        do iatom = 1, natoms
         in1 = imass(iatom)
! Now loop over all neighbors jatom of iatom
         do ineigh = 1, neighPPn(iatom)
          jatom = neighPP_j(ineigh,iatom)
          in2 = imass(jatom)
          do inu = 1, num_orb(in2)
           do imu = 1, num_orb(in1)
           comp_EBS = comp_EBS                                              &
      &      + rhoPP(imu,inu,ineigh,iatom)*vnl(imu,inu,ineigh,iatom)

           comp_VNL = comp_VNL                                              &
      &      + rhoPP(imu,inu,ineigh,iatom)*vnl(imu,inu,ineigh,iatom)

           end do ! do imu
          end do ! do inu
         end do ! do inegh
        enddo ! do iatom

        write (*,*) ' comp_KE     = ', comp_KE
        write (*,*) ' comp_VNA    = ', comp_VNA
        write (*,*) ' comp_VXC    = ', comp_VXC
        write (*,*) ' comp_VNL    = ', comp_VNL
        write (*,*) ' comp_VXC_1C = ', comp_VXC_1C
        write (*,*) ' comp_VCA    = ', comp_VCA
        write (*,*) ' comp_VXC_CA = ', comp_VXC_CA
        write (*,*) ' comp_EWDLR  = ', comp_EWDLR
        write (*,*) ' comp_EWDSR  = ', comp_EWDSR
        write (*,*) '   '
        write (*,*) ' comp_EBS    = ', comp_EBS, ' EBS = ', ebs
        write (*,*) '               error = ', abs(ebs - comp_EBS)


! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
      end subroutine writeout_comph
