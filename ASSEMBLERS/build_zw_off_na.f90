! copyright info:
!
!
!                             @Copyright 2004
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Universidad de Madrid - Pavel Jelinek
! Brigham Young University - Hao Wang

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio State University - Dave Drabold
! University of Regensburg - Juergen Fritsch

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


! build_olsxc_off.f90
! Program Description
! ===========================================================================
!       This routine takes the input density and determines the
! exchange-correlation potential and energy atom-like matrix elements 
! in the average-density approx. (see PRB 40, 3979 (1989) ).

! Formula:
! <i,mu|V_xc(n)|j,nu> =
! = V_xc(n^a)<i,mu|j,nu> + V'_xc(n^a)[<i,mu|n|j,nu> - n^a<i,mu|j,nu>] -
!   - <i,mu|V_xc(n_i+n_j)|j,nu> + V_xc({n_i_n_j}^a)<i,mu|j,nu> +
!   + V'_xc({n_i+n_j}^a)[<i,mu|n_i+n_j|j,nu> - ({n_i+n_j}^a)<i,mu|j,nu>] }
! where
!  V'_xc = d(V_xc)/dn
!  n^a ... average density
!  n_i ... density on i-site
!  n_j ... density on j-site
! 
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162  
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
!
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine build_zw_off_na (in1, in2, den1x, denx, sx, ineigh,       &
     &                              iatom, bcxcx)
        use charges
        use density
        use interactions

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: in1
        integer, intent (in) :: in2
        integer, intent (in) :: ineigh
        integer, intent (in) :: iatom

        real, intent (in), dimension (numorb_max, numorb_max) :: denx
        real, intent (in), dimension (numorb_max, numorb_max) :: den1x
        real, intent (in), dimension (numorb_max, numorb_max) :: sx

! Output
        real, intent (out), dimension (numorb_max, numorb_max) :: bcxcx

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer ind1
        integer ind2
        integer inu
        integer issh
        integer jssh
        integer l1
        integer l2
        integer n1
        integer n2

        real dexc
        real d2exc
        real dmuxc
        real d2muxc
        real exc
        real muxc
        real dmuxcij
        real muxcij
        
        real, dimension (nsh_max, nsh_max) :: dens
        real, dimension (nsh_max, nsh_max) :: densij

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Initialize
        bcxcx = 0.0d0

! set average density
        do issh = 1, nssh(in1)
           do jssh = 1, nssh(in2)
              dens(issh,jssh) = arho_off(issh,jssh,ineigh,iatom)
              densij(issh,jssh) = arhoij_off(issh,jssh,ineigh,iatom)
           enddo
        enddo

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CALCULATE AVERAGE DENSITY
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        n1 = 0

! Loop over shells of i-atom
        do issh = 1, nssh(in1)

! Number of orbitals per the shell
         l1 = lssh(issh,in1)
         n1 = n1 + l1 + 1
         n2 = 0

! Loop over shells of j-atom 
         do jssh = 1, nssh(in2)

! Number of orbitals per the shell 
          l2 = lssh(jssh,in2)
          n2 = n2 + l2 + 1
 
! Calculate XC potentials with average density 
          call cepal (dens(issh,jssh), exc, muxc, dexc, d2exc, dmuxc, d2muxc)
          call cepal (densij(issh,jssh), exc, muxcij, dexc, d2exc, dmuxcij,     &
     &                d2muxc)


! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! BUILD XC-MATRIX
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Calculate the 1st order approximation for the vxc-matrix elements
! (and 2nd order corrections for the "diagonal" terms)
! Comment: it can be moved in previous cycle, but for clarity ...
!
!  Formula:
!  <i,mu|V_xc(n)|j,nu>  
!        V_xc(n^a)<i,mu|j,nu> + V'_xc(n^a)[<i,mu|n|j,nu> - n^a<i,mu|j,nu>] -
!         - <i,mu|V_xc(n_i+n_j)|j,nu> + V_xc([n_i+n_j]^a)<i,mu|j,nu> +
!         + V'_xc([n_i+n_j]^a)[<i,mu|n_i+n_j|j,nu> - [n_i+n_j]^a<i,mu|j,nu>] }
!  where
!        V_xc(n^a)    ... muxc
!        V'_xc(n^a)   ... dmuxc
!        V_xc([n_i+n_j]^a)  ... muxcij
!        V'_xc([n_i+n_j]^a) ... dmuxcij
!        <i,mu|j,nu>  ... sx 
!        n^a          ... dens
!        [n_i+n_j]^a  ... densij
!
! keeping in mind that <i,mu|V_xc(n_i+n_j)|j,nu> is done in assemble_2c() !!!  

! Set the XC-submatrices (in molecular coord)
          ! loop over orbitals in the x-shell
          do ind1 = -l1, l1
           imu = n1 + ind1
           ! loop over orbitals in the y-shell
           do ind2 = -l2, l2
            inu = n2 + ind2

            ! V_xc(n^a)<i,mu|j,nu>
            bcxcx(imu,inu) = muxc*sx(imu,inu)
            ! V'_xc(n^a)[<i,mu|n|j,nu> - n^a<i,mu|j,nu>] 
            bcxcx(imu,inu) = bcxcx(imu,inu)                                  &
     &       + dmuxc*(denx(imu,inu) - dens(issh,jssh)*sx(imu,inu))
            ! - V_xc([n_i+n_j]^a)<i,mu|j,nu>  
            ! - V'_xc([n_i+n_j]^a){<i,mu|n_i+n_j|i,nu> 
            !     - [n_i+n_j]^a<i,mu|i,nu>]}
            bcxcx(imu,inu) = bcxcx(imu,inu) - muxcij*sx(imu,inu)              &
     &       - dmuxcij*(den1x(imu,inu) - densij(issh,jssh)*sx(imu,inu)) 
           end do ! do ind2 = -l2, l2
          end do ! do ind1 = -l1, l1

          ! increment indicator 
          n2 = n2 + l2
         end do ! do jssh = 1, nssh(in2)

         ! increment indicator 
         n1 = n1 + l1
        end do !do issh = 1, nssh(in1)

! Deallocate Arrays
! ===========================================================================
! Format Statements
! ===========================================================================

        return
      end subroutine build_zw_off_na


