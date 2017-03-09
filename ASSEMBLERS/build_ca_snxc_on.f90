! copyright info:
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

 
! build_ca_snxc_on.f90
! Program Description
! ===========================================================================
!       This routine takes the input density and determines the
! exchange-correlation potential and energy atom-like matrix elements 
! in the average-density approx. (see PRB 40, 3979 (1989) ).
!
! The double-counting xc-contribution is also calculated here
! Formula:
! <i,mu|V_xc(n)|i,nu> =
! = V_xc(n^a)<i,mu|i,nu> + V'_xc(n^a)[<i,mu|n|i,nu> - n^a<i,mu|i,nu>] 
! where
!  V'_xc = d(V_xc)/dn
!  n^a ... average density
!  n_i ... density on i-site
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
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine build_ca_snxc_on (in1, iatom, bcxcx, xc)
        use charges
        use density
        use interactions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: in1
        integer, intent (in) :: iatom
  
! Output
        real, intent (out) :: xc

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
        real q_mu
  
        real, dimension (nsh_max,nsh_max) :: arho
        real, dimension (numorb_max, numorb_max) :: denx

! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Initialize
        xc = 0.0d0
        bcxcx = 0.0d0

! Restore average density for input atom into smaller array.
        do inu = 1, nssh(in1)
         do imu = 1, nssh(in1)
          arho(imu,inu) = arho_on(imu,inu,iatom)
         end do
        end do

! Restore true density for input atom into smaller array.
        do inu = 1, num_orb(in1)
         do imu = 1, num_orb(in1)
          denx(imu,inu) = rho_on(imu,inu,iatom) 
         end do
        end do

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! DIAGONAL TERMS
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Find average density for each shell, and
! Calculate the 1st order approx. for the vxc-matrix elements
! (and 2nd order corrections for the "diagonal" terms)
! This part can be change depending on the used definition of the average density 
        n1 = 0
        do issh = 1, nssh(in1)
         l1 = lssh(issh,in1)
         n1 = n1 + l1 + 1

         ! evaluate V_xc for given density
         call cepal (arho(issh,issh), exc, muxc, dexc, d2exc,dmuxc, d2muxc)

! Now calculate different parts of the XC-matrix 
! Formula:
!  <i,mu|V_xc(n)|i,nu> 
!      = V_xc(n^a)<i,mu|i,nu> + V'_xc(n^a)[<i,mu|n|i,nu> - n^a<i,mu|i,nu>]
! where
!        V_xc(n^a)    ... muxc
!        V'_xc(n^a)   ... dmuxc
!
!        n^a          ... density
!
!        mu /= nu :   <mu|nu> = 0 (off-diagonal term) 
!        mu = nu  :   <mu|nu> = 1 (diagonal term)

! Loop over each orbital
         do ind1 = -l1, l1
          imu = n1 + ind1
              
          ! V_xc(n^a)<i,mu|i,nu>
          bcxcx(imu,imu) = muxc

          ! V'_xc(n^a)[<i,mu|n|i,nu> - n^a<i,mu|i,nu>] 
          bcxcx(imu,imu) = bcxcx(imu,imu)                                    &
     &                    + dmuxc*(denx(imu,imu) -  arho(issh,issh))
         end do

! Contribution to the double-counting correction
         ! average charge per shell
         q_mu = Qin(issh,iatom) / (2*l1 + 1)

! Loop over each orbital
         do ind1 = -l1, l1
          imu = n1 + ind1
          xc = xc + q_mu*(exc - muxc + (dexc - dmuxc)                        &
     &                  *(denx(imu,imu)-arho(issh,issh)))
         end do
         n1 = n1 + l1
        end do ! end 'do issh = 1, nssh(in1)'


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! OFF-DIAGONAL TERMS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Now, find average density for each pair of shells, and
! Calculate the vxc-matrix elements for the "non-diagonal" terms
        n1 = 0

! Loop over shells i x-dimension
        do issh = 1, nssh(in1)

! Number of orbitals per the shell in x-dim
         l1 = lssh(issh,in1)
         n1 = n1 + l1 + 1
         n2 = 0

! Loop over shells i y-dimension
         do jssh = 1, nssh(in1)
          ! number of orbitals per the shell in y-dim
          l2 = lssh(jssh,in1)
          n2 = n2 + l2 + 1
          ! calculate XC potentials 
          call cepal (arho(issh,jssh), exc, muxc, dexc, d2exc, dmuxc, d2muxc)
                 
! Set the XC-submatrices
          ! loop over orbitals in the x-shell
          do ind1 = -l1, l1
           imu = n1 + ind1
           ! loop over orbitals in the y-shell
           do ind2 = -l2, l2
            inu = n2 + ind2
            ! only off-diagonal terms !!!
            if (imu .ne. inu) then
             ! V_xc(n^a)<i,mu|i,nu> = 0
             ! V'_xc(n^a) n^a <i,mu|i,nu> = 0
             ! V_xc(n_i^a)<i,mu|i,nu> = 0
             ! V'_xc(n_i) n_i^a <i,mu|i,nu> = 0
             ! V'_xc(n^a) <i,mu|n|i,nu>
             bcxcx(imu,inu) =  dmuxc*denx(imu,inu) 
            end if
           end do !do ind2 = -l2, l2
          end do !do ind1 = -l1, l1
          n2 = n2 + l2
         end do !do jssh = 1, nssh(in1)
         n1 = n1 + l1
        end do !do issh = 1, nssh(in1)

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================

        return
        end subroutine build_ca_snxc_on
