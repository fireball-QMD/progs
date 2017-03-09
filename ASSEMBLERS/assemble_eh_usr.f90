! copyright info:
!
!                             @Copyright 2006
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold
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

 
! assemble_eh_usr.f90
! Program Description
! ===========================================================================
!       This is a subroutine for calculating energy and energy correction    &
! terms for extended hubbard.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658 
! Provo, UT 841184602-4658
! FAX 801-422-2265
! Office telephone 801-422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_eh_usr (iforce, ioff2c, ehxcc, ehcoolc)
        use charges
        use configuration
        use constants_fireball
        use dimensions
        use integrals
        use interactions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iforce

        integer, intent (in), dimension (1:24) :: ioff2c

! Ouput
        real, intent (out) :: ehxcc
        real, intent (out) :: ehcoolc
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ideriv
        integer in1, in2
        integer index
        integer index_coulomb
        integer ineigh
        integer interaction
        integer issh
        integer jatom
        integer jssh
        integer matom
        integer mbeta
        integer n1, n2

        real distance
        real dq1

        real, dimension (nsh_max, nsh_max) :: coulomb
        real, dimension (nsh_max, nsh_max) :: coulombD
        real, dimension (ME2c_max) :: dslist
        real, dimension (3) :: r1, r2, r21
        real, dimension (ME2c_max) :: slist
        real, dimension (natoms) :: sub_ewald
        real, dimension (nsh_max, nsh_max) :: xcnu
        real, dimension (nsh_max, nsh_max) :: xcnuD
 
! Procedure
! ===========================================================================
! Initialize 
        write (*,*) '  '
        write (*,*) ' Welcome to assemble_eh_usr.f! '
        write (*,*) '  '
        
        ehxcc = 0.0d0
        ehcoolc = 0.0d0

        do iatom = 1, natoms
         in1 = imass(iatom)
         do issh = 1, nssh(in1)
          ehxcc = ehxcc - Vxcnu(issh,iatom)*Qin(issh,iatom)
          ehcoolc = ehcoolc - eq2*Vcoulomb(issh,iatom)*Qin(issh,iatom)
         end do
        end do
 
! Now we add the 1/2 nu_ij dq1_i dq1_j term.
! First calculate V_alpha using Qin (which is now the mixed Q).  (Note I 
! intentionally computed the q_i V_i term using V_i which was generated with 
! Qin before mixing... Otherwise q_i V_i will not exactly equal compHubb 
! unless the charge has converged perfectly).
 
! The following is just from assemble_eh_2c.f.
 
! (See Eq. 15 of hubbard.ps)
! V_alpha = sum_beta Interaction_alpha,beta * dq1_beta
! There are two types of interactions. Coulomb and xcnu.
! Note that alpha, beta means atomshells.
 

! First a preliminary quantity. ewald(i,j) = sum(L) 1/(| bi-bj +/- L|
! We need to calculate SUMS of ewald sums in which the charge is included.
! These will also be used later in Pnanl3pP.f.  The details are on page 10 of
! the "Derivative of the long range matrix element". In the nutshell we are
! calculating:
!              Sum_(i,L) q(i)/|b(i)-b(alpha)+L|

        do iatom = 1, natoms
         sub_ewald(iatom) = 0.0d0
         do jatom = 1, natoms
          in2 = imass(jatom)

! Calculate the charge on jatom
          dq1 = 0.0d0
          do issh = 1, nssh(in2)
           dq1 = dq1 + (Qin(issh,jatom) - Qneutral(issh,in2))
          end do
          sub_ewald(iatom) = sub_ewald(iatom) + dq1*ewald(iatom,jatom)
         end do
        end do

! Initialize Vcoulomb, Vxcnu
        Vcoulomb = 0.0d0
        Vewaldsr = 0.0d0
        Vxcnu = 0.0d0

! ****************************************************************************
! Loop over the atoms in the central cell.
        do iatom = 1, natoms                          ! <==== loop over iatom
         matom = neigh_self(iatom)
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)
 
! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)     ! <==== loop over iatom's neighbors
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)

! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
          r21(:) = r2(:) - r1(:)
          distance = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))


! ***************************************************************************
!
! GET EXTENDED-HUBBARD COULOMB INTERACTIONS - STORE IN VCOULOMB.
! ***************************************************************************
! Now find the coulomb integrals need to evaluate the extended-hubbard
! interaction.
! Loop over all the non-zero integrals for this interaction:
          if (ioff2c(24) .eq. 1) then 
           index_coulomb = nssh(in1)*nssh(in2)
           interaction = 12
           ideriv = 0
           do index = 1, index_coulomb
            call interpolate_1d (interaction, ideriv, in1, in2, index,      &
     &                           iforce, distance, slist(index), dslist(index))
           end do

! We have the data, it is stored in the following way: v(1,1), v(1,2),
! v(1,3)...v(1,n2max), v(2,1), v(2,2), ..., but it is a 1D array,
! ordered according to the index_coulomb. Restore this to true 2x2 format:
           n1 = nssh(in1)
           n2 = nssh(in2)
           call recoverC (n1, n2, slist, dslist, coulomb, coulombD)

           do issh = 1, nssh(in1)
            do jssh = 1, nssh(in2)
             Vcoulomb(issh,iatom) = Vcoulomb(issh,iatom)                    &
     &        +  coulomb(issh,jssh)*(Qin(jssh,jatom) - Qneutral(jssh,in2))
            end do
           end do
          end if

! ****************************************************************************
!
! GET EXTENDED-HUBBARD SHORT-RANGE EWALD INTERACTIONS - STORE IN EWALDSR
! ****************************************************************************
          if (matom .ne. ineigh .and. ioff2c(22) .eq. 1) then
           do issh = 1, nssh(in1)
            do jssh = 1, nssh(in2)
             Vewaldsr(issh,iatom) = Vewaldsr(issh,iatom)                     &
     &        + (Qin(jssh,jatom) - Qneutral(jssh,in2))/distance
            end do
           end do
          end if


! ****************************************************************************
!
! GET EXTENDED-HUBBARD EXCHANGE-CORRELATION INTERACTIONS - STORE IN VXCNU
! ****************************************************************************
! Now find the exchange-correlation integrals need to evaluate the
! extended-hubbard interaction.
! Loop over all the non-zero integrals for this interaction:
          interaction = 14
          ideriv = 0
          do index = 1, index_coulomb
           call interpolate_1d (interaction, ideriv, in1, in2, index,  &
     &                          iforce, distance, slist(index), dslist(index))
          end do

! We have the data, it is stored in the following way: v(1,1), v(1,2),
! v(1,3)...v(1,n2max), v(2,1), v(2,2), ..., but it is a 1D array,
! ordered according to the index_coulomb. Restore this to true 2x2 format:
          call recoverC (n1, n2, slist, dslist, xcnu, xcnuD)

          if (matom .ne. ineigh) then
           do issh = 1, nssh(in1)
            do jssh = 1, nssh(in2)
             Vxcnu(issh,iatom) = Vxcnu(issh,iatom)                           &
     &        + xcnu(issh,jssh)*(Qin(jssh,jatom) - Qneutral(jssh,in2))
            end do
           end do
          else
           do issh = 1, nssh(in1)
            do jssh = 1, nssh(in2)
             Vxcnu(issh,iatom) = Vxcnu(issh,iatom)                           &
     &        + xcnu1c(issh,jssh,in1)*(Qin(jssh,jatom) - Qneutral(jssh,in2))
            end do
           end do
          end if


! ****************************************************************************
! End loop over iatom and its neighbors - jatom.
         end do
        end do

! Subtract the short-range and add the long-range electrostatics into Vcoulomb 
! Calculate ehxcc and ehcoolc => 1/2 V_alpha *dq1_i*dq1_j
        do iatom = 1, natoms
         in1 = imass(iatom)
         do issh = 1, nssh(in1)
          ehxcc = ehxcc                                                      &
     &     + Vxcnu(issh,iatom)*(Qin(issh,iatom) - Qneutral(issh,in1))/2.0d0
          ehcoolc = ehcoolc                                                  &
     &     + (eq2/2.0d0)*(Vcoulomb(issh,iatom) - Vewaldsr(issh,iatom)        &
     &                    + sub_ewald(iatom))                                &
     &                  *(Qin(issh,iatom) - Qneutral(issh,in1))
         end do
        end do 

! ===========================================================================
        return
        end
