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

 
! allocate_h.f90
! Program Description
! ===========================================================================
!       This routine allocates the arrays which store the interactions
! of the Hamiltonian matrix.
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
        subroutine allocate_h (natoms, neigh_max, neighPP_max, itheory,      &
     &                         itheory_xc, igauss, iwrtdos, iwrthop, iwrtatom)

        use interactions
        use options, only : idipole, V_intra_dip, iks
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: igauss
        integer, intent (in) :: itheory
        integer, intent (in) :: itheory_xc
        integer, intent (in) :: natoms
        integer, intent (in) :: neigh_max
        integer, intent (in) :: neighPP_max
 ! CGP
        integer, intent (in) :: iwrtdos
        integer, intent (in) :: iwrthop
        integer, intent (in) :: iwrtatom
!end CGP

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        allocate (h_mat (numorb_max, numorb_max, neigh_max, natoms))
        allocate (s_mat (numorb_max, numorb_max, neigh_max, natoms))
        allocate (t_mat (numorb_max, numorb_max, neigh_max, natoms))



! PP part
        allocate (sVNL (numorb_max, numorb_max, neighPP_max, natoms))
        allocate (vnl (numorb_max, numorb_max, neighPP_max**2, natoms))

        allocate (vna (numorb_max, numorb_max, neigh_max, natoms))
        allocate (vxc (numorb_max, numorb_max, neigh_max, natoms))
        allocate (vxc_1c (numorb_max, numorb_max, neigh_max, natoms))

! Interactions needed for gaussian approximation to three-center
! exchange-correlation interactions.        
        if (igauss .eq. 1) then
         allocate (bar_density_2c (numorb_max, numorb_max, neigh_max, natoms))
         allocate (bar_density_3c (numorb_max, numorb_max, neigh_max, natoms))
         allocate (density_2c (numorb_max, numorb_max, neigh_max, natoms))
         allocate (density_3c (numorb_max, numorb_max, neigh_max, natoms))
         allocate (nuxc_3c (numorb_max, numorb_max, neigh_max, natoms))
         allocate (nuxc_total (numorb_max, numorb_max, neigh_max, natoms))
         allocate (vxc_3c (numorb_max, numorb_max, neigh_max, natoms))
        end if 

! Interactions needed for SCF algorithms - either DOGS or extended-Hubbard
        if (itheory .ne. 0) then
         allocate (ewald (natoms, natoms))
         allocate (ewaldlr (numorb_max, numorb_max, neigh_max, natoms))
         allocate (ewaldsr (numorb_max, numorb_max, neigh_max, natoms))
         allocate (vca (numorb_max, numorb_max, neigh_max, natoms))
         allocate (vxc_ca (numorb_max, numorb_max, neigh_max, natoms))
         allocate (ewaldqmmm (numorb_max, numorb_max, neigh_max,natoms))
        end if

! zw mcweda second order
        if (itheory_xc .eq. 4) then
         allocate (g2nu(nsh_max,nsh_max,neigh_max,natoms))
                              
         allocate (g2nup(3,nsh_max,nsh_max,neigh_max,natoms))
        end if !end if itheory_xc .eq. 4

! Interactions needed only for DOGS
        if (itheory .eq. 1 .or. idipole .eq. 1 .or. iks .eq. 1) then
         allocate (dip (numorb_max, numorb_max, neigh_max, natoms))
! JIMM
         allocate (dipcm (3, numorb_max, numorb_max))
         allocate (dipc (3, numorb_max, numorb_max, neigh_max, natoms))
        endif
!Intra-atomic dipolar potential
        if (V_intra_dip .eq. 1) then
          allocate(Vdip_1c(numorb_max,numorb_max,natoms))
        end if
! Interactions needed only for extended-Hubbard
        if (itheory .eq. 2) then
         allocate (Vcoulomb (nsh_max, natoms))
         allocate (Vewaldsr (nsh_max, natoms))
         allocate (Vxcnu (nsh_max, natoms))
        end if
 
! Interactions needed for Sankey-Niklewski type average densities.
        if (itheory_xc .eq. 1 .or. itheory_xc .eq. 2 .or. itheory_xc .eq. 4) then
         allocate (sm_mat (nsh_max, nsh_max, neigh_max, natoms))
        end if

! CGP
! allocations for the dos calculation
        if (iwrtdos.ge.1.or.iwrthop.ge.1.or.iwrtatom.ge.1) then
           allocate (hamk(norbitals,norbitals))
           hamk = 0.0d0
        end if
!end CGP

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
      end subroutine allocate_h
