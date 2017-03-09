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


! allocate_f.f90
! Program Description
! ===========================================================================
!       This routine allocates the arrays which store the derivatives of the
! interactions of the Hamiltonian matrix.
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
        subroutine allocate_f (natoms, neigh_max, neighPP_max, numorb_max,  &
     &                         nsh_max, itheory, itheory_xc, igauss, ivdw,  &
     &                         iharmonic, ibias)
        use forces
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: igauss
        integer, intent (in) :: iharmonic
        integer, intent (in) :: itheory
        integer, intent (in) :: itheory_xc
        integer, intent (in) :: ivdw
        integer, intent (in) :: ibias
        integer, intent (in) :: natoms
        integer, intent (in) :: neigh_max
        integer, intent (in) :: neighPP_max
        integer, intent (in) :: numorb_max
        integer, intent (in) :: nsh_max
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
 
! Allocate Arrays
! ===========================================================================
! Allocate derivatives of interactions
        allocate (sp_mat (3, numorb_max, numorb_max, neigh_max, natoms))
        allocate (spVNL (3, numorb_max, numorb_max, neighPP_max, natoms))
        allocate (tp_mat (3, numorb_max, numorb_max, neigh_max, natoms))

! Allocate components of the forces
        allocate (dusr (3, natoms))
        allocate (dxcv (3, natoms))
        allocate (fro (3, natoms))
        allocate (ft (3, natoms))
        allocate (ftot (3, natoms))
        allocate (ftotold (3, natoms))
        allocate (ftotnew (3, natoms))
        allocate (fana (3, neigh_max, natoms))
        allocate (faxc (3, neigh_max, natoms))
        allocate (f3naa (3, natoms))
        allocate (f3nab (3, natoms))
        allocate (f3nac (3, natoms))
        allocate (f3nla (3, natoms))
        allocate (f3nlb (3, natoms))
        allocate (f3nlc (3, natoms))
        allocate (f3xca (3, natoms))
        allocate (f3xcb (3, natoms))
        allocate (f3xcc (3, natoms))
        allocate (fotxc (3, neigh_max, natoms))
        allocate (fotna (3, neigh_max, natoms))
 
! Procedure
! ===========================================================================
! NOTE: here we can allocate matrix of size only neighPP_max, because 
! ontol and atomic cases are just  2 center interaction, but 
! pure 3 center interactions is different story 
        allocate (fanl (3, neighPP_max, natoms))
        allocate (fotnl (3, neighPP_max, natoms))

        if (igauss .eq. 1) allocate (fxcro (3, neigh_max, natoms))

! Allocate components of the forces - needed for either DOGS or extended-Hubbard
        if (itheory .ne. 0) then
         allocate (dewald (3, natoms, natoms))
         allocate (fewald (3, natoms))
         allocate (flrew (3, natoms))
         allocate (flrew_qmmm (3, natoms))
        end if

! Allocate components of the forces - needed for DOGS   
        if (itheory .eq. 1) then
         allocate (dipp (3, numorb_max, numorb_max, neigh_max, natoms))

         allocate (faca (3, neigh_max, natoms))
         allocate (faxc_ca (3, neigh_max, natoms))
         allocate (f3caa (3, natoms))
         allocate (f3cab (3, natoms))
         allocate (f3cac (3, natoms))
         allocate (f3xca_ca (3, natoms))
         allocate (f3xcb_ca (3, natoms))
         allocate (f3xcc_ca (3, natoms))
         allocate (fotca (3, neigh_max, natoms))
         allocate (fotxc_ca (3, neigh_max, natoms))
        end if

! Allocate components of the forces - needed for extended Hubbard
        if (itheory .eq. 2) then
         allocate (fcoulomb (3, natoms))
         allocate (fxcnu (3, natoms))
        end if

! Allocate snxc forces
        if (itheory_xc .eq. 1 .or. itheory_xc .eq. 2) then
         allocate (spm_mat (3, nsh_max, nsh_max, neigh_max, natoms))
         allocate (arhop_off (3, nsh_max, nsh_max, neigh_max, natoms)) 
         allocate (arhopij_off (3, nsh_max, nsh_max, neigh_max, natoms)) 
         allocate (rhop_off (3, numorb_max, numorb_max, neigh_max, natoms)) 
         allocate (rhopij_off (3, numorb_max, numorb_max, neigh_max, natoms)) 
         allocate (arhop_on (3, nsh_max, nsh_max, neigh_max, natoms)) 
         allocate (rhop_on (3, numorb_max, numorb_max, neigh_max, natoms)) 
         ! OLSXC double count corr forces
         allocate (dxcdcc (3, neigh_max, natoms))
        end if

! Allocate vdw forces if requested.
        if (ivdw .eq. 1) allocate (fvdw (3, natoms))

! Allocate external field forces for thermodynamic integration if requested.
        if (iharmonic .eq. 1) allocate (fharmonic (3, natoms)) 
 
! Allocate bias forces for bias voltage field if requested.
        if (ibias .eq. 1) allocate (fbias (3, natoms)) 

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================

        return
        end subroutine allocate_f
