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

 
! reallocate_f.f90
! Program Description
! ===========================================================================
!       This routine reallocates the arrays which store the derivatives of 
! the interactions of the Hamiltonian matrix.  These arrays need to be 
! reallocated if the maximum number of neighbors changes.
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
        subroutine reallocate_f (natoms, neigh_max, neighPP_max, itheory,   & 
     &                           itheory_xc, igauss )
        use forces
        use density
        use interactions
        use options, only : idipole
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

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Reallocate derivatives of interactions
        deallocate (sp_mat)
        deallocate (tp_mat)

        allocate (sp_mat (3, numorb_max, numorb_max, neigh_max, natoms))
        allocate (tp_mat (3, numorb_max, numorb_max, neigh_max, natoms))

! Reallocate components of the forces
        deallocate (fana)
        deallocate (faxc)
        deallocate (fotna)

        deallocate (fotxc)
        if (igauss .eq. 1) deallocate (fxcro)
 
        allocate (fana (3, neigh_max, natoms))
        allocate (faxc (3, neigh_max, natoms))
        allocate (fotxc (3, neigh_max, natoms))
        allocate (fotna (3, neigh_max, natoms))


        if (igauss .eq. 1) allocate (fxcro (3, neigh_max, natoms))

! Reallocate components of the forces needed for DOGS
        if (itheory .eq. 1 .or. idipole .eq. 1) then
!JIMM
         deallocate (dipp)
         deallocate (dippcm)
         deallocate (dippc)
         allocate (dipp (3, numorb_max, numorb_max, neigh_max, natoms))
         allocate (dippcm (3, 3, numorb_max, numorb_max))
         allocate (dippc (3, 3, numorb_max, numorb_max, neigh_max, natoms))

         deallocate (faca)
         deallocate (faxc_ca)
         deallocate (fotca)
         deallocate (fotxc_ca)

         allocate (faca (3, neigh_max, natoms))
         allocate (faxc_ca (3, neigh_max, natoms))
         allocate (fotca (3, neigh_max, natoms))
         allocate (fotxc_ca (3, neigh_max, natoms))
        end if

! Allocate snxc forces
        if (itheory_xc .eq. 1 .or. itheory_xc .eq. 2 .or. itheory_xc .eq. 4) then
         deallocate (spm_mat)
         deallocate (arhop_off)
         deallocate (arhopij_off)
         deallocate (rhop_off)
         deallocate (rhopij_off)
         deallocate (arhop_on)
         deallocate (rhop_on)
         deallocate (dxcdcc)
         allocate (spm_mat (3, nsh_max, nsh_max, neigh_max, natoms))
         allocate (arhop_off (3, nsh_max, nsh_max, neigh_max, natoms))
         allocate (arhopij_off (3, nsh_max, nsh_max, neigh_max, natoms))
         allocate (rhop_off (3, numorb_max, numorb_max, neigh_max, natoms))
         allocate (rhopij_off (3, numorb_max, numorb_max, neigh_max, natoms))
         allocate (arhop_on (3, nsh_max, nsh_max, neigh_max, natoms))
         allocate (rhop_on (3, numorb_max, numorb_max, neigh_max, natoms))
         allocate (dxcdcc (3, neigh_max, natoms))
        end if

! Allocate xczw forces (double countig correction)
         if (itheory_xc .eq. 4) then 
             deallocate (dxcdcc_zw)
             allocate (dxcdcc_zw (3, neigh_max,natoms))
         end if
 
! PP part
! Deallocate
        deallocate (spVNL)
        deallocate (fanl)
        deallocate (fotnl)

! Allocate 
        allocate (spVNL (3, numorb_max, numorb_max, neighPP_max, natoms))
        allocate (fotnl (3, neighPP_max, natoms))
        allocate (fanl (3, neighPP_max, natoms))


! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
      end subroutine reallocate_f
