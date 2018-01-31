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

 
! reallocate_h.f90
! Program Description
! ===========================================================================
!       This routine reallocates the arrays which store the interactions
! of the Hamiltonian matrix.  These arrays need to be reallocated if the 
! maximum number of neighbors changes.
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
        subroutine reallocate_h (natoms, neigh_max, neighPP_max, itheory, &
     &                           itheory_xc, igauss)
        use interactions
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
        deallocate (h_mat)
        deallocate (s_mat)
        deallocate (sVNL)
        deallocate (t_mat)
        deallocate (vna)
        deallocate (vnl)
        deallocate (vxc)
        deallocate (vxc_1c)

        allocate (h_mat (numorb_max, numorb_max, neigh_max, natoms))
        allocate (s_mat (numorb_max, numorb_max, neigh_max, natoms))
        allocate (t_mat (numorb_max, numorb_max, neigh_max, natoms))


! PP part
        allocate (sVNL (numorb_max, numorb_max, neighPP_max, natoms))
        allocate (vnl (numorb_max, numorb_max, neighPP_max**2, natoms))

        allocate (vna (numorb_max, numorb_max, neigh_max, natoms))
        allocate (vxc (numorb_max, numorb_max, neigh_max, natoms))
        allocate (vxc_1c (numorb_max, numorb_max, neigh_max, natoms))

        if (igauss .eq. 1) then
         deallocate (bar_density_2c)
         deallocate (bar_density_3c)
         deallocate (density_2c)
         deallocate (density_3c)
         deallocate (nuxc_3c)
         deallocate (nuxc_total)
         deallocate (vxc_3c)

         allocate (bar_density_2c (numorb_max, numorb_max, neigh_max, natoms))
         allocate (bar_density_3c (numorb_max, numorb_max, neigh_max, natoms))
         allocate (density_2c (numorb_max, numorb_max, neigh_max, natoms))
         allocate (density_3c (numorb_max, numorb_max, neigh_max, natoms))
         allocate (nuxc_3c (numorb_max, numorb_max, neigh_max, natoms))
         allocate (nuxc_total (numorb_max, numorb_max, neigh_max, natoms))
         allocate (vxc_3c (numorb_max, numorb_max, neigh_max, natoms))
        end if

        if (itheory .ne. 0) then
         deallocate (ewaldlr)
         deallocate (ewaldsr)
         deallocate (vca)
         deallocate (vxc_ca)
         deallocate (ewaldqmmm)

         allocate (ewaldlr (numorb_max, numorb_max, neigh_max, natoms))
         allocate (ewaldsr (numorb_max, numorb_max, neigh_max, natoms))
         allocate (vca (numorb_max, numorb_max, neigh_max, natoms))
         allocate (vxc_ca (numorb_max, numorb_max, neigh_max, natoms))
         allocate (ewaldqmmm (numorb_max, numorb_max, neigh_max,natoms))
        end if

        if (itheory .eq. 1) then
         deallocate (dip)
!JIMM
         deallocate (dipcm)
         deallocate (dipc)

         allocate (dip (numorb_max, numorb_max, neigh_max, natoms))
         allocate (dipcm (3, numorb_max, numorb_max))
         allocate (dipc (3, numorb_max, numorb_max, neigh_max, natoms))
        end if
 
! Interactions needed for Sankey-Niklewski type average densities.
        if (itheory_xc .eq. 1 .or. itheory_xc .eq. 2) then
         deallocate(sm_mat)
         allocate (sm_mat (nsh_max, nsh_max, neigh_max, natoms))
        end if

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
      end subroutine reallocate_h
