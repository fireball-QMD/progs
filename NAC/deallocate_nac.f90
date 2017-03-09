! copyright info:
!
!                             @Copyright 2006
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
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

 
! deallocate_nac.f90
! Program Description
! ===========================================================================
!       This routine dallocates the arrays which store the interactions
! for the calculation of the non-adiabatic couplings.  
! These arrays need to be reallocated if the 
! maximum number of neighbors changes.
!
! ===========================================================================
! Code written by:
! Jose Ortega Mateo
! Departmento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine deallocate_nac (natoms)
!
        use nonadiabatic
        use interactions
        use neighbor_map
        use options
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        deallocate (gover)
        deallocate (gover1c)
        deallocate (gh_2c)
        deallocate (gh_atm)
        deallocate (gh_3c)
        deallocate (gh_lrew_qmmm)
! merge large arrays VLADA
!     deallocate (gh_xc_3c)
! JOM-q or norbitals_new
!       deallocate (gks)


! PP part
!       allocate (gh_pp_2c (3, numorb_max, numorb_max, neighPP_max, natoms))
        deallocate (gh_pp_otr)
        deallocate (gh_pp_otl)
        deallocate (gh_pp_atm)
        deallocate (gh_pp_3c)

!        if (itheory .eq. 1) then
!                deallocate (gh_2c_ca)
!                deallocate (gh_atm_ca)
!                deallocate (gh_3c_ca)
!                deallocate (gh_lrew)
!        end if

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
      end subroutine deallocate_nac
