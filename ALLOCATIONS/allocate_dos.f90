! copyright info:
!
!                             @Copyright 2006
!                           Fireball Committee
! West Virginia University - James P. Lewis, Chair
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


! allocate_dos.f90
! Program Description
! ===========================================================================
! The subroutine allocates the DOS variables. Those variables we need for 
! the dos, hopping or atom calculation (except hamk that it is included in 
! allocate_h).
!
! ===========================================================================
! Code written by:
! Cesar Gonzalez Pascual
! Departamento de Fisica Teorica de la Materia Condensada
! Facultad de Ciencias. Universidad Autonoma de Madrid
! E28049 Madird SPAIN    
! FAX +34-91 3974950
! Office Telephone +34-91 3978648
! ==========================================================================
!
! Program Declaration
! ===========================================================================
      subroutine allocate_dos (natoms, iwrtdos, iwrthop)
      use module_dos
      use neighbor_map
      use interactions
      implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
      integer, intent(in)    :: natoms
      integer, intent(in)    :: iwrtdos
      integer, intent(in)    :: iwrthop

! Procedure
! ===========================================================================
! allocations for the dos calculation
      if (iwrtdos .ge. 1) then 
       allocate (green(norb_act, norb_act, nener))
       green = 0.0d0
      end if 

      return
      end subroutine allocate_dos
