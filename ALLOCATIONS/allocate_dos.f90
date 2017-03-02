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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

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
