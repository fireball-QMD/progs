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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
! allocate_nac.f90
! Program Description
! ===========================================================================
!       This routine allocates the arrays which store the interactions
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
        subroutine allocate_nac (natoms)
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
!       deallocate (gover)
!       deallocate (gover1c)
!       deallocate (gh_2c)
!       deallocate (gh_atm)
!       deallocate (gh_3c)
!       deallocate (gh_pp_2c)
!       deallocate (gh_pp_atm)
!       deallocate (gh_pp_3c)
        if (allocated (gover)) then
          deallocate (gover)
          deallocate (gover1c)
          deallocate (gh_2c)
          deallocate (gh_atm)
          deallocate (gh_3c)
! merge large arrays VLADA 
     !    deallocate (gh_xc_3c)

! PP part
!       allocate (gh_pp_2c (3, numorb_max, numorb_max, neighPP_max, natoms))
          deallocate (gh_pp_otr)
          deallocate (gh_pp_otl)
          deallocate (gh_pp_atm)
          deallocate (gh_pp_3c)
! merge large arrays VLADA
      !   if (itheory .eq. 1) then
      !      deallocate (gh_2c_ca)
      !      deallocate (gh_atm_ca)
      !      deallocate (gh_3c_ca)
      !      deallocate (gh_lrew)
      !   end if
        endif ! if allocated 
        
        if (.not. allocated (gover)) then
          allocate (gover (3, numorb_max, numorb_max, neigh_max, natoms))
          allocate (gover1c (3, numorb_max, numorb_max))
          allocate (gh_2c    (3, numorb_max, numorb_max, neigh_max, natoms))
          allocate (gh_atm   (3, numorb_max, numorb_max, neigh_max, natoms))
          allocate (gh_3c    (3, natoms, numorb_max, numorb_max, neigh_max, natoms))
          allocate (gh_lrew_qmmm (3, natoms, numorb_max, numorb_max, neigh_max,natoms))
! merge large arrays VLADA
       !  allocate (gh_xc_3c (3, natoms, numorb_max, numorb_max, neigh_max, natoms))

! PP part
!       allocate (gh_pp_2c (3, numorb_max, numorb_max, neighPP_max, natoms))
          allocate (gh_pp_otr (3, numorb_max, numorb_max, neighPP_max, natoms))
          allocate (gh_pp_otl (3, numorb_max, numorb_max, neighPP_max, natoms))
          allocate (gh_pp_atm (3, numorb_max, numorb_max, neighPP_max, natoms))
          allocate (gh_pp_3c (3, natoms, numorb_max, numorb_max, neighPP_max**2, natoms))
! merge large arrays VLADA
        !  if (itheory .eq. 1) then
        !    allocate (gh_2c_ca (3, numorb_max, numorb_max, neigh_max, natoms))
        !    allocate (gh_atm_ca (3, numorb_max, numorb_max, neigh_max, natoms))
        !    allocate (gh_3c_ca (3, natoms, numorb_max, numorb_max, neigh_max, natoms))
        !    allocate (gh_lrew  (3, natoms, numorb_max, numorb_max, neigh_max, natoms))
        !  end if
!
        end if ! not. allocate
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
      end subroutine allocate_nac
