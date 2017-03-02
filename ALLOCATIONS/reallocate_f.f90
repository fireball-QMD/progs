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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
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
        if (itheory .eq. 1) then
         deallocate (dipp)

         allocate (dipp (3, numorb_max, numorb_max, neigh_max, natoms))

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
        if (itheory_xc .eq. 1 .or. itheory_xc .eq. 2) then
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
