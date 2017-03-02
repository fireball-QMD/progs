! copyright info:
!
!                             @Copyright 2009
!                FAST (Fireball Atomic Simulation Techniques)
! West Virginia University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
! Institute of Physics, Czech Republic - Pavel Jelinek

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

! save_mdetstuff.f90
! Program Description
! ===========================================================================
!       This routine saves some variables needed for the time
!       integration of the TD-wfs
!
! ===========================================================================
! Code written by Jose Ortega Mateo
! Departmento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine save_mdetstuff ()

        use configuration
        use nonadiabatic
        use density

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================


! Procedure
! ===========================================================================
        ratom_old = ratom
        vatom_old = vatom
        gks_old = gks
        eigen_old = eigen_k
        bbnkre_old = bbnkre 
        blowre_old = blowre 

! Deallocate Arrays
! ===========================================================================
!       deallocate (bbnkre)
!       deallocate (blowre)
!       deallocate (eigen_k)


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine save_mdetstuff

