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

! allocate_rho.f90
! Program Description
! ===========================================================================
!       This routine allocates the arrays used within TDSE module.
!
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine allocate_tdse ()

        use density
        use interactions
        use configuration
        use kpoints
        use options
        use tdse
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================

! allocate array of eigenvalues and eigenfunctions
        if (iqout .ne. 2) allocate (blowre (norbitals, norbitals, nkpoints))
        if (iqout .ne. 2 .and. icluster .ne. 1)                           &
     &      allocate (blowim (norbitals, norbitals, nkpoints))
        allocate (bbnkre (norbitals, norbitals, nkpoints))
        if (icluster .ne. 1)                                              &
     &      allocate (bbnkim (norbitals, norbitals, nkpoints))
        allocate (eigen_k (norbitals, nkpoints))

        allocate (Hlow (norbitals, norbitals, nkpoints))    ! orthogonal Hamiltonian in MO basis set
        allocate (psi (norbitals, nelec, nkpoints))         ! TD-wf in MO basis set
        allocate (psiAO (norbitals, nelec, nkpoints))       ! TD-wf in AO basis set
        allocate (sm12 (norbitals, norbitals, nkpoints))    ! S^(-1/2)
        allocate (Enev (nelec, nkpoints))                   ! Energy expectation value of electron
        allocate (psi2es (norbitals, nelec, nkpoints))      ! projection psi to eigenfuction

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================

        return
        end subroutine allocate_tdse
