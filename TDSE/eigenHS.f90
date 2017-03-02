! copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
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

! eigenHS.f90
! Program Description
! ===========================================================================
!       This routine gets noSCF eigenvalues & eigenvectors of H,S
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine eigenHS (itime_step)

        use options
        use tdse

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itime_step

! Local Parameters and Data Declaration
! ===========================================================================


! Local Variable Declaration and Description
! ===========================================================================
        integer iclock
        logical isH

! Procedure
! ===========================================================================
! do we dump psi in this time step?
        iclock = mod(itime_step,np2es)
        write (*,*) 'iclk',itime_step,np2es,iclock
        if ((iclock .eq. 0) .or. (itime_step .eq. 1)) then
          isH = .true.
          write (*,*) ' In this time step we project actual psi onto eigenstate;'
          write (*,*) ' full diagonalization will be performed.'
        else
          isH = .false.
        endif

! assemble H,S
        if (isH) then

! assemble Hamiltonian
         call assemble_h ()
! diagonalize H,S; get S^(1/2),S^(-1/2)
         write (*,*)  ' Call tddiag_k'
         call tddiag_k (isH)

! project psi onto set of eigenstates
         write (*,*) 'call psi2es'
         call get_psi2es (itime_step)

! get charges and rho (??)
!         call build_rho (0)

        else

! assemble S
         call assemble_h ()
! diagonalize S; get S^(1/2),S^(-1/2) and orthogonalize H
         write (*,*)  ' Call tddiag_k'
         call tddiag_k (isH)

        endif

! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine eigenHS

