! copyright info:
!
!                             @Copyright 2003
!                            Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega

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
 
! initamat.f90
! Program Description
! ===========================================================================
! Sets up the matrix amat, which is used by D-orbitals in twister(d)
! If we do not have any D-orbitals, then just return
!
! ===========================================================================
! Code written by:
! Kurt R. Glaesemann
! LLNL
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine initamat(nspecies)
        use constants_fireball
        use interactions
        implicit none

! Argument Declaration and Description
! ===========================================================================
        integer, intent(in) :: nspecies
 
! Local Variable Declaration and Description
! ===========================================================================
        integer in1
        integer issh

! Procedure
! ===========================================================================
! Initialize the a-matrices:  
        haveDorbitals = .false.
        do in1 = 1, nspecies
         do issh = 1, nssh(in1)
          if (lssh(issh,in1) .eq. 2) haveDorbitals = .true.
         end do
         do issh = 1, nsshPP(in1)
          if (lsshPP(issh,in1) .eq. 2) haveDorbitals = .true.
         end do
        end do

        amat(:,:,:) = 0.0d0

        if (.not. haveDorbitals) return

        amat(2,1,1) = 0.5d0
        amat(1,2,1) = 0.5d0
  
        amat(3,2,2) = 0.5d0
        amat(2,3,2) = 0.5d0
   
        amat(1,1,3) = - 1.0d0/(2.0d0*sqrt(3.0d0))
        amat(2,2,3) = - 1.0d0/(2.0d0*sqrt(3.0d0))
        amat(3,3,3) = 2.0d0/(2.0d0*sqrt(3.0d0))
   
        amat(3,1,4) = 0.5d0
        amat(1,3,4) = 0.5d0
 
        amat(1,1,5) = 0.5d0 
        amat(2,2,5) = - 0.5d0
 
        return
        end
