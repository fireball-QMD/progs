! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
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

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
! initAtomicE.f90
! Program Description
! ===========================================================================
!       This calculates the atomic energy.
!
! ===========================================================================
! Code written by:
! Kurt R. Glaesemann
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine initatomicE (natoms, etotatom, imass, atomic_energy)
        use dimensions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: natoms

        integer, intent(in), dimension (natoms) :: imass

        real, intent(in), dimension (nspec_max) :: etotatom
 
! Output
        real, intent(out) :: atomic_energy
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer in1
 
! Procedure
! ===========================================================================
! Determine the atomic energy for the system. This way we can get cohesive
! energies more easily. 
        atomic_energy = 0.0d0
        do iatom = 1, natoms
         in1 = imass(iatom)
         atomic_energy = atomic_energy + etotatom(in1)
        end do
 
! Format Statements
! ===========================================================================
        return
        end
