! copyright info:
!                             @Copyright 2001
!                            Fireball Committee
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
 
! get_exfield.f90
! Program Description
! ===========================================================================
! This routine adds any external constraint energies to the main energy
! such as harmonic oscillators. This is done for every atom in the system.
! ===========================================================================
! Code written by:
! J. B. Keith
! jbrkeith@gmail.com
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine getHarmonic ()
        use configuration
        use constants_fireball
        use forces
        use interactions
!        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
! Output
! Local Parameters and Data Declaration
! ===========================================================================
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom,i
        real, dimension(3) :: R
        real :: omega = 8.751663874677414 !Debye frequency for Titanium in 1/fs
        real :: conversion = 0.000104403692 !convert to from au,fs,Ang to eV
        
! Allocate Arrays
! ===========================================================================
! Procedure
! ===========================================================================

! use the original sites after having been shifted (just so it can be used with xdot)
        enHarmonic = 0.0d0
        fharmonic = 0.0d0
        do iatom = 1, natoms
         R = nowMinusInitialPos(:,iatom)
!         do i=1,3
!          R(i) = original_shift(i,iatom) - ratom(i,iatom)
!         enddo
! calculate energies of current position compared to original site
         enHarmonic = enHarmonic + 0.5*xmass(iatom)*omega*omega*dot_product(R,R)*conversion
! calculate forces of current position compared to original site
         fharmonic(:,iatom) = -xmass(iatom)*omega*omega*R*conversion
!         write(*,*)'xmass',xmass(iatom),'omega',omega
!         write(*,*)'R_original',original_shift(:,iatom),'R_new',ratom(:,iatom)
!         write(*,*)'R',R
!         write(*,*)'for',iatom,'the harmonic osciallator force is',fharmonic(:,iatom)
        end do
!        write(*,*)'the harmonic oscillator energy is',enHarmonic
! Deallocate Arrays
! ===========================================================================
! Format Statements
! ===========================================================================
        return
        end
