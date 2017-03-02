! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Motorola, Physical Sciences Research Labs - Alex Demkov
! University of Regensburg - Juergen Fritsch
! Universidad de Madrid - Jose Ortega

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola - Jun Wang
! Ohio University - Dave Drabold

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! gaussT.f90
! Program Description
! ===========================================================================
!       This routine renormalizes the forces for each atom, such that 
! a constant temperature in the system is maintained.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine gaussT (natoms, vatom, xmass, T_want, ftot) 
        use constants_fireball
        use dimensions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms

        real, intent (in) :: T_want

        real, intent (in), dimension (3, natoms) :: vatom
        real, intent (in), dimension (natoms) :: xmass

! Output
        real, intent (inout), dimension (3, natoms) :: ftot
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom

        real denominator
        real xksi
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! We are doing gaussian constant T dynamics.
        write (*,*) '  '
        write (*,*) ' We are doing constant temperature dynamics! '

! Calculate xksi (See Evans et al. P.R.A 28, 1016 (1983).
! xksi = SUM(1:3,1:natoms) ftot(ix,iatom)*v(ix,iatom)/(3/2*n*kB*T)
        denominator = 2.0d0*(3.0d0/2.0d0)*natoms*T_want/kconvert
        xksi = 0.0d0
        do iatom = 1, natoms
         xksi = xksi + ftot(1,iatom)*vatom(1,iatom)                          &
     &               + ftot(2,iatom)*vatom(2,iatom)                          &
     &               + ftot(3,iatom)*vatom(3,iatom)
        end do
        if (denominator .lt. 1.0d-3) then
         write (*,*) ' T_want = ', T_want
         write (*,*) ' This value will not work. Reset! '
         stop
        end if
        xksi = xksi/denominator
        write (*,*) ' xksi = ', xksi, ' RENORMALIZE ftot! '

! fovermp is a conversion factor; velocity in A/fs.
        do iatom = 1, natoms
         ftot(:,iatom) =                                                     &
     &    ftot(:,iatom) - xksi*xmass(iatom)*vatom(:,iatom)/fovermp
        end do
        do iatom = 1, natoms
         write (*,100) iatom, ftot(:,iatom)
        end do
        write (*,*) '  '

! Deallocate Arrays
! ===========================================================================
100     format (2x, ' iatom = ', i4, ' ftot   = ', 3e14.6)
 
! Format Statements
! ===========================================================================
 
        return
        end
