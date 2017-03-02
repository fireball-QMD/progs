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

! tclmdtransfer.f90
! Program Description
! ===========================================================================
!       TclMD utility for FIREBALL data transfer to VMD; Molecule number 
! (molNumber) is not used.
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
        subroutine tclmdatomtransfer (natoms)
        use dimensions
        use configuration
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
 
! Input
        integer, intent(in) :: natoms
 
! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer molNumber
        integer natoms_tclmd
 
        real*8 x, y, z

! Procedure
! ===========================================================================
! Body - calls a c function
! "tclmdsendmol_(*int,*int,*double,*double,*double)" that should be provided
        molNumber = 0
        natoms_tclmd = natoms
        do iatom = 1, natoms
         x = dble(ratom(1,iatom) + ximage(1,iatom))
         y = dble(ratom(2,iatom) + ximage(2,iatom))
         z = dble(ratom(3,iatom) + ximage(3,iatom))
 
         call tclmdsendmol(molNumber, natoms_tclmd, x,y,z)
 
!end loop to step through all atoms of molecule
        end do
 
! Format Statements
! ===========================================================================

        return
        end
