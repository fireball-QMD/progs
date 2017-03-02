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

! writeout_xv.f90
! Program Description
! ===========================================================================
!       This routine writes out the positions and velocities to the xvfile.
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
        subroutine writeout_xv (xvfile, itime_step, time, imass, nzx, iwrtvel)
        use dimensions
        use configuration
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itime_step
        integer, intent (in) :: iwrtvel
 
        integer, intent (in), dimension (natoms) :: imass
        integer, intent (in), dimension (nspecies) :: nzx
        
        real, intent (in) :: time

        character (len=30), intent (in) :: xvfile
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer in1
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        write (*,*) '  '
        write (*,*) ' Writing positions and velocities to the xvfile: '
! xv-file
        open (unit = 63, file = xvfile, status = 'unknown')
        write (63,100) natoms
        write (63,101) itime_step, time
        do iatom = 1, natoms
         in1 = imass(iatom)
         write (63,102) nzx(in1), ratom(:,iatom) + ximage(:,iatom),     &
     &                  vatom(:,iatom)
        end do
        close (unit=63)

! velocity file
        if (iwrtvel .eq. 1) then
         open (unit = 62, file = "VELOCITY.dat", status = 'unknown',     &
     &        position = 'append')
         if (itime_step .eq. 1) write (62,100) natoms
         write (62,101) itime_step, time
         do iatom = 1, natoms
          in1 = imass(iatom)
          write (62,102) nzx(in1), ratom(:,iatom) + ximage(:,iatom),     &
     &                  vatom(:,iatom)
         end do
! close velocity.dat file
         close (unit=62)
 
        endif

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (2x, i4)
101     format (2x, i6, 2x, f9.3)
102     format (2x, i2, 3(2x,f16.9), 3(2x,f16.9))
 
        close (unit = 63)
        return
      end subroutine writeout_xv
