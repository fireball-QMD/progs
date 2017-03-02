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

! Program Description
! ===========================================================================
!       This routine writes out the acceleration to the file - acfile.
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
        subroutine writeout_ac (acfile, itime_step, time, imass, nzx)
        use dimensions
        use configuration
        use constants_fireball
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itime_step

        integer, intent (in), dimension (natoms) :: imass
        integer, intent (in), dimension (nspecies) :: nzx
 
        real, intent (in) :: time

        character (len=30), intent (in) :: acfile
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer in1
        integer iorder
        integer i
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        write (*,*) ' Writing accelerations to the acfile. '
        open (unit = 64, file = acfile, status = 'unknown')
        write (64,100) natoms
        do iorder = 2, 5
         write (64,101) iorder, itime_step, time
         do iatom = 1, natoms
          in1 = imass(iatom)
          if (iorder .le. gear_order) then
            write (64,102) nzx(in1), (xdot(iorder,i,iatom),i=1,3)
          else
            write (64,102) nzx(in1), 0.0d0, 0.0d0, 0.0d0
          end if
         end do
        end do
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (2x, i4)
101     format (2x, i1, 2x, i6, 2x, f9.3)
102     format (2x, i2, 3(2x,f16.9))

 
        close (unit = 64)
        return
      end subroutine writeout_ac
