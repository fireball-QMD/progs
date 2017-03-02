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

! make_mu2shell.f90
! Program Description
! ===========================================================================
!       This routine determines the shell number of the in1'th atomtype and
! stores that information in the array mu2shell.
!
! ===========================================================================
! Code written by:
! Otto F. Sankey (visiting)
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
        subroutine make_mu2shell (nspecies)
        use dimensions
        use interactions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nspecies
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer in1
        integer issh
        integer imumin, imumax
        integer Lvalue
        integer num_Mvalues
 
! Allocate Arrays
! ===========================================================================
        allocate (mu2shell(numorb_max, nspecies))
 
! Procedure
! ===========================================================================
! The extended hubbard subroutines need array mu2shell(imu,in1) which tells
! the shell number of the imu'th orbital of the in1'th atomtype.
        do in1 = 1, nspecies
         imumin = 0
         imumax = 0
         do issh = 1, nssh(in1)
          Lvalue = lssh(issh,in1)
          num_Mvalues = 2*Lvalue + 1
          imumax = imumin + num_Mvalues
          imumin = imumin + 1
          do imu = imumin, imumax
           mu2shell(imu,in1) = issh
          end do
          imumin = imumax
         end do
        end do
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
