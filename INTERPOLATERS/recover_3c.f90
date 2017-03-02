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

! recover_3c.f90
! Program Description
! ===========================================================================
!       This subroutine calculates the BOX Hbox (num_orb(in1) x num_orb(in2))
! of matrix-elements from the list of matrix elements stored in Hlist.
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
        subroutine recover_3c (in1, in2, hlist, hbox)
        use dimensions
        use interactions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: in1, in2
 
        real, intent(in) :: hlist (ME3c_max)
 
! Output
        real, intent(out) :: hbox (numorb_max, numorb_max)
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu, inu
        integer index
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Initialize hbox
        do inu = 1, num_orb(in2)
         do imu = 1, num_orb(in1)
          hbox(imu,inu) = 0.0d0
         end do
        end do
 
! Now, construct Hbox
        do index = 1, index_max3c(in1,in2)
         imu = mu(index,in1,in2)
         inu = nu(index,in1,in2)
         hbox(imu,inu) = hlist(index)
        end do
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
