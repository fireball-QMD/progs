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

! readdata_2c.f90
! Program Description
! ===========================================================================
!       This routine reads the data from the 2-center integral files. When
! read, the information is stored in the array xintegral_2c.  This array
! is the field that stores all non-vanishing matrix elements for a general
! 2-center integral.  There are maximal ME2c_max non-vanishing matrix
! elements given on a grid of maximal nfofx data points.  The exact dimensions
! for a given interaction, and a given pair of atoms are numz and num_nonzero.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine readdata_2c (interaction, iounit, num_nonzero, numz, zmax,&
     &                          itype, in1, in2, ioff2c)
        use dimensions
        use integrals
        use interactions
        use constants_fireball
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: in1, in2
        integer, intent (in) :: interaction
        integer, intent (in) :: ioff2c
        integer, intent (in) :: iounit
        integer, intent (in) :: itype
        integer, intent (in) :: num_nonzero
        integer, intent (in) :: numz
        real, intent (in) :: zmax
 
! Output is xintegral_2c
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ipoint
        integer integral
 
        real xoff
 
        real, dimension (ME2c_max, nfofx) :: gstore
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! This set to zero unset elements (see 'diagnostics.input')
        xoff = real(ioff2c)
        if (interaction .ne. 8) then
         do ipoint = 1, numz
          read (iounit,*) (gstore(integral,ipoint), integral = 1, num_nonzero)
         end do
         do ipoint = 1, numz
          do integral = 1, num_nonzero
           xintegral_2c(integral,ipoint,itype,in1,in2) =                     &
     &      gstore(integral,ipoint)*xoff
          end do
         end do
         if (superspline) then
          do integral = 1, num_nonzero
           call buildspline_1d (integral, numz, itype, in1, in2, zmax,       &
     &                          interaction)
          end do
         end if
        else
         do ipoint = 1, numz
          read (iounit,*) gstore(1,ipoint)
          xintegral_2c(1,ipoint,itype,in1,in2) = gstore(1,ipoint)*xoff
         end do
         if (superspline) then
          call buildspline_1d (1, numz, itype, in1, in2, zmax, interaction)
         end if
        end if

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end subroutine readdata_2c
