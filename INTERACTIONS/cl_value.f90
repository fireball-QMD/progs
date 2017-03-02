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

! cl_value.f90
! Program Description
! ===========================================================================
!       This routine returns the Kleinman Bylander cl values for atom itype.
! The "raw" date is read in vnl.z1.z2.dat by  program read2c. We include up to
! 5 non-local (L values) of the pseudopotential. Usually, you will have 2
! (L = 0, 1) and sometimes 3 (L = 2).
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
        subroutine cl_value (itype, cl)
        use dimensions
        use interactions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itype
 
! Output
        real, intent (out), dimension (numorb_max) :: cl
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer index
        integer issh
        integer Lvalue
        integer Lmax
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Initialize to zero
        cl(1:numorb_max) = 0.0d0
 
! We now loop though all shells, and create cl for each orbital.  For example,
! sp^3 has two shells; cl(1) = cl_PP(0) and cl(2) = cl(3) = cl(4) = cl_PP(1).
        index = 0
        do issh = 1, nsshPP(itype)
         Lvalue = lsshPP(issh,itype)
         Lmax = (2*Lvalue + 1)
         do imu = 1, Lmax
          index = index + 1
          cl(index) = cl_PP(issh,itype)
         end do
        end do
 
! Sanity check.
        if (index .ne. num_orbPP(itype)) then
         write (*,*) ' itype = ', itype
         write (*,*) ' index of orbitals for pseudopotential = ',index
         write (*,*) ' Program has num_orbPP = ', num_orbPP(itype)
         write (*,*) ' cl_value: index and num_orbPP DO NOT agree. '
         stop
        end if
 
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
