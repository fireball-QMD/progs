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

! readlvs.f90
! Program Description
! ===========================================================================
!       This program read in the lattice vectors from the *.lvs file.
!
! ===========================================================================
! Code rewritten by:
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
        subroutine readlvs (lvsfile, a1vec, a2vec, a3vec, icluster, rescal)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        character (len = 40), intent (in) :: lvsfile
        integer, intent(in) :: icluster
        real, intent(in) :: rescal
 
! Output
        real, intent (out), dimension (3) :: a1vec, a2vec, a3vec
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        real a1, a2, a3
        integer ilvs
 
! Procedure
! ===========================================================================
! Open the lattice vector file.
        if (icluster .ne. 0) then
          write (*,*) '  '
          write (*,*) ' You are doing a cluster calculation '
          write (*,*) ' Lattice vectors are not needed or read in '
          write (*,*) '  '
          a1vec(1) = 100
          a1vec(2) = 0
          a1vec(3) = 0
          a2vec(1) = 0 
          a2vec(2) = 100
          a2vec(3) = 0
          a3vec(1) = 0 
          a3vec(2) = 0
          a3vec(3) = 100
          return
        end if

        open (unit = 72, file = lvsfile, status = 'old')
 
        write (*,*) '  '
        write (*,*) ' Lattice Vectors from File: '
        write (*,100)
        write (*,*) ' (a1x,a1y,a1z) (in angstroms): '
        read (72,*) a1vec(:)
        do ilvs=1,3
          a1vec(ilvs)=a1vec(ilvs)*rescal
        end do
        write (*,101) a1vec(:)
        a1 = sqrt(a1vec(1)**2 + a1vec(2)**2 + a1vec(3)**2)
        write (*,*) ' Magnitude of a1 = ', a1
 
        write (*,*) '  '
        write (*,*) ' (a2x,a2y,a2z) (in angstroms): '
        read (72,*) a2vec(:)
        do ilvs=1,3
          a2vec(ilvs)=a2vec(ilvs)*rescal
        end do
        write (*,101) a2vec(:)
        a2 = sqrt(a2vec(1)**2 + a2vec(2)**2 + a2vec(3)**2)
        write (*,*) ' Magnitude of a2 = ', a2
 
        write (*,*) '  '
        write (*,*) ' (a3x,a3y,a3z) (in angstroms): '
        read (72,*) a3vec(:)
        do ilvs=1,3
          a3vec(ilvs)=a3vec(ilvs)*rescal
        end do
        write (*,101) a3vec(:)
        a3 = sqrt(a3vec(1)**2 + a3vec(2)**2 + a3vec(3)**2)
        write (*,*) ' Magnitude of a3 = ', a3
 
        close (unit = 72)
 
        write (*,100)
        write (*,*) '  '
 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
101     format (2x, 3f10.4) 
 
        return
        end
