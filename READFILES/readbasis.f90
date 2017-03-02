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

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! readbasis.f90
! Program Description
! ===========================================================================
!       This routine reads in the basis file. The format of the basis file
! has been changed compared to the previous Fireball program.  Now the
! basis file only contains the xyz coordinates along with Z of the atom.
! For example, the Si dimer file would now look like:
!
!     2 (number of atoms)
!  14   0.000  0.000  0.000
!  14   2.350  0.000  0.000
!
! ===========================================================================
! Original code written by Otto F. Sankey
!
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
        subroutine readbasis ( nzx, imass)
        use configuration
        use dimensions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
 
        integer, intent (in), dimension (nspecies) :: nzx
 
 
! Output
        integer, dimension (natoms), intent(out) :: imass
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer in1
        integer ispec
        integer nucz
        integer natoms_again
 
        logical zindata

! Procedure
! ===========================================================================
! Initialize all the positions to zero.
        ratom = 0.0d0

! Open the basis file.
        open (unit = 69, file = basisfile, status = 'old')
        read (69, *) natoms_again
        if (natoms .ne. natoms_again) then
          write(*,*) ' Strange error in readbasis'
          stop
        end if
 
! Loop over the number of atoms
        do iatom = 1, natoms
         read (69,*) nucz, ratom(:,iatom)
         zindata = .false.
         do ispec = 1, nspecies
          if (nucz .eq. nzx(ispec)) then
           zindata = .true.
           imass(iatom) = ispec
          end if
         end do
         if (.not. zindata) then
          write (*,*) ' The atomic number, nucz = ', nucz
          write (*,*) ' that is contained in your basis file is '
          write (*,*) ' not contained in your data files - info.dat '
          write (*,*) ' Remake your create data files or fix your basis file.'
          stop
         end if
        end do


        do iatom = 1, natoms
          ratom(:,iatom)=ratom(:,iatom)*rescal
        end do 

! Now write out the basis file information.
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' Atom Coordinates from Basis File: '
        write (*,200)
        write (*,201)
        write (*,200)
        do iatom = 1, natoms
         in1 = imass(iatom)
         symbol(iatom) = symbolA(in1)
         write (*,202) iatom, symbol(iatom), ratom(:,iatom), imass(iatom)
        end do
        write (*,200)
        write (*,*) '  '
 
! Format Statements
! ===========================================================================
200     format (2x, 70('='))
201     format (2x, ' Atom # ', 2x, ' Type ', 5x,   &
     &              ' x ', 8x, ' y ', 8x, ' z ', 6x, ' Species # ')
202     format (3x, i5, 7x, a2, 3(2x,f9.3), 7x, i2)
 
        close (unit = 69)
        return
        end
