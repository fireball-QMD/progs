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

! readephc.f90
! Program Description
! ===========================================================================
!       This reads in the epch.optional file
!
! ===========================================================================
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine readephc (natoms)
        use dynamo
        use interactions
        use kpoints

        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in)        :: natoms

! Ouput
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer i
        character (len = 30) filein

! Procedure
! ===========================================================================
        write (*,*) ' ----------------------------------------'
        write (*,*) ' You want to calculate e-ph coupling'
        write (*,*) ' Read the parameters from ephc.optional'

! open input file
        open (file ='ephc.optional', unit=17, status='old')

! read name of output file 
        read (17,*) fileephc
        write (*,200) fileephc

! read name of file containing reference eigenvalues 
        read (17,*) filein
        write (*,300) filein

! read effective temperature 
        read (17,*) temp_ephc
        write (*,310) temp_ephc

! read effective temperature 
        read (17,*) mass_ephc
        write (*,320) mass_ephc

! read number of e-ph couplings
        read (17,*) nephc
        write (*,100) nephc

! 1. remember we have two displacement ('+' and '-')for each elementary 
! 2. 0-term contains reference eigenvalue
        allocate (deigen(nephc,0:3*natoms))
        allocate (eiglist(nephc))
        allocate (eigref(nephc))

! read ID eigenvalues for which will be e-ph evaluated
        read (17,*)  (eiglist(i), i =1,nephc)
        write (*,350)
        write (*,*) eiglist(:)

! read number of vibrational modes for e-ph couplings
        read (17,*) nvmodes
        write (*,110) nvmodes
        if (nvmodes .eq. 0) then
          nvmodes = 3*natoms
          write (*,*) ' All vib. modes will be analyzed'
          allocate (idvmode(nvmodes))
          do i = 1, nvmodes
            idvmode(i) = i
          enddo
        else  
          allocate (idvmode(nvmodes))
! read vibrational modes for which will be e-ph evaluated
          read (17,*)  (idvmode(i), i=1,nvmodes)
          write (*,360)
          write (*,*) idvmode(:)
        endif 
! close input file
        close (unit = 17)

! openfile containing reference eigenvalues
        open (file =filein, unit=18, status='old')
! read reference iegenvalues
	do i = 1, nephc
         read (18,*) eigref(i)
	enddo 
        close (18)

        write (*,400)

! Format Statements
! ===========================================================================
100    format ('   Number of e-ph couplings to be evaluate:      ', i4)
110    format ('   Number of vibrational modes to be considered: ', i4)
200    format ('   e-ph couplings will be written into file:     ',a30)
300    format ('   File with reference eigenvalues:              ',a30)
310    format ('   Effective temperature:                        ',f12.8)
320    format ('   Effective mass:                               ',f12.8)
350    format ('   Eigenvalues to be considered:                 ')
360    format ('   Vibrational modes to be considered:           ')
400    format ('  ===========================================================')

  return
end subroutine readephc
