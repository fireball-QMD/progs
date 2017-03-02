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
 
! diagnostics.f90
! Program Description
! ===========================================================================
!       This routine reads in the diagnostics.input file. This file allows
! the user to turn interactions on or off.  Usually for debugging purposes.
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
        subroutine diagnostics (ioff2c, ioff3c, itestrange, testrange)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Output
        integer, intent (out) :: itestrange
        integer, intent (out), dimension (1:24) :: ioff2c

 
        integer, intent (out), dimension (1:4) :: ioff3c
 
        real, intent (out) :: testrange
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer interaction
        logical readdiag

! Procedure
! ===========================================================================
! Note ioff=1 means the file is ON. ioff=0 means the file is OFF.
! We multiply the data file by 0 if ioff is off.
        inquire (file = 'diagnostics.input', exist = readdiag)
        if (.not. readdiag) then
          write (*,*) ' No diagnostics.input file to read. Will run normally '
          ioff2c(:) = 1
          ioff3c(:) = 1
          itestrange= 1
          testrange = 0
          return
        end if
        write (*,*) ' Reading diagnostics.input '
        open (unit = 26, file = 'diagnostics.input', status = 'old')
 
! ****************************************************************************
! 2 Center !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           interaction    subtype
!               1            0            overlap
!               2            0            vna   ontopl
!               2            isorp        vna   ontopl shell isorp
!               3            0            vna   ontopr
!               3            isorp        vna   ontopr shell isorp
!               4            0            vna   atom
!               4            isorp        vna_  atom  shell isorp
!               5            0            non-local
!               6            0,4          xc ontop
!               7            0,4          xc atom-atom
!               8            0,4          xc correction
!               9            0            z-dipole
!               10           0            y-dipole
!               11           0            x-dipole
!               12           0            coulomb
!               13           0            kinetic
!               14           0            extended hubbard
!               15           isorp        den_ontopl
!               16           isorp        den_ontopr
!               17           isorp        den_atom
!               18           isorp        denS_ontopl
!               19           isorp        denS_ontopr
!               20           isorp        denS_atom
!               21           isorp        overlapS
! 3 Center !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!               1            0            neutral atom
!               2            0            exchange-correlation Horsfield
!               3            0            average density OLSXC 
!               4            0            average density OLSXC (spher approx)
! ****************************************************************************
        read (26,*)
        read (26,*)
        read (26,*)

        do interaction = 1, 23

! Read ioff for the two-center interactions.
         read (26,*) ioff2c(interaction)
         if (ioff2c(interaction) .eq. 0) then
          write (*,*) ' ********************************************** '
          write (*,*) ' *******************WARNING******************** '
          write (*,*) ' *******************WARNING******************** '
          write (*,*) ' *******************WARNING******************** '
          write (*,*) '  '
          write (*,*) ' You are setting some two-center interactions   '
          write (*,*) ' to zero! Make sure you know what you are doing!'
          write (*,*) '  '
          write (*,*) ' *******************WARNING******************** '
          write (*,*) ' *******************WARNING******************** '
          write (*,*) ' *******************WARNING******************** '
          write (*,*) ' ********************************************** '
         end if
        end do
 
        do interaction = 1, 4
         read (26,*) ioff3c(interaction)
         if (ioff3c(interaction) .eq. 0) then
          write (*,*) ' ********************************************** '
          write (*,*) ' *******************WARNING******************** '
          write (*,*) ' *******************WARNING******************** '
          write (*,*) ' *******************WARNING******************** '
          write (*,*) '  '
          write (*,*) ' You are setting some three-center interactions '
          write (*,*) ' to zero! Make sure you know what you are doing!'
          write (*,*) '  '
          write (*,*) ' *******************WARNING******************** '
          write (*,*) ' *******************WARNING******************** '
          write (*,*) ' *******************WARNING******************** '
          write (*,*) ' ********************************************** '
         end if
        end do
 
        read (26,*) itestrange, testrange
        if (itestrange .eq. 0) then
         write (*,*) ' ********************************************** '
         write (*,*) ' *******************WARNING******************** '
         write (*,*) ' *******************WARNING******************** '
         write (*,*) ' *******************WARNING******************** '
         write (*,*) '  '
         write (*,*) ' You are testing the range in fireball.f! This  '
         write (*,*) ' is controlled by diagnostics.input! '
         write (*,*) ' The range you chose is ', testrange
         write (*,*) '  '
         write (*,*) ' *******************WARNING******************** '
         write (*,*) ' *******************WARNING******************** '
        end if
 
! Format Statements
! ===========================================================================
 
        return
        end
