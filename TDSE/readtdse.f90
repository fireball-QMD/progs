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

! readtdse.f90
! Program Description
! ===========================================================================
!       This routine reads input parameters for TDSE simulation
!
! ===========================================================================
! Code written by:
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
        subroutine readtdse ( )

        use MD
        use configuration
        use tdse
        use charges
        use options

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

! Local Parameters and Data Declaration
! ===========================================================================
        Namelist /tds/ netime, nexcite, idelec, eband, hoccup, np2es

! Local Variable Declaration and Description
! ===========================================================================
         logical isfile
         character (len = 30) str
         logical isstr
         integer i

! Procedure
! ===========================================================================

! set default options
        netime = 1000  ! no. of electron time steps within one MD step
        nexcite = 0  ! no. of excitations
        idelec = 0
        hoccup = 0.0d0
        eband = 0
        np2es = -1
        isp2es = .false.

        write (*,*) '  '
        write (*,100)
        write (*,*) ' You have chosen itdse = ', itdse
        write (*,*) ' Optionally read information from the fireball.in file. '
        write (*,*) '  '

! checkout existence of fireball.in file
        inquire (file = 'fireball.in', exist = isfile)
! file fireball.in exists so let's read it
        if (isfile) then

          write (*,*) '  '
          write (*,100)
          write (*,*) ' Now reading file fireball.in '
          write (*,*) '  '
! open param file
          open (unit = 60, file = initfile, status = 'old')

! section TDSE
! check if the section TDSE is listed
         str = '&TDS'

         call issection (str,isstr)
         if (isstr) then
          write (*,*) ' ---- The section &TDS found!  ----'
          read (60, NML=tds)
         endif
        endif

! set electron time step
        dte = dt/netime
! set nelec
        if (abs(qstate) .gt. 0.00001) then
         write (*,*) ' qstate =',qstate
         write (*,*) ' ztot = ',ztot
         write (*,*) ' We cannot deal with charged systems at the moment'
         write (*,*) ' Must stop, sorry! '
         stop
        else
         nelec = ztot
        endif

        write (*,*) ' Number of excitations = ', nexcite
        if (nexcite .gt. 0) then
         do i = 1, nexcite
          write (*,*) ' --- excitation no.', i
          write (*,*) ' elec no. ',idelec(i),' to band no. ',eband(i),' h-occup =',hoccup(i)
         end do ! do nexcite
        endif ! nexcite


! write out information into param.dat
        open (unit = 50, file = 'param.dat', position='append' ,status='old')
        write (50, *) ' TDS:'

        write (50, *) '  netime            : ',netime
        write (50, *) '  dte               : ',dte
        write (50, *) '  nelectrons        : ',nelec
        write (50, *) '  nexcite           : ',nexcite
        do i = 1,nexcite
         write (50, *) '   pair e-h          : ',idelec(i),eband(i),hoccup(i)
        enddo
        if (np2es .lt. 0) then
         write (50, *) '  np2es           :    OFF'
        else
         write (50, *) '  np2es           : ',np2es
        endif
        close (50)


! Format Statements
! ===========================================================================
100     format (2x, 70('='))

        return
        end
