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

! readquench.f90
! Program Description
! ===========================================================================
!       This routine reads in the information from the quench.optional file.
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
        subroutine readquench (iquench, dt, energytol, forcetol,             &
     &                         iensemble, T_initial, T_want, taurelax)

        use configuration
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iquench
        integer, intent (in) :: iensemble
 
        real, intent (in) :: dt
        real, intent (in) :: T_initial
 
! Output
        real, intent (out) :: energytol
        real, intent (out) :: forcetol
        real, intent (out) :: taurelax
        real, intent (out) :: T_want
 
! Local Parameters and Data Declaration
! ===========================================================================
        Namelist /quench/ energytol, forcetol, T_want, taurelax
 
! Local Variable Declaration and Description
! ===========================================================================
         logical isfile
         character (len = 30) str
         logical isstr
         
! Procedure
! ===========================================================================
        if (iquench .eq. 0 .or. iquench .eq. -1 .or. iquench .eq. -2 .or. iquench .eq.-3 .or. iquench .eq.-6  ) then

! set default options
         energytol = 0.0001  ! eV
         forcetol  = 0.05    ! eV/A
         T_want    = 0.0     ! K
         taurelax  = 50.0    
         
         write (*,*) '  '
         write (*,100)
         write (*,*) ' You have chosen iquench = ', iquench
         write (*,*) ' Optionally read information from the fireball.in file. '
         write (*,*) '  '
         
! checkout existence of fireball.in file
         inquire (file = initfile, exist = isfile)
! file fireball.in exists so let's read it 
         if (isfile) then 

           write (*,*) '  '
           write (*,100)
           write (*,*) ' Now reading file fireball.in '
           write (*,*) '  '
! open param file
           open (unit = 60, file = initfile, status = 'old')

! section QUENCH
! check if the section QUENCH is listed
          str = '&QUENCH'
          call issection (str,isstr)
          if (isstr) read (60, NML=quench)
          
         endif

! write out resume 
         write (*,*) ' In the event that you are searching for an energy '
         write (*,*) ' minimum, you can also choose energy and force '
         write (*,*) ' tolerances - energy_tol and force_tol. These tolerances '
         write (*,*) ' will stop execution after they are achieved. Input a '
         write (*,*) ' number for all of the following, even if option is not '
         write (*,*) ' being used! '
 
         write (*,*) ' Input energy_tol, force_tol, T_want, taurelax: '
         write (*,*) '  '
         write (*,*) ' energy tolerance = ', energytol,' [eV]'
         write (*,*) ' force tolerance  = ', forcetol,' [eV/A]'
         write (*,*) ' T_want = ', T_want,' [K]'
         if (iensemble .gt. 0 .and. T_want .eq. 0.0d0) then
          write (*,*) ' ******************** NOTE ********************* '
          write (*,*) ' You have chosen to do constant temperature MD, '
          write (*,*) ' but T_want = ', T_want
          write (*,*) ' Setting T_want = T_initial! '
          write (*,*) ' ******************** NOTE ********************* '
          T_want = T_initial
         end if
         write (*,*) ' taurelax = ', taurelax
 
         if (iquench .eq. -2) then
          write (*,*) '  '
          write (*,*) ' You have chosen the annealing option. '
          write (*,*) ' The quench rate = ', (dt/taurelax)
          write (*,*) ' The annealing temperature = ', T_want
         end if
         write (*,100)
         
! write out information into param.dat
         open (unit = 50, file = 'param.dat', position='append' ,status='old')
         write (50, *) ' QUENCH:'
       
         write (50, *) '  energytol         : ',energytol
         write (50, *) '  forcetol          : ',forcetol
         write (50, *) '  T_want            : ',T_want
         write (50, *) '  taurelax          : ',taurelax
         close (50)
        end if
 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
 
        return
        end
