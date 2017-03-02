! copyright info:
!
!                             @Copyright 2002
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
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
! University of Regensburg - Juergen Fritsch

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
! get_umbrella.f90 
! Program Description
! ===========================================================================
!       This is the main subroutine for the umbrella sampling potential option
! At each call the value of the reaction coordinate, umbrella energy and 
! total energy (before umbrella) is written... The reaction coordinate
! value is sufficient for the WHAM analysis...
!
! ===========================================================================
! Code written by:
! Hilaire Chevreau
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! (801) 585-1078 (office)      email: lewis@hec.utah.edu
! (801) 581-4353 (fax)         Web Site: http://
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine get_umbrella (nstepi, nstepf, itime_step, time, natoms,   &
     &                           iwrtfpieces, ratom, etot, ftot)
        use umbrella
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itime_step          
        integer, intent (in) :: iwrtfpieces          
        integer, intent (in) :: nstepi          
        integer, intent (in) :: nstepf         
        integer, intent (in) :: natoms     

        real, intent (in) :: time    ! Current time of the MD simulation

        real, intent (in), dimension (3, natoms) :: ratom

! Output
        real, intent (inout) :: etot                         ! total energy

        real, intent (inout), dimension (3, natoms) :: ftot  ! total force  

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        real coord_reac_inter                ! Reaction coordinate (MD step i)  
        real etot_initial                    ! total energy before umbrella
        real umb_e                           ! umbrella energy
 
! Procedure
! ===========================================================================
! Initialize 

! If doing umbrella sampling, read in umbrella.optional file. This file 
! contains all the information needed for the umbrella sampling MD runs.
        if (itime_step .eq. nstepi) then          ! for restart option...
         write (*,*) ' Umbrella Sampling Option '
         write (*,*) '   '
 
         call readumbrella (natoms, ratom, nstepf)
        endif                                               

! Now, calculation of the umbrella sampling energy due to the umbrella sampling
! constraint between pair of atoms.  The constraint is only on the distance at 
! this moment.
        call assemble_umbrella (umb_e, coord_reac_inter, natoms, ratom)
        write (*,*) ' The umbrella energy is : ', umb_e
        write (*,*) '   ' 
        etot_initial = etot
        etot = etot + umb_e

! Umbrella sampling forces contributions
        write (*,*) ' Determine umbrella sampling forces contributions '
        write (*,*) '   ' 
        call Dassemble_umbrella (natoms, iwrtfpieces, ratom, ftot)

        if (time .ge. umb_time_start) then
         open (unit = 10, file = umb_output_file, status = 'unknown',        &
     &         position = 'append')
         write (10,*) coord_reac_inter, umb_e, etot_initial 
         close (unit = 10)
        end if

! Format Statements
! ===========================================================================
100     format (2x, 70('='))  
 
! ===========================================================================
        return
        end
