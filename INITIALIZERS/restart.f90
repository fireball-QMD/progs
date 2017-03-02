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

! restart.f90
! Program Description
! ===========================================================================
!       This routine reads in the necessary files - *.ac and *.xv in order
! to restart the MD simulation in progress.  Call this routine if the 
! first step in the script.input file is not equal to 1!
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
        subroutine restart (nstepi, dt, acfile, xvfile,             &
     &                      T_average, T_previous, time)
        use configuration
        use constants_fireball
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (inout) :: nstepi

        real, intent (in) :: dt
        real, intent (inout) :: time

        character (len = 30), intent (in) :: acfile
        character (len = 30), intent (in) :: xvfile

! Output
        real, intent (out) :: T_average
        real, intent (out) :: T_previous

 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description 
! ===========================================================================
        integer iatom
        integer iorder
        integer iorder_in
        integer itime_step_in
        integer itime_step
        integer nzx
        integer num_atoms
        integer nstepi_ac

        real delta
        real tkinetic
        real time_save
        real T_instantaneous
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        write (*,*) '  '
        write (*,100)
        write (*,*) ' Reading restart information! '
        write (*,100)

! Read in the acfile.
        write (*,*) ' Reading accelerations from the acfile. '
        open (unit = 64, file = acfile, status = 'old')
        read (64,*) num_atoms
        if (num_atoms .ne. natoms) then
         write (*,*) ' The acfile that you are reading from does not '
         write (*,*) ' correspond to the chosen basis file. ' 
         write (*,*) ' Synchronize everything and start again! '
         stop
        end if
        do iorder = 2, min(gear_order,5) ! ac file only has up to fifth-order
         read (64,101) iorder_in, nstepi_ac, time
!         nstepi_ac = nstepi_ac + 1
         if (iorder_in .ne. iorder) then
          write (*,*) '  '
          write (*,*) ' Something is wrong here! '
          write (*,*) ' The iorder numbers are not monotonically increasing. '
          stop
         end if
         do iatom = 1, natoms
          read (64,102) nzx, xdot(iorder,:,iatom)
         end do
        end do
        close (unit = 64)
 
! Read in the xvfile.
! vel-note: only last xv-record is saved for restart (see UTIL/writeout_xv.f90)
        write (*,*) ' Reading positions and velocities from the xvfile. '
        open (unit = 63, file = xvfile, status = 'old')
        read (63,*) num_atoms
        if (num_atoms .ne. natoms) then
         write (*,*) ' The acfile that you are reading from does not '
         write (*,*) ' correspond to the chosen basis file. ' 
         write (*,*) ' Synchronize everything and start again! '
         stop
        end if
!        time = - dt
!        do itime_step = 1, nstepi - 1
!         time_save = time 
         read (63,201) itime_step_in, time
         write (*,*) ' itime_step_in, nstepi = ', itime_step_in, nstepi
         if ( itime_step_in .ne. nstepi ) then
          write (*,*) '  '
          write (*,*) ' The time step that you have chosen is not consistent '
          write (*,*) ' with the time step of the previous run. '
          write (*,*) ' itime_step_in, nstepi = ', itime_step_in, nstepi
          write (*,*) ' Synchronize time steps and start again! '
          stop
         endif 
!         if (itime_step_in .ne. itime_step) then
!          write (*,*) '  '
!          write (*,*) ' Something is wrong here! '
!          write (*,*) ' The step numbers are not monotonically increasing. '
!          stop
!         end if
         do iatom = 1, natoms
          read (63,202) nzx, ratom(:,iatom), vatom(:,iatom)
          write(*,202) nzx, ratom(:,iatom), vatom(:,iatom)
         end do
!        end do
        close (unit = 63)

         nstepi =nstepi + 1

!        delta = time - time_save
!        if (abs(delta - dt) .gt. 1.0d-3) then
!         write (*,*) '  '
!         write (*,*) ' The time step that you have chosen is not consistent '
!         write (*,*) ' with the time step of the previous run. '
!         write (*,*) ' time, time_save = ', time, time_save
!         write (*,*) ' delta, dt = ', delta, dt
!         write (*,*) ' Synchronize time steps and start again! '
!         stop
!        end if
        xdot(0,:,1:natoms) = ratom(:,1:natoms)
        xdot(1,:,1:natoms) = vatom(:,1:natoms)

! Update time
!        time = dt*(nstepi - 1) 

! Calculate the average temperature to date.
!        T_average = 0.0d0
!        do itime_step = 1, nstepi - 1
         tkinetic = 0.0d0
         do iatom = 1, natoms
          tkinetic = tkinetic                                                &
     &     + (0.5d0/fovermp)*xmass(iatom)                                    &
     &      *(vatom(1,iatom)**2 + vatom(2,iatom)**2 + vatom(3,iatom)**2)
         end do
         T_instantaneous = (2.0d0/3.0d0)*tkinetic*kconvert/natoms
         T_previous = T_instantaneous
!         T_average = T_average + T_instantaneous
         T_average = T_instantaneous
!        end do
!        T_average = T_average/(nstepi - 1)
        write (*,*) ' T_previous = ', T_previous, ' T_average = ', T_average
        write (*,100)
       
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
101     format (2x, i1, 2x, i6, 2x, f9.3)
102     format (2x, i2, 3(2x,f16.9))
201     format (2x, i6, 2x, f9.3)
202     format (2x, i2, 3(2x,f16.9), 3(2x,f16.9)) 

        return
      end subroutine restart
