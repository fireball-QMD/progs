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
! fireball-qmd is a free (GPLv3) open project.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


! main_loop_socket.f90
! Program Description
! ===========================================================================
!       This routine performs socket loop
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine main_loop_socket ()

        use options
        use configuration
        use options
        use MD
        use forces
        use constants_fireball
        use energy
        use optimization
        use scf
        use fb_sockets
        use f90sockets, only : open_socket, readbuffer, writebuffer

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================
      integer :: bas1, bas2
      integer :: numfrags
      integer :: ifrag
      integer :: fatom
      integer :: fx, fy, fz
      logical :: readsup
! Local Variable Declaration and Description
! ===========================================================================
       integer itime_step

! Procedure
! ===========================================================================

! in case of constrain DFT, first we perform SCF loop for the ground state;
! This will be used later as a reference state
       if(icDFT .eq. 1 ) then
        call initcDFT ()
        call getenergy (0)
        cDFT_active = .true.
       endif  ! end icDFT 

! ===========================================================================


       call open_socket(socket, inet, port, host)

       itime_step = 1

       do while (.true.)

         call readbuffer(socket, header, MSGLEN)

         if (trim(header) == "STATUS") then

            ! The wrapper is inquiring on what we are doing
            if (.not. isinit) then
               call writebuffer(socket,"NEEDINIT    ",MSGLEN)  ! Signals that we need initialization data "
            elseif (hasdata) then
               call writebuffer(socket,"HAVEDATA    ",MSGLEN)  ! Signals that we are done computing and can return forces"
            else
               call writebuffer(socket,"READY       ",MSGLEN)  ! We are idling and eager to compute something
            endif

         elseif (trim(header) == "INIT") then     ! The driver is kindly providing a string for initialization

            call readbuffer(socket, rid)
            call readbuffer(socket, cbuf)
            call readbuffer(socket, initbuffer, cbuf)
            isinit=.true. ! We actually do nothing with this string, thanks anyway. Could be used to pass some information (e.g. the input parameters, or the index of the replica, from the driver

         elseif (trim(header) == "POSDATA") then ! The driver is sending the positions of the atoms. Here is where we do the calculation!

            call coordsFromSocket

            call scf_loop (itime_step)

            call postscf ()

            call getenergy (itime_step)
        
            call getforces_socket (itime_step)

            !.........FRAGMENTS

         !    inquire (file = 'FRAGMENTS', exist = readsup)
         !    if(.not. readsup) then
         !    else
         !       open(unit = 33, file = 'FRAGMENTS', status = 'old')
         !       read(33,*) bas1
         !       read(33,*) bas2
         !       read(33,*) numfrags
         !       do ifrag = 1,numfrags
         !          read(33,*) fatom, fx, fy, fz
         !          write(*,*) 'Ankais : ',fatom,fx,fy,fz
         !          if (fx .eq. 1) then
         !             ftot(fx,fatom) = 0.0d0
         !             vatom(fx,fatom) = 0.0d0
         !          end if !end if fx .eq. 1
         !          if (fy .eq. 1) then
         !             ftot(fy,fatom) = 0.0d0
         !             vatom(fy,fatom) = 0.0d0
         !          end if !end if fy .eq. 1
         !          if (fz .eq. 1) then
         !             ftot(fz,fatom) = 0.0d0
         !             vatom(fz,fatom) = 0.0d0
         !          end if !end if fz .eq. 1
         !       end do !end do ifrag
         !      close(33)
         !    end if !end if FRAGMENTS exists


            !.........END FRAGMENTS

            hasdata = .true. ! Signal that we have data ready to be passed back to the wrapper

         elseif (trim(header) == "GETFORCE") then  ! The driver calculation is finished, it's time to send the results back to the wrapper


            call forcesToSocket

            itime_step = itime_step + 1
            hasdata = .false.

         elseif (trim(header) == "EXIT") then
           exit
         else
            write(*,*) " Unexpected header ", trim(header)
            STOP "ENDED"
         endif


       end do


! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))

        return
        end subroutine main_loop_socket

