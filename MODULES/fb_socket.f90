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

! fb_sockets.f90
! Program Description
! ===========================================================================
!       This module add socket variables and soubroutines
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ======================================================================

module fb_sockets
         
         
         ! SOCKET VARIABLES
         integer, parameter :: MSGLEN=12   ! length of the headers of the driver/wrapper communication protocol
         integer socket, inet, port        ! socket ID & address of the server
         character(LEN=1024) :: host

         ! SOCKET COMMUNICATION BUFFERS
         character(LEN=12) :: header
         logical :: isinit=.false., hasdata=.false.
         integer cbuf, rid, nat, i
         character(LEN=2048) :: initbuffer      ! it's unlikely a string this large will ever be passed...
         double precision, allocatable :: msgbuffer(:)
         double precision cell_h(3,3), cell_ih(3,3), virial(3,3), mtxbuf(9), dip(3), charges(3), dummy(3,3,3), vecdiff(3)


!*****************************************************************************
contains
         subroutine coordsFromSocket

	   use configuration
           use constants_fireball
           use f90sockets, only : readbuffer

           call readbuffer(socket, mtxbuf, 9)  ! Cell matrix
           cell_h = RESHAPE(mtxbuf * abohr, (/3,3/))
           call readbuffer(socket, mtxbuf, 9) ! Inverse of the cell matrix (so we don't have to invert it every time here)
           cell_ih = RESHAPE(mtxbuf * abohr, (/3,3/))

           ! The wrapper uses atomic units for everything, and row major storage.
           ! At this stage one should take care that everything is converted in the
           ! units and storage mode used in the driver.
           cell_h = transpose(cell_h)
           cell_ih = transpose(cell_ih)
           ! We assume an upper triangular cell-vector matrix
           !volume = cell_h(1,1)*cell_h(2,2)*cell_h(3,3)
           call readbuffer(socket, cbuf)       ! The number of atoms in the cell
           nat = cbuf

           allocate (msgbuffer(3*nat))          

           call readbuffer(socket, msgbuffer, nat*3)

           do i = 1, nat
               ratom(:,i) = msgbuffer(3*(i-1)+1:3*i)*abohr
           enddo

         end subroutine

         subroutine forcesToSocket

           use energy, only : etot
           use forces, only : ftot
           use constants_fireball
           use configuration, only : vatom
           use f90sockets, only : writebuffer
           integer :: bas1, bas2
           integer :: numfrags
           integer :: ifrag
           integer :: fatom
           integer :: fx, fy, fz
           logical :: readsup


           do i = 1, nat
             msgbuffer(3*(i-1)+1:3*i) = ftot(:,i) / eq2
           enddo

           pot = etot / Hartree

           virial = 0.0d0 ! No viral tensor in fireball, transpose(virial)

           call writebuffer(socket,"FORCEREADY  ",MSGLEN)
           call writebuffer(socket,pot) ! Writing the potential 
           call writebuffer(socket,nat)  ! Writing the number of atoms
           call writebuffer(socket,msgbuffer,3*nat) ! Writing the forces
           call writebuffer(socket,reshape(virial,(/9/)),9)  ! Writing the virial tensor, NOT divided by the volume
           cbuf = 7 ! Size of the "extras" string
           call writebuffer(socket,cbuf)  ! This would write out the "extras" string, but in this case we only use a dummy string.
           call writebuffer(socket,"nothing",7)

           deallocate (msgbuffer)

         end subroutine

end module fb_sockets
