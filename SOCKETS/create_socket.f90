! copyright info:
!
!                             @Copyright 2003
!                            Fireball Committee
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

 
! create_socket.f90
! Program Description
! ===========================================================================
! Creates a socket to talk to another program that is driving fireball.
! Actual work is done in C code.
! ===========================================================================
! Code written by:
! Kurt R. Glaesemann
! LLNL
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine create_socket()
        implicit none
 
! Local Variable Declaration and Description
! ===========================================================================
        integer port ! The port number used to talk to other code
        integer pid  ! My process ID
        integer, external :: getpid
        character (len=10) my_proc_ascii
 
! Procedure
! ===========================================================================
        pid = getpid()
        write(my_proc_ascii, '(i10.10)') pid
        write(*,*) ' Getting the port number '
        open (unit=44, file='pimc.port'//my_proc_ascii, status='old')
        read (44,*) port
        close (unit=44)
        write(*,*) ' Using port number ', port
        call soc_init(port)
        write(*,*) ' Done creating socket '

        return
        end
