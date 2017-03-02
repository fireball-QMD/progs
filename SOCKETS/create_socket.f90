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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
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
