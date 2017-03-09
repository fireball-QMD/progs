!                             @@Copyright 2005
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek
! Brigham Young University - Hao Wang

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang 
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

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


! fireball.f90
! Program Description
! ===========================================================================
!       This is the main driver program for the FIREBALL package.
!
! ===========================================================================
! Code written by:
! Jesús I. Mendieta & Daniel G. Trabada
! Departamento de físca teórica de la materia condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================


        module fmpi
         include 'mpif.h'
         integer rank, ierr_mpi, intercomm, status, size,mpi
         character*(256) serv_name
         character*(MPI_MAX_PORT_NAME) port_name
        end module fmpi

        subroutine connect()
         use fmpi
         implicit none
         serv_name = 'fireball_server'
         call MPI_INIT ( ierr_mpi )
         call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr_mpi)
        !call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr_mpi)
         if(rank==0) then
          call MPI_OPEN_PORT(MPI_INFO_NULL, port_name, ierr_mpi)
          print *,"port_name =",port_name
          call MPI_PUBLISH_NAME(serv_name, MPI_INFO_NULL, port_name, ierr_mpi)
          call MPI_COMM_ACCEPT(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF,intercomm, ierr_mpi)
         end if
       end subroutine



        program fireball_server
        use mpi_main
        use interactions
        use options
        use configuration
        use charges
        use fmpi

  
! Parameters and Data Declaration
! ===========================================================================
 
! Variable Declaration and Description
! ===========================================================================
 
! --------------------------------------------------------------------------
! Timer (Intel Fortran)
! --------------------------------------------------------------------------
        real time_begin
        real time_end
        !real,allocatable :: R(:,:)
        real,allocatable,dimension(:,:) :: R

! Procedure
! ===========================================================================

! ***************************************************************************
!                           MPI Specific Stuff
! ***************************************************************************
        print *,'PROC',nprocs, my_proc
        call connect()
        if(rank==0)then
         print *,'first time : Client should be created fireball.in,'
         print *,' Fdata.optional and input.bas'
         print *,'needed to initbasis() and readdata()'
         print *,'natoms not change from input.bas'

         call MPI_RECV(mpi,1,MPI_INTEGER,  0, 0,intercomm, status, ierr_mpi) 
         call cpu_time (time_begin)
         call initbasics () 
         call readdata ()
         allocate(R(3,natoms))
        !basisfile read from fireball.in

         if (iqmmm .eq. 1) then
           !!   call main_loop_qmm_server ()
           !Send Client ftot(3xnatom), vatom(3xnatom), etot, QLowdin_TOT(natom), ftot_rij_qmmm_atoms(nqmmm_atoms)
           !! SEND....
         else
           !**************** RECV **************************   
           ! natoms == natoms read in input.bas, Z == Z read input.bas
           ! call MPI_RECV(natoms,1, MPI_INTEGER, 0, 0, intercomm, status,ierr_mpi)
           ! call MPI_RECV(S,2*natoms, MPI_CHARACTER, 0, 0, intercomm,status,ierr_mpi)
           call MPI_RECV(R,3*natoms, MPI_REAL, 0, 0, intercomm,status,ierr_mpi)
 
           print *,'recv'
           do iatom = 1, natoms
            ratom(:,iatom)=R(:,iatom)+shifter
            write (*,'(2x, 3(2x,f12.6))') R(:,iatom)
           end do
           xdot(0,:,1:natoms) = ratom(:,1:natoms)        
  
           ! *************** RUN Fireball ********************
 
           call main_loop ()
 
           !**************** SEND **************************   
           print *,'send'
           do iatom = 1, natoms
            R(:,iatom)=ratom(:,iatom)-shifter 
            write (*,'(2x, 3(2x,f12.6))') R(:,iatom)
           end do
         
           call cpu_time (time_end)
           write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'
           call MPI_SEND(R,3*natoms, MPI_REAL, 0, 0, intercomm, ierr_mpi)
 
         end if
         call MPI_UNPUBLISH_NAME(serv_name, MPI_INFO_NULL, port_name, ierr_mpi)
         call MPI_CLOSE_PORT(port_name, ierr_mpi)
        end if
        call MPI_FINALIZE(ierr_mpi)
 
! timer
!       write (*,1001) omp_get_max_threads()

! Format Statements
! ===========================================================================
100     format (2x, 70('='))
1001    format (2x, '    Number OpenMP threads: ', i3)
 
        stop
      end program fireball_server
