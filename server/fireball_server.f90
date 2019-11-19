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
         call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr_mpi)
         if(rank==0) then
          call MPI_OPEN_PORT(MPI_INFO_NULL, port_name, ierr_mpi)
          print *,"port_name =",port_name
          call MPI_PUBLISH_NAME(serv_name, MPI_INFO_NULL, port_name, ierr_mpi)
          call MPI_COMM_ACCEPT(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF,intercomm, ierr_mpi)
         end if
       end subroutine



        program fireball_server
        use mpi_main
        use fmpi
        use qmmm_module
        use energy
        use forces
        use configuration
        use interactions
  
! Parameters and Data Declaration
! ===========================================================================
 
! Variable Declaration and Description
! ===========================================================================
 
! --------------------------------------------------------------------------
! Timer (Intel Fortran)
! --------------------------------------------------------------------------
        real time_begin
        real time_end
        integer :: itime_step
        logical :: finish_qmmm=.true.

! Procedure
! ===========================================================================

! ***************************************************************************
!                           MPI Specific Stuff
! ***************************************************************************
        call init_MPI (iammaster, iammpi, my_proc, nprocs)
        wrtout = .true.
        if (.not. iammaster) then
          write (*,*) '------ SLAVE:  CALL KSPACE_SLAVE ------'
          call kspace_slave (nprocs, my_proc)
        end if
        print *,'PROC',nprocs, my_proc
        call connect()
         print*,' rank=',rank
        if(rank==0)then
         print *,'first time : Client should be created fireball.in,'
         print *,' Fdata.optional and input.bas'
         print *,'needed to initbasis() and readdata()'
         print *,'natoms not change from input.bas'

         call MPI_RECV(mpi,1,MPI_INTEGER,  0, 0,intercomm, status, ierr_mpi) 
         call cpu_time (time_begin)
        
         call initbasics ()    
         call readdata ()
 
         itime_step = 1
        
        do while (finish_qmmm)

            call MPI_RECV(ratom,3*natoms, MPI_DOUBLE_PRECISION,0,0,intercomm,status,ierr_mpi)
            print*,'ratoms, ' ,ratom
            if(ratom(1,1).ne.null) then
              call MPI_RECV(qmmm_struct%qm_mm_pairs,1, MPI_INTEGER,0,0, intercomm,status,ierr_mpi)
              allocate(qmmm_struct%qm_xcrd(4,qmmm_struct%qm_mm_pairs))
              allocate(qmmm_struct%dxyzcl(3,qmmm_struct%qm_mm_pairs))
              call MPI_RECV(qmmm_struct%qm_xcrd,4*qmmm_struct%qm_mm_pairs,MPI_DOUBLE_PRECISION,0,0,intercomm,status,ierr_mpi)



              call scf_loop (itime_step)
              call postscf ()
              call getenergy (itime_step)
              call getforces (itime_step)

              call MPI_SEND((etot*23.061d0),1, MPI_DOUBLE_PRECISION,0,0,intercomm,ierr_mpi)
              call MPI_SEND(-ftot*23.061d0,3*natoms, MPI_DOUBLE_PRECISION,0,0,intercomm,ierr_mpi)
              call MPI_SEND(qmmm_struct%dxyzcl,3*qmmm_struct%qm_mm_pairs, MPI_DOUBLE_PRECISION,0,0,intercomm,ierr_mpi)
              itime_step = itime_step +1
           else
              finish_qmmm=.false.
           end if

        end do

         
         call MPI_UNPUBLISH_NAME('fireball_server', MPI_INFO_NULL, port_name, ierr_mpi)
         call MPI_CLOSE_PORT(port_name, ierr_mpi)
         end if 

         print*,'mpi_f2' 
  !    call MPI_FINALIZE(ierr_mpi)
         print*,'mpi_f3' 
 
! timer

! Format Statements
! ===========================================================================
100     format (2x, 70('='))
 
        stop
      end program fireball_server
