! copyright info:
!
!                             @Copyright 2005
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

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! diag_k.f90
! Program Description
! ===========================================================================
!       This is  a subroutine of k-loop
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
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine kspace_slave (nprocs, my_proc)

        use configuration  
        use density
        use dimensions
        use interactions
        use neighbor_map
        use kpoints
        use mpi_declarations
        implicit none

        include 'mpif.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: nprocs
        integer, intent (in) :: my_proc


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
! Map of global variables
        integer natoms
        integer iqout
        integer Kscf
        integer icluster
        integer iwrteigen
        integer iwrtdos
        integer iwrtatom
        integer iwrthop
        integer norbitals_new

! Local variables
        integer ikpoint
        integer imu
        integer istart
        integer nkpts
        integer ierror

        real, dimension (3) :: k_temp

! Procedure
! ===========================================================================
!        write (*,*) 'SLAVE: Call Bcast_Data_0_slave'
        call Bcast_data_0_slave (natoms, iqout, Kscf)
        icluster = 0
        iwrteigen = 0
        iwrtdos = 0
        iwrtatom = 0
        iwrthop = 0
        iwrteigen = 0

! Slave trapped in infinity loop
        do while (1) 

          blowre = 0.0d0
          blowim = 0.0d0
          eigen_k = 0.0d0
          bbnkim = 0.0d0
          bbnkre = 0.0d0

! Send Kscf status
          call MPI_BCAST (Kscf, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
          if (Kscf .eq. 1) then     
!            write (*,*) 'SLAVE: Call Bcast_Data_1_slave'
            call Bcast_data_1_slave (natoms)
          else 
!            write (*,*) 'SLAVE: Call Bcast_Data_2_slave'
            call Bcast_data_2_slave (natoms)
          endif
! define rank of the k-loop
          call MPI_COMM_RANK (MPI_COMM_WORLD, my_proc, ierror)
          nkpts = nkpoints/nprocs
          if (my_proc .lt. mod(nkpoints,nprocs)) then
          nkpts = nkpts + 1
          istart = nkpts*my_proc + 1
         else
          istart = (nkpts + 1)*mod(nkpoints,nprocs)                 &
                      + nkpts*(my_proc - mod(nkpoints,nprocs)) + 1
         end if
          do ikpoint = istart, istart - 1 + nkpts

! The subroutine kspace wants the k-vector in inverse angstrom units.
! NOT pi/alat units.
           k_temp(:) = special_k(:,ikpoint)

           call kspace (natoms, nprocs, my_proc, Kscf, iqout, icluster,  &
     &                iwrteigen, ikpoint, k_temp,         &
     &                nkpoints, iwrtdos, iwrthop, iwrtatom, itrans)


          end do ! do ikpoint
! Gather Data
        call Gather_k (natoms, iqout, norbitals)
        end do ! end while


! Format Statements
! ===========================================================================

        return
      end subroutine kspace_slave
