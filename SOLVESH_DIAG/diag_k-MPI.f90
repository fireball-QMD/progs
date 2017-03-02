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
        subroutine diag_k (natoms, nprocs, my_proc, Kscf, iqout, icluster,  &
     &                     iwrteigen, iwrtdos, iwrthop, iwrtatom)

        use configuration
        use mpi_declarations  
        use density
        use dimensions
        use interactions
        use neighbor_map
        use kpoints
        implicit none

        include 'mpif.h'

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icluster
        integer, intent (in) :: iqout
        integer, intent (in) :: iwrteigen
        integer, intent (in) :: Kscf
        integer, intent (in) :: natoms
        integer, intent (in) :: nprocs
        integer, intent (in) :: my_proc
! CGP
        integer, intent (in) :: iwrtdos
        integer, intent (in) :: iwrthop
        integer, intent (in) :: iwrtatom
! end CGP
! MPI
        integer, dimension (5) :: buffer_int

! Output

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ikpoint
        integer imu
        integer istart
        integer nkpts
        integer ierror

        real, dimension (3) :: k_temp

! Procedure
! ===========================================================================
! Init MPI; broadcast variables
        if (mpi_on .eq. 0) then
!          write (*,*) 'MASTER: Call Bcast_Data_0' 
          call Bcast_data_0 (natoms, iqout, Kscf)
          mpi_on = 1
        endif

! Send Kscf status 
        call MPI_BCAST (Kscf, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
        if (Kscf .eq. 1) then 
!          write (*,*) 'MASTER: Call Bcast_data_1'
          call Bcast_data_1 (natoms)
        else
!          write (*,*) 'MASTER: Call Bcast_data_2'
          call Bcast_data_2 (natoms)
        endif

! Now we have the real space hamiltonian and overlap. Compute the k-space
! Hamiltonian and overlap and diagonalize. First put the k-points into the

! The eigenvectors of the overlap are always computed.
! Allocate the arrays for wavefunction coefficients
        if (iqout .ne. 2) then 
          allocate (blowre (norbitals, norbitals, nkpoints))
          allocate (blowim (norbitals, norbitals, nkpoints))
          blowre = 0.0d0
          blowim = 0.0d0
        endif
        allocate (bbnkre (norbitals, norbitals, nkpoints))
        allocate (bbnkim (norbitals, norbitals, nkpoints))
        allocate (eigen_k (norbitals, nkpoints))
        eigen_k = 0.0d0
        bbnkim = 0.0d0
        bbnkre = 0.0d0

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
     &                iwrteigen, ikpoint, k_temp,        &
     &                nkpoints, iwrtdos, iwrthop, iwrtatom, itrans)


        end do ! do ikpoint

! Gather data
        call Gather_k (natoms, iqout, norbitals)

! Format Statements
! ===========================================================================
100     format (2x, 4(2x,f11.5))
101     format (i4, 8f11.5)


        return
      end subroutine diag_k
