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

 
! kspace_ordern_init.f90
! Program Description
! ===========================================================================
!       This routine initializes a bunch of stuff before kspace_ordern.f90 
! is called. 
!
! ===========================================================================
! Original Order-N compilation by Spencer Shellman

! Code rewritten by:
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

! simple bcast function so that we need only one dummy function in kspace_fk

        subroutine bcast (buf, len, type)
          use mpi_declarations
          use ordern
          implicit none

          include 'mpif.h'

          integer*1, dimension (*), intent (inout) :: buf
          integer, intent (in) :: len, type

          integer ierror

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

          call MPI_BCAST (buf, len, type, 0, MPI_BTN_WORLD, ierror)

          return
        end subroutine bcast



! functions that communicate module variables to all processors

        subroutine Qin_bcast (natoms)
          use mpi_declarations
          use ordern
          use charges
          use interactions
          implicit none

          integer, intent (in) :: natoms

          call bcast (Qin, nsh_max*natoms, mpi_whatever_real)

          return
        end subroutine Qin_bcast

        subroutine Qneutral_bcast (nspecies)
          use mpi_declarations
          use ordern
          use charges
          use interactions
          implicit none

          integer, intent (in) :: nspecies

          call bcast (Qneutral, nsh_max*nspecies, mpi_whatever_real)

          return
        end subroutine Qneutral_bcast

        subroutine ewald_bcast (natoms, itheory)
          use mpi_declarations
          use ordern
          use interactions
          use neighbor_map
          implicit none

          integer, intent (in) :: natoms, itheory

          if (itheory .eq. 1 .or. itheory .eq. 2) then 
             call bcast (ewaldlr, (numorb_max**2)*natoms*neigh_max, mpi_whatever_real)
             call bcast (ewaldsr, (numorb_max**2)*natoms*neigh_max, mpi_whatever_real)
             call bcast (vca, (numorb_max**2)*natoms*neigh_max, mpi_whatever_real)
             call bcast (vxc_ca, (numorb_max**2)*natoms*neigh_max, mpi_whatever_real)
          end if

          return
        end subroutine ewald_bcast

        subroutine scf_bcast (scf_achieved)
          include 'mpif.h'

          logical, intent(in) :: scf_achieved

          integer ierror

          call bcast (scf_achieved, 1, MPI_LOGICAL)

          return
        end subroutine scf_bcast


        subroutine ME_max_bcast
          use interactions
          include 'mpif.h'

          integer ierror

          call bcast (ME2c_max, 1, MPI_INTEGER)
          call bcast (ME2cPP_max, 1, MPI_INTEGER)
          call bcast (ME3c_max, 1, MPI_INTEGER)

          return
        end subroutine ME_max_bcast


        subroutine ratom_bcast (natoms, ratom)
          use mpi_declarations
          use ordern
          implicit none

          integer, intent (in) :: natoms
          real, intent (in), dimension (3, natoms) :: ratom

          integer ierror

          call bcast (ratom, 3*natoms, mpi_whatever_real)

          return
        end

        subroutine kspace_ordern_init (natoms, nspecies, nprocs, my_proc,    &
     &                                 Kscf, itime_step, ratom, ztot,        &
     &                                 itheory, icluster, iforce)
        use charges
        use interactions
        use neighbor_map
        use forces
        use ordern
        implicit none

        include 'mpif.h'
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itime_step
        integer, intent (in) :: Kscf
        integer, intent (in) :: natoms
        integer, intent (in) :: nspecies
        integer, intent (in) :: nprocs
        integer, intent (in) :: my_proc
        integer, intent (in) :: itheory, icluster
        integer, intent (in) :: iforce

        real, intent (in) :: ztot

        real, intent (in), dimension (3, natoms) :: ratom

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ibsize
        integer ierror
        integer iproc
        integer irbsize
        integer isendrows
        integer pksize,pkpos
        integer ioptionlwf

        integer, dimension (:), allocatable :: ibcastbuffer

        integer*1, dimension (:), allocatable :: irbcb

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ===========================================================================
        ioptionlwf = 1

! Initialize the dimensions for the sparse calculations.
        if (Kscf .eq. 1) then
         call set_dimensions (natoms, ioptionlwf, ratom, ztot)

! Now that nbands is set...
         ncrowsmax = nbands/nactualprocs
         if (mod(nbands,nactualprocs) .gt. 0) ncrowsmax = ncrowsmax + 1
         do iproc = 1, nactualprocs - 1
          isendrows = nbands/nactualprocs
          if (iproc .lt. mod(nbands,nactualprocs)) isendrows = isendrows + 1 
          if (ncrowsmax .lt. isendrows) ncrowsmax = isendrows
         end do
        end if

        ibsize = 9
        allocate (ibcastbuffer(ibsize))

        ibcastbuffer(1) = itime_step
        ibcastbuffer(2) = nbands 
        ibcastbuffer(3) = ncrowsmax 
        ibcastbuffer(4) = nhmax 
        ibcastbuffer(5) = ncmax 
        ibcastbuffer(6) = nctmax 
        ibcastbuffer(7) = nFmax 
        ibcastbuffer(8) = nFtmax 
        ibcastbuffer(9) = iforce
        call MPI_BCAST (ibcastbuffer, ibsize, MPI_INTEGER, 0, MPI_BTN_WORLD, &
     &                  ierror)
        deallocate (ibcastbuffer)
  
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================

        return
      end subroutine kspace_ordern_init
