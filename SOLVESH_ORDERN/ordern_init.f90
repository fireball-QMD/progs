! copyright info:
!
!                             @Copyright 2001
!                            Fireball Committee
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

 
! ordern_init.f90
! Program Description
! ===========================================================================
!       This routine initializes a bunch of stuff necessary for running the 
! ordern - parallel version of the code. 
!
! ===========================================================================
! Original Order-N compilation written by Spencer Shellman.

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
        subroutine ordern_init (nprocs, my_proc,           &
     &                          itheory, itheory_xc, icluster, a1vec, a2vec, &
     &                          a3vec, ifixcharge, ifixneigh, iqout,         &
     &                          itestrange, testrange, rcutoff, nstepi,      &
     &                          max_scf_iterations, igauss, ivdw)
        use charges
        use interactions
        use dimensions
        use neighbor_map
        use configuration
        use forces
        use mpi_declarations
        use ordern
        implicit none

        include 'mpif.h'

! The arguments should be passed in from the main routine.
! The purpose of ordern_init is to distribute some global variables across all
! processors so that they can run some routines in parallel that depend on the
! variables.
 
! Argument Declaration and Description
! ===========================================================================
! Input 
        integer, intent (in) :: icluster
        integer, intent (in) :: ifixcharge
        integer, intent (in) :: ifixneigh
        integer, intent (in) :: iqout
        integer, intent (in) :: itestrange
        integer, intent (in) :: itheory
        integer, intent (in) :: itheory_xc
        integer, intent (in) :: max_scf_iterations
        integer, intent (in) :: my_proc
        integer, intent (inout) :: nprocs
        integer, intent (in) :: nstepi
        integer, intent (in) :: igauss, ivdw
 
        real, dimension (3), intent (in) :: a1vec, a2vec, a3vec
        real, intent (in), dimension (nspec_max, nsh_max) :: rcutoff
        real, intent (in) :: testrange

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ibsize
        integer irbsize
        integer ierror
        integer iproc
        integer isendrows
        integer nprocs_max
        integer packsize,packpos

        integer, dimension (:), allocatable :: ibcastbuffer
 
        integer*1, dimension (:), allocatable :: ibcb
        integer*1, dimension (:), allocatable :: irbcastbuffer

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! The actual number of processors that will work on the matrix.
! For best results the number of actual processors used should be something
! that evenly divides the Hamiltonian matrix into equal sized blocks. 
        do nprocs_max = natoms, 1, -1
         if (mod(natoms,nprocs_max) .eq. 0) exit
        end do
        nactualprocs = nprocs_max
        do nprocs_max = nprocs, 1, -1
         if (mod(norbitals,nprocs_max) .eq. 0) exit
        end do
        if (nactualprocs .gt. nprocs_max) nactualprocs = nprocs_max 
        nprocs = nactualprocs

! Find the maximum number of rows that would exist over all of the processors.
! This is stored in ncrowsmax and nprowsmax.
        nprowsmax = norbitals/nactualprocs
        if (mod(norbitals,nactualprocs) .gt. 0) nprowsmax = nprowsmax + 1
        do iproc = 1, nactualprocs - 1
         isendrows = norbitals/nactualprocs
         if (iproc .lt. mod(norbitals,nactualprocs)) isendrows = isendrows + 1 
         if (nprowsmax .lt. isendrows) nprowsmax = isendrows
        end do

        ibsize = 19
        allocate (ibcastbuffer(ibsize))

! Broadcast the number of orbitals, etc. to all processors.
! The information will be received in kspace_ordern_slave.
        ibcastbuffer(1) = numorb_max
        ibcastbuffer(2) = natoms
        ibcastbuffer(3) = nspecies
        ibcastbuffer(4) = norbitals
        ibcastbuffer(5) = nprowsmax
        ibcastbuffer(6) = nstepi
        ibcastbuffer(7) = itheory
        ibcastbuffer(8) = itheory_xc
        ibcastbuffer(9) = icluster
        ibcastbuffer(10) = ifixcharge
        ibcastbuffer(11) = ifixneigh
        ibcastbuffer(12) = iqout
        ibcastbuffer(13) = isorpmax
        ibcastbuffer(14) = ideriv_max
        ibcastbuffer(15) = itestrange
        ibcastbuffer(16) = max_scf_iterations
        ibcastbuffer(17) = igauss
        ibcastbuffer(18) = ivdw
        ibcastbuffer(19) = interactions2c_max

        call MPI_BCAST (ibcastbuffer, ibsize, MPI_INTEGER, 0, MPI_COMM_WORLD,&
    &                   ierror)
        deallocate (ibcastbuffer)

! Create the communication domain consisting of relevant processors.
        call MPI_COMM_SPLIT (MPI_COMM_WORLD, min(my_proc/nactualprocs,1),    &
    &                        my_proc, MPI_BTN_WORLD, ierror)

        irbsize = 11 + 3*natoms + 3*125 + nspec_max*nsh_max
        call MPI_PACK_SIZE(irbsize,mpi_whatever_real,MPI_BTN_WORLD,packsize,ierror)
        allocate (irbcastbuffer(packsize))

! Send more stuff to all processors.
        packpos = 0
        call MPI_PACK(a1vec,3,mpi_whatever_real,irbcastbuffer,packsize,packpos,MPI_BTN_WORLD,ierror)
        call MPI_PACK(a2vec,3,mpi_whatever_real,irbcastbuffer,packsize,packpos,MPI_BTN_WORLD,ierror)
        call MPI_PACK(a3vec,3,mpi_whatever_real,irbcastbuffer,packsize,packpos,MPI_BTN_WORLD,ierror)
        call MPI_PACK(testrange,1,mpi_whatever_real,irbcastbuffer,packsize,packpos,MPI_BTN_WORLD,ierror)
        call MPI_PACK(range_vdw,1,mpi_whatever_real,irbcastbuffer,packsize,packpos,MPI_BTN_WORLD,ierror)
        call MPI_PACK(ratom,3*natoms,mpi_whatever_real,irbcastbuffer,packsize,packpos,MPI_BTN_WORLD,ierror)
        call MPI_PACK(xl,3*125,mpi_whatever_real,irbcastbuffer,packsize,packpos,MPI_BTN_WORLD,ierror)
        call MPI_PACK(rcutoff,nspec_max*nsh_max,mpi_whatever_real,irbcastbuffer,packsize,packpos,MPI_BTN_WORLD,ierror)

        call MPI_BCAST(packpos,1,MPI_INTEGER,0,MPI_BTN_WORLD,ierror)
        call MPI_BCAST (irbcastbuffer, packpos, MPI_PACKED, 0, MPI_BTN_WORLD, ierror)
        deallocate (irbcastbuffer)

        ibsize = natoms + nsh_max*nspecies + natoms + nspecies + nspecies    &
     &          + natoms + nspecies
        call MPI_PACK_SIZE(ibsize,MPI_INTEGER,MPI_BTN_WORLD,packsize,ierror)
        allocate (ibcb(packsize))

        packpos = 0
        call MPI_PACK(imass,natoms,MPI_INTEGER,ibcb,packsize,packpos,MPI_BTN_WORLD,ierror)
        call MPI_PACK(nssh,nspecies,MPI_INTEGER,ibcb,packsize,packpos,MPI_BTN_WORLD,ierror)
        call MPI_PACK(degelec,natoms,MPI_INTEGER,ibcb,packsize,packpos,MPI_BTN_WORLD,ierror)
        call MPI_PACK(lssh,nsh_max*nspecies,MPI_INTEGER,ibcb,packsize,packpos,MPI_BTN_WORLD,ierror)
        call MPI_PACK(nelectron,natoms,MPI_INTEGER,ibcb,packsize,packpos,MPI_BTN_WORLD,ierror)
        call MPI_PACK(num_orb,nspecies,MPI_INTEGER,ibcb,packsize,packpos,MPI_BTN_WORLD,ierror)
        call MPI_PACK(nzx,nspecies,MPI_INTEGER,ibcb,packsize,packpos,MPI_BTN_WORLD,ierror)

        call MPI_BCAST(packpos,1,MPI_INTEGER,0,MPI_BTN_WORLD,ierror)
        call MPI_BCAST (ibcb, packpos, MPI_PACKED, 0, MPI_BTN_WORLD, ierror)
        deallocate (ibcb)

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end subroutine ordern_init
