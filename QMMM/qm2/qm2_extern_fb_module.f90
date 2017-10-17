#include "../include/dprec.fh"
module qm2_extern_fb_module
! ----------------------------------------------------------------
! Interface for TeraChem based QM and QM/MM MD 
!
! Currently supports:
! pure QM
! QM/MM with cutoff for QM-MM electrostatics under periodic
! boundary conditions
!
! Author: Andreas Goetz (agoetz@sdsc.edu)
!
! Date: November 2010
!
! ----------------------------------------------------------------

  use qm2_extern_util_module, only: debug_enter_function, debug_exit_function

  implicit none

  private
  public :: get_fb_forces, fb_finalize
  logical, save :: do_mpi = .false.  ! Used in finalize subroutine

  character(len=*), parameter :: module_name = "qm2_extern_fb_module"

  type fb_nml_type
        ! Poner los input de fireball
     character(len=20) :: executable
     integer :: debug
     integer :: mpi
  end type fb_nml_type

  integer, save         :: newcomm ! Initialized in mpi_init subroutine

contains

  ! --------------------------------------
  ! Get QM energy and forces from TeraChem
  ! --------------------------------------
  subroutine get_fb_forces( do_grad, nstep, ntpr_default, id, nqmatoms, qmcoords,&
       qmtypes, nclatoms, clcoords, escf, dxyzqm, dxyzcl, charge, spinmult )

    use qm2_extern_util_module, only: print_results, check_installation, write_dipole, write_charges
    use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS, ZERO

    logical, intent(in) :: do_grad              ! Return gradient/not
    integer, intent(in) :: nstep                ! MD step number
    integer, intent(in) :: ntpr_default         ! frequency of printing
    character(len=3), intent(in) :: id          ! ID number for PIMD or REMD
    integer, intent(in) :: nqmatoms             ! Number of QM atoms
    _REAL_,  intent(in) :: qmcoords(3,nqmatoms) ! QM coordinates
    integer, intent(in) :: qmtypes(nqmatoms)    ! QM atom types (nuclear charge in au)
    integer, intent(in) :: nclatoms             ! Number of MM atoms
    _REAL_,  intent(in) :: clcoords(4,nclatoms) ! MM coordinates and charges in au
    _REAL_, intent(out) :: escf                 ! SCF energy
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)   ! SCF QM force
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)   ! SCF MM force
    integer, intent(in) :: charge, spinmult     ! Charge and spin multiplicity

    _REAL_              :: dipmom(4,3)          ! Dipole moment {x, y, z, |D|}, {QM, MM, TOT}
    _REAL_              :: qmcharges(nqmatoms)  ! QM charges from population analysis

    type(fb_nml_type), save :: fb_nml
    logical, save :: first_call = .true.
    integer :: i
    integer :: printed =-1 ! Used to tell if we have printed this step yet 
                           ! since the same step may be called multiple times
    character(len=150) :: call_buffer

    ! assemble input - / output data filenames

    
    ! Setup on first program call
    if ( first_call ) then
      first_call = .false.
      write (6,*) '>>> Running QM calculation with Fireball  <<<'
      write(6,*) 'call get_namelist( ntpr_default, fb_nml )',ntpr_default,fb_nml
      call get_namelist( ntpr_default, fb_nml )
      call check_installation( trim(fb_nml%executable), id, .true., fb_nml%debug )
     ! call print_namelist(fb_nml)
    end if

      write (6,*)'3dgt fb_nml%mpi,nqmatoms,,',fb_nml%mpi,nqmatoms, qmcoords, &
         qmtypes, nclatoms, clcoords, &
        fb_nml, escf, dxyzqm, dxyzcl, dipmom, qmcharges, do_grad, id, charge,&
        spinmult
#ifdef MPI
# ifndef MPI_1
    if (fb_nml%mpi==1 ) then ! Do mpi (forced to 0 ifndef MPI)
      print *,'mpi_hook'
      call mpi_hook( nqmatoms, qmcoords, qmtypes, nclatoms, clcoords,&
        fb_nml, escf, dxyzqm, dxyzcl, dipmom, qmcharges, do_grad, id, charge, spinmult )
      print *,'mpi_hook'
    else
      write (6,'(/,a,/)')'1dgt'
# else
    ! If we are using MPI 1.x the code will not compile since
    ! MPI_LOOKUP_NAME is part of the MPI 2 standard, so  just quit
    if (fb_nml%mpi==1 ) then 
    call sander_bomb('(qm2_extern_fb_module)', &
      '&unsupported MPI version', &
      'Will quit now.')
    else
# endif
#endif
      


#ifdef MPI
    end if
#endif

    if ( do_grad ) then
       ! Convert Hartree/Bohr -> kcal/(mol*A)
       dxyzqm(:,:) = dxyzqm(:,:) * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
       if ( nclatoms > 0 ) then
          dxyzcl(:,:) = dxyzcl(:,:) * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
       end if
    else
       dxyzqm = ZERO
       if ( nclatoms > 0 ) dxyzcl = ZERO
    end if

    escf = escf * CODATA08_AU_TO_KCAL

    !call print_results( 'qm2_extern_fb_module', escf, nqmatoms, dxyzqm, nclatoms, dxyzcl )


  end subroutine get_fb_forces


  subroutine get_namelist( ntpr_default, fb_nml)

    use UtilitiesModule, only: Upcase
    implicit none
    character(len=20) :: executable
    integer, intent(in) :: ntpr_default
    type(fb_nml_type), intent(out) :: fb_nml
    integer :: mpi
    namelist /fb/ mpi
    integer :: i, ierr
     executable      = 'fireball_server'
   
    mpi          = 1 ! Default to using MPI if available

    ! Read namelist
    rewind 5
    read(5,nml=fb,iostat=ierr)

#ifndef MPI
        if ( fb_nml%mpi == 1 ) then
          write(6,'(a)') '| Warning: mpi=1 selected but sander was not compiled with MPI support.'
          write(6,'(a)') '| Continuing with mpi=0'
        end if
        fb_nml%mpi         = 0 ! Can't pick MPI if not available 
#else
        fb_nml%mpi         = mpi
#endif


        fb_nml%executable   = executable
        fb_nml%debug        =  0
    ! Need this variable so we don't call MPI_Send in the finalize subroutine
    if (mpi==1 ) then
      do_mpi=.true.
    end if

  end subroutine get_namelist


#if defined(MPI) && !defined(MPI_1)
  ! Perform MPI communications with terachem. Requires MPI 2.0 or above to use
  subroutine mpi_hook( nqmatoms, qmcoords, qmtypes, nclatoms, clcoords,&
       fb_nml, escf, dxyzqm, dxyzcl, dipmom, qmcharges, do_grad, id, charge, spinmult )
    
    use ElementOrbitalIndex, only : elementSymbol
    use qm2_extern_util_module, only: check_installation    
    
    implicit none
    include 'mpif.h'

    integer, intent(in) :: nqmatoms
    _REAL_,  intent(in) :: qmcoords(3,nqmatoms) 
    integer, intent(in) :: qmtypes(nqmatoms)
    integer, intent(in) :: nclatoms
    _REAL_,  intent(in) :: clcoords(4,nclatoms)
    type(fb_nml_type), intent(in) :: fb_nml
    _REAL_, intent(out) :: escf
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)
    _REAL_, intent(out) :: dipmom(4,3)
    _REAL_, intent(out) :: qmcharges(nqmatoms)
    logical, intent(in) :: do_grad
    character(len=3), intent(in) :: id
    integer         , intent(in) :: charge, spinmult

    character(len=2)    :: atom_types(nqmatoms)
    _REAL_              :: coords(3,nqmatoms+nclatoms)
    _REAL_              :: charges(nclatoms)
    _REAL_              :: dxyz_all(3,nclatoms+nqmatoms)

    logical,save        :: first_call=.true.
    integer             :: i, status(MPI_STATUS_SIZE)
    integer             :: ierr

    integer :: system
    integer :: stat

    call debug_enter_function( 'mpi_hook', module_name, fb_nml%debug )


    ! ---------------------------------------------------
    ! Initialization: Connect to "terachem_port", set    
    ! newcomm (global), send relevant namelist variables.
    ! ---------------------------------------------------

    write(6,*)'dgt connect', first_call


    if (first_call) then 
      first_call=.false.
      write (6,*) 'connect_to_fireball only first time'

      write (6,'(/,a,/)') '   >>> Running QM calculation with Fireball <<<'
      !call get_namelist( ntpr_default, fb_nml )
!      call check_installation( trim(fb_nml%executable), id, .true., fb_nml%debug)
      !call print_namelist(fb_nml)
      print *,'mpiexec -n 1 -env I_MPI_DEBUG 2 -env I_MPI_DEVICE sock ./fireball_server'
!      stat = system('mpiexec -n 1 -env I_MPI_DEBUG 2 -env I_MPI_DEVICE sock ./fireball_server &')
      print *,'jesus'
      call connect_to_fireball( fb_nml, nqmatoms, atom_types, do_grad, id, charge, spinmult )
    end if

       print *,'send'
       !call MPI_Send( nqmatoms, 1, MPI_INTEGER, 0, 0, newcomm, ierr )
       !write(*,*) 'qmcoords'
       !write(*,*) qmcoords
       call MPI_SEND(qmcoords,3*nqmatoms, MPI_DOUBLE_PRECISION, 0, 0, newcomm, ierr) 
       call MPI_SEND(nclatoms, 1, MPI_INTEGER, 0, 0, newcomm, ierr )
       call MPI_SEND(clcoords,4*nclatoms, MPI_DOUBLE_PRECISION, 0, 0, newcomm,ierr)

print *,'qmcoords,nclatoms,clcoords', qmcoords,nclatoms,clcoords
       !call MPI_Send( qmcoords, 3*nqmatoms, MPI_REAL, 0, 2, newcomm, ierr )
       !call MPI_SEND(qmcoords,3*nqmatoms, MPI_REAL, 0, 0, newcomm, ierr)
       !call MPI_SEND(S,2*natoms, MPI_CHARACTER, 0, 0, intercomm, ierr)
       !call MPI_RECV(qmcoords,3*nqmatoms, MPI_DOUBLE_PRECISION, 0, 0, newcomm, status,ierr)
       call MPI_Recv(escf, 1, MPI_DOUBLE_PRECISION, 0, 0, newcomm, status, ierr )
       call MPI_RECV(dxyzqm,3*nqmatoms, MPI_DOUBLE_PRECISION, 0, 0, newcomm,status,ierr)
       call MPI_RECV(dxyzcl,3*nclatoms, MPI_DOUBLE_PRECISION, 0, 0, newcomm,status,ierr)

print *,'escf,dxyzqm,dxyzcl,',escf,dxyzqm,dxyzcl
       !call MPI_RECV(qmcharges,nqmatoms, MPI_DOUBLE_PRECISION, 0, 0, newcomm,status,ierr)
       !write(*,*) 'qmcoords'
       !write(*,*) qmcoords
       !call MPI_RECV(S,2*natoms, MPI_CHARACTER, 0, 0, intercomm,status,ierr)
       !print *,'rec'


    !call MPI_SEND(R,3*natoms, MPI_REAL, 0, 0, intercomm, ierr)
!    call MPI_Send( qmcoords, 3*nqmatoms, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr ) 

    !call MPI_Send( nclatoms, 1, MPI_INTEGER, 0, 2, newcomm, ierr ) 


    !call MPI_Send( charges, nclatoms, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr ) 
    ! call MPI_Send( coords, 3*nclatoms, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr ) 

    !call MPI_Recv( escf, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
    !call MPI_Recv( qmcharges(:), nqmatoms, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
   ! if ( do_grad ) then
   !   call MPI_Recv( dxyz_all, 3*(nqmatoms+nclatoms), MPI_DOUBLE_PRECISION, &
   !         MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
       
       ! Poplulate our output arrays with gradients from terachem
    !   do i=1, nqmatoms
    !      dxyzqm(:,i)=dxyz_all(:,i)
    !   end do
    !   do i=1, nclatoms
    !      dxyzcl(:,i)=dxyz_all(:,i+nqmatoms)
    !   end do

    !end if

!    call MPI_RECV(iR,3*natoms, MPI_REAL, 0, 0, intercomm, status,ierr)
!    call MPI_Recv( qmcoords, 3*nqmatoms, MPI_DOUBLE_PRECISION, 0, 2, newcomm,status, ierr ) 

    do i=1,nqmatoms
     write (*,'(2x, 3(2x,f12.6))') qmcoords(:,i)
    end do



  end subroutine mpi_hook

  ! -------------------------------------------------
  ! Search for name published by TeraChem and connect
  ! (this step initializes newcomm)
  ! Send relevant namelist variables to terachem
  ! -------------------------------------------------
  subroutine connect_to_fireball( fb_nml, nqmatoms, atom_types, do_grad, id, charge, spinmult )

    implicit none
    include 'mpif.h'
    
    type(fb_nml_type), intent(in) :: fb_nml
    integer          , intent(in) :: nqmatoms
    character(len=2) , intent(in) :: atom_types(nqmatoms)
    logical          , intent(in) :: do_grad
    character(len=3) , intent(in) :: id
    integer          , intent(in) :: charge, spinmult

    character(len=17) :: server_name="fireball_server"
    integer, parameter  :: clen=128 ! Length of character strings we are using
    character(255) :: port_name
    character(len=clen) :: dbuffer(2,32)
    _REAL_          :: timer
    integer         :: ierr, i, j, irow
    logical         :: done=.false.
 
    integer rank,  intercomm, status, size, mp


    !call debug_enter_function( 'connect_to_fireball', module_name, fb_nml%debug )

    ! -----------------------------------
    ! Look for server_name, get port name
    ! After 60 seconds, exit if not found
    ! -----------------------------------

    print *,'fb_nml, nqmatoms, atom_types, do_grad, id, charge, spinmult', fb_nml, nqmatoms, atom_types, do_grad, id, charge, spinmult
    print *,'server_name=',server_name
    print *,'id =', id 
    if ( trim(id) /= '' ) then
      server_name = trim(server_name)//'.'//trim(id)
    end if
    !if ( fb_nml%debug > 1 ) then
    !  write(6,'(2a)') 'Looking up server under name:', trim(server_name)
    !  call flush(6)
    !end if
    timer = MPI_WTIME(ierr)
    do while (done .eqv. .false.)
   !port_name='4160094208.0;tcp://10.1.255.253:54296+4160094209.0;tcp://10.1.255.253:43662:300'
   

 !   print *, '-------------------'
 !      server_name = 'fireball_server'
 !      call MPI_INIT ( ierr )
 !      print *,'ierr =',ierr
 !      call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
 !      call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
 !      print*, 'rank', rank, ': size ',size
 !      if(rank == 0) then
 !       call MPI_LOOKUP_NAME (server_name, MPI_INFO_NULL, port_name, ierr)
 !       call MPI_COMM_CONNECT(port_name, MPI_INFO_NULL,0,MPI_COMM_SELF,intercomm, ierr)
 !      end if



 
    print *,'server_name', server_name, 'port_name', port_name 
    call MPI_LOOKUP_NAME(trim(server_name), MPI_INFO_NULL, port_name, ierr)
      if (ierr == MPI_SUCCESS) then
        if ( fb_nml%debug > 1 ) then
          write(6,'(2a)') 'Found port: ', trim(port_name)
          call flush(6)
        end if
        done=.true.

      end if

      if ( (MPI_WTIME(ierr)-timer) > 60 ) then ! Time out after 60 seconds
        call sander_bomb('connect_to_fireball() ('//module_name//')', &
          '"'//trim(server_name)//'" not found. Timed out after 60 seconds.', &
          'Will quit now')
      end if

    end do

    ! ----------------------------------------
    ! Establish new communicator via port name
    print *,'test'
    ! ----------------------------------------


       call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
       print*, 'rank', rank, ': size ',size


    print *,'port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, newcomm, ierr', port_name,MPI_INFO_NULL, 0, MPI_COMM_SELF, newcomm, ierr

    call MPI_COMM_CONNECT(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, newcomm, ierr)

    print *,'test'
  !  if ( fb_nml%debug > 1 ) then
  !    write(6,'(a,i0)') 'Established new communicator:', newcomm
    !print *,'test'

!    call SendRecv(natoms,ratom)


   !   call flush(6)
   ! end if

    print *,'test'
    mp=0
    print *,'(mp,1, MPI_INTEGER, 0, 0, intercomm, ierr',mp, MPI_INTEGER,  intercomm, ierr
    call MPI_SEND(mp,1, MPI_INTEGER, 0, 0, newcomm, ierr)



  end subroutine connect_to_fireball

#endif


  subroutine fb_finalize()

    implicit none
#ifdef MPI
    include 'mpif.h'

    integer :: ierr
    _REAL_  :: empty
    if (do_mpi) then
      call MPI_Send( empty, 1, MPI_DOUBLE_PRECISION, 0, 0, newcomm, ierr )
    end if
#endif

  end subroutine fb_finalize
      

end module qm2_extern_fb_module
