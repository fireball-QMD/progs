#include "../include/dprec.fh"
module qm2_extern_fb_module
! ----------------------------------------------------------------
! Interface for Fireball based QM and QM/MM MD 
!
! Currently supports:
! pure QM
! QM/MM with cutoff for QM-MM electrostatics under periodic
! boundary conditions
!
! Author: Jesus I. Mendieta-Moreno (jesus.mendieta@uam.es)
!
! Date: November 2010
!
! ----------------------------------------------------------------

  use qm2_extern_util_module, only: debug_enter_function, debug_exit_function

  implicit none

  private
  public :: get_fb_forces

  character(len=*), parameter :: module_name = "qm2_extern_fb_module"

  type fb_nml_type
        ! Poner los input de fireball
     character(len=20) :: executable
     integer :: max_scf_iterations
     integer :: qmmm_int
     integer :: idftd3
     integer :: idipole
     integer :: debug
     integer :: mpi
     integer :: iqout
     integer :: imcweda
     integer :: ihorsfield
     integer :: iensemble
     integer :: imdet
     integer :: nddt
     integer :: icluster
     integer :: iwrtpop
     integer :: iwrtvel
     integer :: iwrteigen
     integer :: iwrtefermi
     integer :: iwrtdos
     integer :: iwrtdosng
     integer :: ifixcharge
     integer :: iwrtewf
     integer :: iwrtatom
     integer :: iewform
     integer :: npbands
     integer :: mix_embedding
     integer :: iwrtcharges
     integer :: verbosity
     real :: cut_embedding
     real :: tempfe
     real :: sigmatol
     character(len=100) :: basis
     character(len = 20) :: dftd3_func
     integer, dimension(200) :: pbands     
  end type fb_nml_type


contains

  ! --------------------------------------
  ! Get QM energy and forces from Fireball
  ! --------------------------------------
  subroutine get_fb_forces( do_grad, nstep, ntpr_default, id, nqmatoms, qmcoords,&
       qmtypes, nclatoms, clcoords, escf, dxyzqm, dxyzcl, charge, spinmult)

    use qm2_extern_util_module, only: print_results, check_installation, write_dipole, write_charges
    use constants, only: ZERO, EV_TO_KCAL
    use ElementOrbitalIndex, only : elementSymbol
        use qmmm_module, only : qmmm_struct

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
    _REAL_              :: mulliken_charge
    _REAL_              :: total_mulliken_charge
    

    type(fb_nml_type), save :: fb_nml
    logical, save :: first_call = .true.
    integer :: i
    integer :: itime_step = 0
    integer :: printed =-1 ! Used to tell if we have printed this step yet 
                           ! since the same step may be called multiple times
    character(len=150) :: call_buffer

    ! assemble input - / output data filenames

    
    ! Setup on first program call
    if ( first_call ) then
      first_call = .false.
      write (6,*) '>>> Running QM calculation with Fireball  <<<'
      call get_namelist( fb_nml )
      call write_inpfile( fb_nml, qmcoords, nqmatoms, charge)
      call initbasics ()
      call readdata ()
      itime_step = 0
    end if


    itime_step = itime_step +1
    call fireball_qmmm_loop(itime_step,qmcoords,nclatoms,clcoords,escf,dxyzqm,dxyzcl,qmcharges)

    if ( do_grad ) then
       ! Convert eV/A -> kcal/(mol*A)
       dxyzqm(:,:) = dxyzqm(:,:)*EV_TO_KCAL
       if ( nclatoms > 0 ) then
          dxyzcl(:,:) = dxyzcl(:,:)
       end if
    else
       dxyzqm = ZERO
       if ( nclatoms > 0 ) dxyzcl = ZERO
    end if

    escf = escf*EV_TO_KCAL
    

!    do i=1,nqmatoms
!     mulliken_charge = qmcharges(i)
!     total_mulliken_charge=total_mulliken_charge+mulliken_charge
!        write(6,'(" ",i5,"      ",A2,"        ",F14.3)') i, &
!        elementSymbol(qmmm_struct%iqm_atomic_numbers(i)), &
!        mulliken_charge
!    end do
!    write(6,'(" Total Mulliken Charge =",F12.3)') total_mulliken_charge





  end subroutine get_fb_forces


  subroutine get_namelist(fb_nml)

    use UtilitiesModule, only: Upcase
    implicit none
    character(len=20) :: executable
    !integer, intent(in) :: ntpr_default
    type(fb_nml_type), intent(out) :: fb_nml
    integer ::  max_scf_iterations, qmmm_int, idftd3, debug, iqout, iensemble,        &
               imcweda, ihorsfield, imdet, nddt, icluster, iwrtpop, iwrtvel, iwrteigen,   &
               iwrtefermi, iwrtdos, iwrtdosng, ifixcharge, iwrtewf, iwrtatom, iewform, idipole,      &
               npbands, mix_embedding, iwrtcharges, verbosity
    real :: cut_embedding, tempfe, sigmatol
    character (len = 100) :: basis
    character (len = 20) :: dftd3_func
    integer, dimension(200) :: pbands
    namelist /fb/ executable, max_scf_iterations, qmmm_int, idftd3, debug, iqout,   &
               imcweda, ihorsfield, iensemble, imdet, nddt, icluster, iwrtpop, iwrtvel,  &
               iwrteigen, iwrtefermi, iwrtdos, iwrtdosng, ifixcharge, iwrtewf, iwrtatom,            &
               iewform, npbands, idipole, pbands, mix_embedding, cut_embedding,          &
               iwrtcharges, tempfe, sigmatol, basis, dftd3_func, verbosity
    integer :: i, ierr

    ! Default values
    executable = "fireball_server"
    max_scf_iterations = 70
    qmmm_int = 1
    idftd3 = 0
    idipole = 1
    debug = 0
    iqout = 1
    ifixcharge = 1
    idipole = 1
    imcweda = 1
    ihorsfield = 0
    iensemble = 0
    imdet = 0
    nddt = 1000
    icluster = 1
    iwrtpop = 0
    iwrtvel = 0
    iwrteigen = 0
    iwrtefermi = 0 
    iwrtdos = 0
    iwrtdosng = 0
    ifixcharge = 0
    iwrtewf = 0
    iwrtatom = 0
    iewform = 0
    npbands = 0
    pbands = 0
    mix_embedding = 0
    cut_embedding = 99.0d0
    iwrtcharges = 0
    tempfe = 100.0d0
    sigmatol = 1.0E-8
    basis = 'Fdata'
    dftd3_func = 'fb_bio'
    verbosity = 0   

    ! Read namelist, 
    rewind 5
    read(5,nml=fb,iostat=ierr)


    ! Assign namelist values to fb_nml data type
    fb_nml%executable           = executable
    fb_nml%max_scf_iterations   = max_scf_iterations
    fb_nml%qmmm_int             = qmmm_int
    fb_nml%idftd3               = idftd3
    fb_nml%idipole              = idipole
    fb_nml%debug                = debug
    fb_nml%iqout                = iqout
    fb_nml%imcweda              = imcweda
    fb_nml%ihorsfield           = ihorsfield
    fb_nml%iensemble            = iensemble
    fb_nml%imdet                = imdet
    fb_nml%nddt                 = nddt
    fb_nml%icluster             = icluster  
    fb_nml%iwrtpop              = iwrtpop
    fb_nml%iwrtvel              = iwrtvel
    fb_nml%iwrteigen            = iwrteigen
    fb_nml%iwrtefermi           = iwrtefermi  
    fb_nml%iwrtdos              = iwrtdos
    fb_nml%iwrtdosng            = iwrtdosng
    fb_nml%ifixcharge           = ifixcharge
    fb_nml%iwrtewf              = iwrtewf
    fb_nml%iwrtatom             = iwrtatom
    fb_nml%iewform              = iewform
    fb_nml%npbands              = npbands
    fb_nml%pbands               = pbands
    fb_nml%mix_embedding        = mix_embedding
    fb_nml%cut_embedding        = cut_embedding
    fb_nml%iwrtcharges          = iwrtcharges
    fb_nml%tempfe               = tempfe
    fb_nml%sigmatol             = sigmatol
    fb_nml%basis                = basis
    fb_nml%dftd3_func           = dftd3_func
    fb_nml%verbosity            = verbosity

    ! Need this variable so we don't call MPI_Send in the finalize subroutine

  end subroutine get_namelist

  ! --------------------------------
  ! Print Fireball namelist settings
  ! --------------------------------
  subroutine write_inpfile( fb_nml, qmcoords, nqmatoms, charge)

    use qmmm_module, only : qmmm_struct

    implicit none
    type(fb_nml_type), intent(in) :: fb_nml
    _REAL_           , intent(in) :: qmcoords(:,:)
    integer          , intent(in) :: charge
    integer          , intent(in) :: nqmatoms
    integer :: k

    open (unit = 226, file = 'fireball.in', status = 'unknown')

    write (226,600) "&option"
    
    if (fb_nml%ihorsfield .eq. 1) then
    write (226,600) "imcweda = 0"
    write (226,600) "ihorsfield = 1"
    end if
    write (226,601) "qstate = ", -1*charge
    write (226,602) "icluster = ", fb_nml%icluster
    write (226,602) "iqout = ", fb_nml%iqout
    write (226,603) "max_scf_iterations = ", fb_nml%max_scf_iterations
    write (226,602) "iensemble = ", fb_nml%iensemble
    write (226,602) "imdet = ", fb_nml%imdet
    write (226,605) "nddt = ", fb_nml%nddt
    write (226,602) "iqmmm = ", fb_nml%qmmm_int
    write (226,602) "ifixcharge =", fb_nml%ifixcharge
    write (226,600) "iquench = 0"
    write (226,602) "idftd3 = ",fb_nml%idftd3
    write (226,*) "dftd3_func = ", fb_nml%dftd3_func
    write (226,606) "tempfe =", fb_nml%tempfe
    write (226,607) "sigmatol =", fb_nml%sigmatol
    write (226,609) "fdataLocation = '", fb_nml%basis,"'"
    write (226,602) "idipole = ",fb_nml%idipole
    write (226,602) "verbosity = ",fb_nml%verbosity    
    write (226,602) "mix_embedding = ", fb_nml%mix_embedding
    write (226,606) "cut_embedding = ", fb_nml%cut_embedding
    write (226,600) "&end"
    write (226,600) "&output"
    write (226,602) "iwrtcharges = ",fb_nml%iwrtcharges 
    write (226,602) "iwrtpop = ", fb_nml%iwrtpop
    write (226,602) "iwrtvel = ", fb_nml%iwrtvel
    write (226,602) "iwrteigen = ", fb_nml%iwrteigen
    write (226,602) "iwrtefermi = ", fb_nml%iwrtefermi
    write (226,602) "iwrtdos = ", fb_nml%iwrtdos
    write (226,602) "iwrtdosng = ", fb_nml%iwrtdosng
    write (226,602) "iwrtewf = ",fb_nml%iwrtewf
    write (226,602) "iwrtatom = ",fb_nml%iwrtatom
    write (226,600) "&end"
    if (fb_nml%iwrtewf .eq. 1) then
    write (226,600) "&mesh"
    write (226,603) "iewform = ", fb_nml%iewform
    write (226,603) "npbands = ", fb_nml%npbands
    write (226,608) "pbands = ", fb_nml%pbands
    write (226,600) "&end"
    end if

    open (unit = 227, file = 'input.bas', status = 'unknown')

    write (227,*) nqmatoms
    do k = 1, nqmatoms
      write (227,700) qmmm_struct%iqm_atomic_numbers(k), qmcoords(:,k)
    end do



! Format Statements
! ===========================================================================
600     format (a)!(a,f8.5)
601     format (a,i2)
602     format (a,i1)
603     format (a,i3)
604     format (a,i1)
605     format (a,i4)
606     format (a,f8.3)
607     format (a,e13.6)
608     format (a,199(i4,','),i3)
609     format (a,a,a)
700     format (i2, 3(2x,f8.4))

  end subroutine write_inpfile
      

end module qm2_extern_fb_module
