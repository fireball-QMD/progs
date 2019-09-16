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

! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in
! subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! mbs_corrections.f90
! Program Description
! ===========================================================================
!     Add MBD many body dispersion corrections with libmbd
!     (https://github.com/jhrmnn/libmbd) 
!
! ===========================================================================
! Code written by:
! Jesús I. Mendieta
! Departamento de físca teórica de la materia condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================

        subroutine mbd_corrections

        use mbd, only: mbd_input_t, mbd_calc_t 
        use configuration
        use interactions
        use neighbor_map
        use energy
        use forces
        use charges
        use density
        use options
       

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        type(mbd_input_t) :: input
        type(mbd_calc_t) :: mbDispersion
        integer :: iatom,jatom,i
        real :: eMbd
        real :: grads(3, natoms)
        real :: stress(3,3)
        real :: lvs33(3,3)
        real :: cpa(natoms) 
        integer imu, inu
        integer in1, in2
        integer issh
        integer mqn
        real, dimension (numorb_max, natoms) :: QMulliken_in
        real, dimension (natoms) :: QMulliken_in_TOT
        real, dimension (natoms) :: Qneutral_TOT

        character(len=30) :: mbd_theory = 'ts'
!        logical :: first_time = .true.
        character(len=30) :: method = 'ts'!'mbd-rsscs'
        real :: ts_ene_acc = 1d-6
        real :: ts_f_acc = 1d-7
        real :: ts_d = 20d0
        real :: ts_sr = -1
        integer :: n_omega_grid = 15
        integer :: k_grid(3) = [1, 1, 1]
        real :: k_grid_shift = 0.5d0

integer :: code
character(200) :: origin, msg


! Procedure
! ===========================================================================
        etot_mbd = 0.0d0
        ftot_mbd = 0.0d0

!         if (first_time) then
!           if (mbd_theory .eq. 'ts') then
             !input%method = 'ts'
             input%ts_ene_acc = ts_ene_acc
             input%ts_f_acc = ts_f_acc
             input%ts_d = ts_d
             input%ts_sr = ts_sr   
             input%vdw_params_kind = 'ts'
!           elseif (mbd_theory .eq. 'mbd') then
             input%method = 'mbd-rsscs'
!             input%n_omega_grid = n_omega_grid
             input%k_grid = k_grid
             input%k_grid_shift = k_grid_shift
             input%xc = 'pbe'
!             input%vdw_params_kind = 'ts'
!           else
!             write(*,*) 'Invalid dispersion model name'
!             stop
!           endif

!do iatom = 1, natoms
!write(*,*) ratom(:,iatom)
!enddo
           input%atom_types = symbol
           input%coords = ratom / 0.52917726
           input%calculate_forces = .true.

!           inp%calculate_forces = tForces
!            inp%atom_types = speciesName(species0)
!            inp%coords = coord0
!            if (tPeriodic) inp%lattice_vectors = latVec
!            call mbDispersion%init(inp)

           call mbDispersion%init(input)
call mbDispersion%get_exception(code, origin, msg)
if (code > 0) then
    print *, msg
    stop
end if

!         else
!         endif

        if (icluster .eq. 0) then
          do i = 1,3
            lvs33(1,i)=a1vec(i)
            lvs33(2,i)=a2vec(i)
            lvs33(3,i)=a3vec(i)
          enddo
        endif
           QMulliken_in = 0.0d0
           QMulliken_in_TOT = 0.0d0
           Qneutral_TOT = 0.0d0
           eMBd = 0
           grads = 0.0d0


          do iatom = 1, natoms
           in1 = imass(iatom)
           jatom = neigh_self(iatom)
           in2 = imass(iatom)

           do imu = 1, num_orb(in1)
            do inu = 1, num_orb(in2)
             QMulliken_in_TOT(iatom) = QMulliken_in_TOT(iatom)                    &
     &       + rho(imu,inu,jatom,iatom)*s_mat(imu,inu,jatom,iatom)
            end do
           end do
! Finally the imu loop.

           do issh = 1, nssh(in1)
             Qneutral_TOT(iatom) = Qneutral_TOT(iatom)  + Qneutral(issh,in1)
           end do
           !write(*,*) 'iatom:', iatom, 'Qneutral:',Qneutral_TOT(iatom)
           cpa(iatom) =  QMulliken_in_TOT(iatom) / Qneutral_TOT(iatom)
           !write(*,*) 'iatom:', iatom, 'cpa:', cpa(iatom)
          end do
          ! End loop over atoms

!write(*,*)'size(ratom)', size(ratom)
!write(*,*)'size(grads)', size(grads)

        !from A to bohr ratom/0.52917726
        if (icluster .eq. 0) call mbDispersion%update_lattice_vectors(lvs33 / 0.52917726)
        call mbDispersion%update_coords(ratom / 0.52917726)
        call mbDispersion%update_vdw_params_from_ratios(cpa)
        call mbDispersion%evaluate_vdw_method(eMbd)
!write(*,*) 'eMbd', eMbd
        call mbDispersion%get_gradients(grads)
!write(*,*) 'calculated mbd gradients', grads
        if (icluster .eq. 0) call mbDispersion%get_lattice_derivs(stress)
        call mbDispersion%destroy

        !energy from a.u. to eV * 27.2107
        etot_mbd = eMbd * 27.2107
write(*,*) 'etot_mbd', etot_mbd
        !forces from a.u. (Eh/ao) to eV/A * 51.421
        do iatom = 1, natoms
          ftot_mbd(:,iatom) = - grads(:,iatom) * 51.421
        enddo

        endsubroutine mbd_corrections
