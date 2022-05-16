subroutine start()
  use mpi_main
  use interactions
  call init_MPI (iammaster, iammpi, my_proc, nprocs)
  call cpu_time (time_begin)
  call initbasics ()
  call readdata ()
  call main_loop ()
  call cpu_time (time_end)
  write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'
end

subroutine f2py_run()
  use mpi_main
  use interactions
  call init_MPI (iammaster, iammpi, my_proc, nprocs)
  call cpu_time (time_begin)
  !call initbasics () = f2py_initbasics(fdatalocation) + fb.f2py_init()
  call main_loop ()
  call cpu_time (time_end)
  write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'
end


subroutine set_icluster(iclusteraux)
  use options
  integer,intent(in)::iclusteraux
  icluster=iclusteraux
end

subroutine set_iquot(iquotaux)
  use options
  integer,intent(in)::iquotaux
  iquot=iquotaux
end

subroutine set_iquench(iquenchaux)
  use options
  integer,intent(in)::iquenchaux
  iquench=iquenchaux
end
subroutine set_dt(idtaux)
  use options
  real,intent(in)::idtaux
  dt=idtaux
end

subroutine set_nstepf(instepfaux)
  use options
  integer,intent(in)::instepfaux
  nstepf=instepfaux
end


subroutine set_iwrtxyz(iwrtxyzaux)
  use options
  integer,intent(in)::iwrtxyzaux
  iwrtxyz=iwrtxyzaux
end


subroutine f2py_natoms(aux)
  use configuration
  implicit none
  integer, intent(inout) :: aux
  print*,'naux',aux
  natoms = aux
end


subroutine f2py_nucz(zaux,nuczaux)
  use configuration
  use interactions
  use charges
  implicit none
  integer, intent(in) :: zaux
  integer, intent(in), dimension(zaux) ::  nuczaux
  integer i,nucz,iatom,ispec
  logical zindata

  allocate (imass (natoms))
  do iatom = 1, natoms
   nucz = nuczaux(iatom)
    do ispec = 1, nspecies
      if (nucz .eq. nzx(ispec)) then
        zindata = .true.
        imass(iatom) = ispec
      end if
    end do
  end do
end

subroutine f2py_ratom(zaux,raux)
  use configuration
  use interactions
  implicit none
  integer, intent(in) :: zaux
  real, intent(in), dimension(zaux,3) ::  raux
  integer i,nucz, iatom,j,in1
  print*,'raux',raux,zaux
  allocate (ratom (3, natoms))
  allocate (symbol (natoms))
  do iatom = 1, natoms
  print*,raux(iatom,:)
  end do
  do iatom = 1, natoms
    in1 = imass(iatom)
    symbol(iatom) = symbolA(in1)
    ratom(:,iatom) = raux(iatom,:)
  end do
  do iatom = 1, natoms
    ratom(:,iatom)=ratom(:,iatom)*rescal
  end do 
  print*,'iatom,symbol,ratom,imass'
  do iatom = 1, natoms
     write (*,202) iatom, symbol(iatom), ratom(:,iatom), imass(iatom)
  end do
202     format (3x, i5, 7x, a2, 3(2x,f9.3), 7x, i2)
end

subroutine f2py_initbasics(f2py_fdataLocation)
   ! es como initbasics (), pero sin cargar posiciones de atomos
   use options
   use configuration
   use interactions
   use scf
   use integrals
   use outputs
   use kpoints
   use optimization
   use md
   use charges
   use barrier
   use transport
   use energy
   use forces
   use mpi_main
   use hartree_fock

   implicit none
   character (len = 200),intent(in) ::  f2py_fdataLocation
   integer iatom
   integer in1
   integer icount
   integer counter
   integer counter_ini
   integer isorp
   integer ideriv
   integer issh
   integer numorbPP_max
   integer numorb
   integer l
   integer imu

   real distance
   real, dimension (3) :: vector

   logical file_exists
   call welcome
   call initconstants (sigma, sigmaold, scf_achieved)
   call diagnostics (ioff2c, ioff3c, itestrange, testrange)
   call readparam ()
   !pero fireball.in no existe
   fdatalocation = f2py_fdatalocation
   verbosity=10
   !cambiamos readinfo para cargar Fdata completa sin necesidad de leer 
   !las posiciones en el bas
   !call readinfo ()
   call readinfoall()

   if (ivdw .eq. 1) call readvdw (nspecies, symbolA, ivdw)

   call make_munu (nspecies)
   call make_munuPP (nspecies)
   call make_munuS (nspecies)
   if (idipole .eq. 1) then
     call make_munuDipY (nspecies)
     call make_munuDipX (nspecies)
   end if

   ! Count the maximum number of orbital interactions between any given two atoms.
   numorb_max = 0
   do in1 = 1, nspecies
    numorb = 0
    do issh = 1, nssh(in1)
     numorb = numorb + 2*lssh(issh,in1) + 1
    end do
    if (numorb .gt. numorb_max) numorb_max = numorb
   end do
   isorpmax = 0
   if (itheory .eq. 1) then
    do in1 = 1, nspecies
     isorpmax = max(isorpmax,nssh(in1))
    end do
   end if
   isorpmax_xc = 0
   do in1 = 1, nspecies
      isorpmax_xc = max(isorpmax_xc,nssh(in1))
   end do

   ideriv_max = 0
   if (itheory .eq. 1) ideriv_max = 6

   icount = 0
   ind2c = 0
   icount = icount + 1
   ind2c(1,0) = icount
   do isorp = 0, isorpmax
    icount = icount + 1
    ind2c(2,isorp) = icount
   end do
   do isorp = 0, isorpmax
    icount = icount + 1
   ind2c(3,isorp) = icount
   end do
   do isorp = 0, isorpmax
    icount = icount + 1
    ind2c(4,isorp) = icount
   end do
   icount = icount + 1
   ind2c(5,0) = icount
   do ideriv = 0, 4
    icount = icount + 1
    ind2c(6,ideriv) = icount
   end do
   do ideriv = 0, 4
    icount = icount + 1
    ind2c(7,ideriv) = icount
   end do
   do ideriv = 0, 4
    icount = icount + 1
    ind2c(8,ideriv) = icount
   end do
   icount = icount + 1
   ind2c(9,0) = icount
   icount = icount + 1
   ind2c(10,0) = icount
   icount = icount + 1
   ind2c(11,0) = icount
   icount = icount + 1
   ind2c(12,0) = icount
   icount = icount + 1
   ind2c(13,0) = icount
   icount = icount + 1
   ind2c(14,0) = icount
   if (itheory_xc .eq. 1 .or. itheory_xc .eq. 2 .or. itheory_xc .eq. 4 ) then
     if (itheory_xc .eq. 4) then
       icount = icount + 1
       ind2c(14,0) = icount
      end if !end if itheory_xc .eq. 4
      do isorp = 1, isorpmax_xc
       icount = icount + 1
       ind2c(15,isorp) = icount
      end do
      do isorp = 1, isorpmax_xc
       icount = icount + 1
       ind2c(16,isorp) = icount
      end do
      do isorp = 1, isorpmax_xc
       icount = icount + 1
       ind2c(17,isorp) = icount
      end do
      do isorp = 1, isorpmax_xc
       icount = icount + 1
       ind2c(18,isorp) = icount
      end do
      do isorp = 1, isorpmax_xc
       icount = icount + 1
       ind2c(19,isorp) = icount
      end do
      do isorp = 1, isorpmax_xc
       icount = icount + 1
       ind2c(20,isorp) = icount
      end do
      do isorp = 1, isorpmax_xc
       icount = icount + 1
       ind2c(21,isorp) = icount
      end do
      do isorp = 1, isorpmax_xc
       icount = icount + 1
       ind2c(22,isorp) = icount
      end do
      icount = icount + 1
      ind2c(23,0) = icount
   end if
   interactions2c_max = icount
   !orbital to shell
   if (itheory_xc .eq. 4) then
     allocate (orb2shell (numorb_max,nspecies) )
     do in1 = 1,nspecies
      counter = 1
      do issh = 1,nssh(in1)
         counter_ini = counter
         l = lssh(issh,in1)
         do imu = counter_ini,counter_ini+2*l
            orb2shell(imu,in1) = issh
            counter = imu+1
         end do !end imu
      end do ! end do issh = 1,nssh(in1)   
     end do ! end do in1 = 1,nspecies
     counter = imu+1
   end if ! end if itheory = 4  

  ! call get_info_orbital (natoms)

   numorbPP_max = 0
   do in1 = 1, nspecies
    numorb = 0
    do issh = 1, nsshPP(in1)
     numorb = numorb + 2*lsshPP(issh,in1) + 1
    end do
    if (numorb .gt.  numorbPP_max) numorbPP_max = numorb
   end do
   if (numorbPP_max .gt.  numorb_max) numorb_max = numorbPP_max
   if (itheory .eq. 2) call make_mu2shell (nspecies)

   call initamat(nspecies)
   call readdata ()
   print*,'load f2py_fdatalocation = ',f2py_fdatalocation

end subroutine


subroutine f2py_init() !zauxf2py) !,pos,Zin)
   use options
   use configuration
   use interactions
   use scf
   use integrals
   use outputs
   use kpoints
   use optimization
   use md
   use charges
   use barrier
   use transport
   use energy
   use neighbor_map
   use forces
   use mpi_main
   use hartree_fock

   implicit none
   integer iatom
   integer in1
   integer icount
   integer counter
   integer counter_ini
   integer isorp
   integer ideriv
   integer issh
   integer numorbPP_max
   integer numorb
   integer l
   integer imu

   real distance
   real, dimension (3) :: vector

   logical file_exists
! Initialize aux. variable
        if (nstepi .eq. 1) then
         T_average = T_initial
         T_previous = 0.0d0
         time = 0.0d0
        end if

! Allocate more arrays.
        allocate (degelec (natoms))
!        allocate (imass (natoms))
!        allocate (ratom (3, natoms))
        allocate (nowMinusInitialPos (3, natoms))
        allocate (initialPosition (3, natoms))
        allocate (vatom (3, natoms))
!        allocate (symbol (natoms))
        allocate (xmass (natoms))
        allocate (ximage (3, natoms))
        allocate (mask (3,natoms))
        mask = 1.0d0
        ximage = 0.0d0
! decide if need to shift atoms since one is situated at origin
        ishiftO = 0
        do iatom = 1, natoms
         distance = ratom(1,iatom)**2 + ratom(2,iatom)**2 + ratom(3,iatom)**2
         distance = sqrt(distance)
         if (distance .lt. 1.0d-4) ishiftO = 1
        end do

! Read data from the lattice vectors file - XXX.lvs.
! Set up the boxes surrounding the central unit cell.
        call readlvs (lvsfile, a1vec, a2vec, a3vec, icluster, rescal)

! Define volume of unit cell
        call cross (a2vec, a3vec, vector)
        Vouc=a1vec(1)*vector(1)+a1vec(2)*vector(2)+a1vec(3)*vector(3)
        call initboxes (1)

! Get kpoints for Brillouin zone integration or bandstructure calculation
		if(iclassicMD /= 1) then !not usefull with empirical potentials
	        call getkpoints(icluster, Vouc, a1vec, a2vec, a3vec, &
    	 &       lvsfile, basisfile, iquench, ireducekpts, rescal)
		end if

        if (iordern .eq. 1 .and. nkpoints .ne. 1) then
         write (*,*) ' Order-N method only works with one k-point! '
         write (*,*) ' Sorry, you must abort your run! '
         stop
        end if

! Read in a FRAGMENTS file (if it exists).  Use it to fix the internal
! geometry of various parts of the system.
        call readfragments ()

! Initialize the masses from the info.dat file.
        call initmasses (natoms, symbol, smass, xmass)

! Read information from the quench.optional file. This file contains all the
! information needed for the quenching cycles.
        call readquench (iquench, dt, energy_tol, force_tol, iensemble,       &
     &                   T_initial, T_want, taurelax)

! Read information from cgmin.input file. This file contains informations
! needed for running conjugate gradient minimization.
        if(iquench == -4 .or. iquench == -5 ) then
           call readcgo ( natoms, iforce )
        endif

! read optionally information for transport calculation
       if (itrans .eq. 1) then
           call readtrans ( natoms )
       endif


! This little section will incrementally increase/decrease the temperature
! over the course of a calculation.
         if (iendtemp .eq. 1 .and. iquench .eq. 0) then

              if (T_initial .eq. T_final) then
                write (*,*) ''
                write (*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                write (*,*) 'T_initial = T_final'
                write (*,*) 'If you wanted this then set iendtemp = 0'
                write (*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                write (*,*) ''
                stop
              end if

           T_increment = (T_final - T_initial)/nstepf
         else if (iendtemp .eq. 0) then
           T_increment = 0.00
         end if

         if (iendtemp .eq. 1 .and. iquench .ne. 0) then
           write (*,*) ''
           write (*,*) ''
           write (*,*) 'STOPPING'
           write (*,*) 'iendtemp = 1  and iquench != 0'
           write (*,*) ''
           write (*,*) ''
           stop
         end if
! ***************************************************************************

! Read information from the vdw.optional file. This file contains all of
! the information needed for adding in the vdw interations.
        if (ivdw .eq. 1) call readvdw (nspecies, symbolA, ivdw)

! Read information from the barrier.optional file.  This file contains all
! of the information needed for calculating the crude energy barrier.
        allocate (ratom_final (3, natoms))
        call readbarrier (natoms, nspecies, imass, nzx)

! Read the dos.input file if iwrtdos is greater than 1 CGP
        if (iwrtdos.ge.1 .or. iwrtatom .ge. 1) then
          call readdos ( )
        end if

! Dynamical matrix calculation
        if (idynmat .eq. 1) then
! read input file
           call readphi (natoms, ratom, nstepi, nstepf)
! perform first displacement
           call bvec ( nstepi, natoms, nstepf, ratom)
! jel-eph
! Electron-phono coupling
       	   if (iephc .eq. 1) call readephc (natoms)
        endif

! write out population analysis of MOs
        if (iwrtpop .eq. 1) then 

          write (*,*) '   ========   Initialize population analysis tool '
          inquire(file="pop.optional", exist=file_exists )
! set energy range from input file 'pop.optional'
          if ( file_exists ) then
           open (unit= 33, file="pop.optional", status='unknown')
           write (*,*) '    File pop.optional exists!'
           write (*,*) '    Reading parameters from pop.optional file'
           read (33,*)  Epop_L
           read (33,*)  Epop_U
           write (*,*) '    Epop_L =',Epop_L
           write (*,*) '    Epop_U =',Epop_U
          else
           Epop_L = - 6.0d0
           Epop_U = 0.0d0 
           write (*,*) '    File pop.optional does not exist!'
           write (*,*) '    Setting default parameters'
           write (*,*) '    Epop_L =',Epop_L
           write (*,*) '    Epop_U =',Epop_U
          endif ! file_exists
           write (*,*) '    The results will be written into populations.dat file '
        endif ! iwrtpop


! Count the orbitals
        norbitals = 0
        do iatom = 1, natoms
         in1 = imass(iatom)
         norbitals = norbitals + num_orb(in1)
        end do

! Count the total number of shells in the system.
        nssh_tot = 0
        do iatom = 1, natoms
         in1 = imass(iatom)
         do issh = 1, nssh(in1)
          nssh_tot = nssh_tot + 1
         end do
        end do


! We originally use the neutral atom charges from the info.dat file.
! If this is a restart run, or the user wishes to start with charges other
! than the neutral atom stuff, then information is read in from a charge file.
        if (iKS .eq. 1) then
         call initcharges_KS (natoms, nspecies, itheory, ifixcharge, symbol)
        else
         call initcharges (natoms, nspecies, itheory, ifixcharge, symbol)
        endif

! Calculate ztot
        ztot = 0.0d0
        nelectron = 0.0d0
        do iatom = 1, natoms
         in1 = imass(iatom)
         do issh = 1, nssh(in1)
          ztot = ztot + Qneutral(issh,in1)
          nelectron(iatom) = nelectron(iatom) + Qneutral(issh,in1)
         end do
        end do

! ADD Feb. 15, 2006, need to be tested.
        if (abs(qstate) .gt. 0.00001) ztot = ztot + qstate

! Calculate degelec.  We only need this once at the beginning of the simulation.
        degelec(1) = 0
        do iatom = 2, natoms
         degelec(iatom) = 0
         in1 = imass(iatom - 1)
         degelec(iatom) = degelec(iatom - 1) + num_orb(in1)
        end do

! NPA initialize auxillary arrays (dani goes to hollywood)
         call get_info_orbital (natoms)

! This is a fix from Jose - evidently, there were problems with transition
! metals where the orbitals are different than the pseudopotential.
! For the dimensions of pp-arrays, we need to fix numorb_max to something
! different. For the moment find numorbPP_max and choose the greater of
! numorb_max and numorbPP_max. This can be improved later.
        numorbPP_max = 0
        do in1 = 1, nspecies
         numorb = 0
         do issh = 1, nsshPP(in1)
          numorb = numorb + 2*lsshPP(issh,in1) + 1
         end do
         if (numorb .gt.  numorbPP_max) numorbPP_max = numorb
        end do
        if (numorbPP_max .gt.  numorb_max) numorb_max = numorbPP_max

! Call make_mu2shell if we are doing extended hubbard.  The variable numorb_max
! is needed here, so call after finding numorb_max.
        if (itheory .eq. 2) call make_mu2shell (nspecies)

! Initialize the amat array for twister routines and set haveDorbitals
        call initamat(nspecies)


! check if we need the grid
       igrid = 0
       if (iwrtden .eq. 1) igrid = 1
       if (iwrtewf .eq. 1) igrid = 1
       if (iks .eq. 1) igrid = 1
       if (iwrtdipole .eq. 1) igrid = 1
       if (igrid .eq. 1) then
! call readgrid before initconstraints subroutine to avoid atom shift
! when we fix the mesh position
         call readgrid (iwrtewf)
       endif

! read bais option
        if (ibias .eq. 1) then
         write (*,*) ' Read bias parameters'
         call readbias (natoms)
         write (*,*) ' Allocate bias arrays'
         call allocate_bias (natoms)
        endif

! Apply constraints to initial velocities.
       call initconstraints (iconstraints, iensemble,   &
     &                        T_initial, ibarrier, ratom_final, imass,   &
     &                        fixCenOfMass, rcmOld, xmassTot)
       initialPosition = ratom !do this POST shift

       if (igrid .eq. 1) then
         call allocate_grid (natoms, nspecies)
         write (*,*) ' call read_wf'
         call read_wf ()
         write (*,*) ' call read_vna'
         call read_vna ()
         write (*,*) ' initialize grid'
         call initgrid (icluster)
        endif

! NEB
! initialize parameters of NEB methof
        call initneb (natoms, nspecies, imass, nzx, ratom)
! end NEB

! Initialize xdot
        allocate (xdot (0:5, 3, natoms)) ! Initialized below

        if (nstepi .ne. 1) then
         call restart (nstepi, dt, acfile, xvfile, T_average, &
     &                 T_previous, time)
         if (ishiftO .eq. 1) then
          do iatom = 1, natoms
           ratom(:,iatom) = ratom(:,iatom) + shifter
           xdot(0,:,iatom) = xdot(0,:,iatom) + shifter
          end do
         end if
        else
         xdot = 0.0d0
         xdot(0,:,1:natoms) = ratom(:,1:natoms)
         xdot(1,:,1:natoms) = vatom(:,1:natoms)
        end if

! Calculate the atomic energy.
        call initatomicE (natoms, etotatom, imass, atomic_energy)

! Set the gear algorithm constants.
        call setgear

! Initialize the thermostat info
        if (iensemble .eq. 2) call initNH(natoms,T_want)

! Allocate the stuff that depends on natoms, neigh_max, and numorb_max
        write (*,*) ' Initiallizing arrays '
        call allocate_neigh (nprocs, my_proc, iordern, icluster,     &
     &                       ivdw, ifixneigh, iwrthampiece,  &
     &                       iwrtatom)
        call allocate_f (natoms, neigh_max, neighPP_max, numorb_max, nsh_max,&
     &                   itheory, itheory_xc, igauss, ivdw, iharmonic, ibias)
        call allocate_h (natoms, neigh_max, neighPP_max, itheory, itheory_xc,&
     &                   igauss, iwrtdos, iwrthop, iwrtatom)
! jel-grid
        call allocate_rho (natoms, neigh_max, neighPP_max, numorb_max,       &
     &                     nsh_max, itheory_xc, igrid)
! end jel-grid
        call allocate_dos (natoms, iwrtdos, iwrthop)
! itrans
        if (itrans .eq. 1) then
         call readbind ()
         call allocate_trans (nkpoints)
         if (ifithop .eq. 1) call readhop ( nspecies )
        endif

        if (nstepi .eq. 1) then
         etotnew = 0.0d0
         if (.not. allocated (ftotnew)) allocate (ftotnew (3, natoms))
         ftotnew = 0.0d0
        end if

! jel-grid
! Initialize the density matrix
        if (igrid .eq. 1 ) then
         call neighbors (nprocs, my_proc, iordern, icluster,      &
     &                      iwrtneigh, ivdw)
         !SFIRE  APRIL 2018
         call neighbors_pairs(icluster)
         !SFIRE  APRIL 2018
         write (*,*) 'Initialize density matrix'
         call initdenmat (natoms)
! end jel-grid
!CHROM
        elseif ( iclassicMD > 0 .and. igrid /= 1 )then
		 	call neighbors (nprocs, my_proc, iordern, icluster, iwrtneigh, ivdw)
                         !SFIRE  APRIL 2018
                        call neighbors_pairs(icluster)
                         !SFIRE  APRIL 2018
		endif
!END CHROM

! initialize time dependent variables
        if (itdse .eq. 1) then
! allocate arrays
		 write (*,*) ' Read TD parameters'
         call readtdse ()
         write (*,*) ' Allocate TD-matrices'
         call allocate_tdse ()
        endif

! initialize cDFT
    !    if (icDFT .eq. 1) then
    !     call initcDFT ()
    !    endif

! GAP ENRIQUE-FF
        if (igap .eq. 1) then
	  call readhartree (nspecies,natoms)
            max_scf_iterations = 2
            sigmatol = 1e-20
        else if (igap .eq. 2) then
          call readhartree (nspecies,natoms)
            max_scf_iterations = 2
            sigmatol = 1e-20
            iwrtcdcoefs = 1
        else if (igap .eq. 3) then
            allocate(hs_mat(norbitals,norbitals))
        end if
! end GAP ENRIQUE-FF


! End Procedure
! ============================================================================
   return


end


subroutine readinfoall ()
  use charges
  use dimensions
  use interactions
  use configuration 
  use integrals
  use options, only : verbosity, inputxyz
  implicit none
  integer ispec
  integer jspec
  integer iatom
  integer iline
  integer issh
  integer imu
  integer ins
  integer nzx_temp
  real, dimension (:,:), allocatable :: cutoff  ! cutoff radius in bohr
  real :: nucz
  character (len=25) :: signature
  character (len=2)  :: symbolA_temp
  character (len=1)  :: charx
  integer nchar
  integer foundspecies
   foundspecies=0
  write (*,*) ' Now we are reading from the info.dat file. '
  write (*,*) '  '
  open (unit = 12, file = trim(fdataLocation)//'/info.dat', status = 'old')
  read (12,101) signature
  read (12,*) nspecies
  allocate (nzx (nspecies))
  allocate (symbolA (nspecies)) 
  allocate (etotatom (nspecies)) 
  allocate (smass (nspecies)) 
  allocate (rc_PP (nspecies)) 
  allocate (rcutoff (nspecies, nsh_max)) 
  allocate (cl_PP (0:nsh_max - 1, nspec_max))
  allocate (nssh (nspec_max))
  allocate (lssh (nsh_max, nspec_max))
  allocate (nsshPP (nspec_max))
  allocate (lsshPP (nsh_max, nspec_max))
  allocate (Qneutral (nsh_max, nspec_max))
  allocate (cutoff (nsh_max, nspec_max))
  allocate (wavefxn (nsh_max, nspec_max))
  allocate (napot (0:nsh_max, nspec_max))
  ispec = 1
  nsup = 0
  do jspec = 1, nspecies
   read (12,*)
   read (12,*)
   read (12,102) symbolA_temp
   read (12,*) nzx_temp
   if (verbosity .ge. 3) write(*,102) symbolA_temp
   if (verbosity .ge. 3) write (*,*) nzx_temp
   foundspecies=foundspecies+1
   symbolA(ispec) = symbolA_temp
   nzx(ispec) = nzx_temp
   read (12,*) smass(ispec)
   read (12,*) nssh(ispec)
   if (nssh(ispec) .gt. nsh_max) then
    write (*,*) ' nssh(ispec) = ', nssh(ispec),' nsh_max = ', nsh_max
    write (*,*) ' Sorry -- redimension nsh_max in MODULES/dimensions.f90'
    stop
   end if
   if (nssh(ispec) .gt. 8) then
    write (*,*) ' nssh(ispec) = ', nssh(ispec)
    write (*,*) ' Sorry -- Currently the basis set cannot be larger '
    write (*,*) ' than s, s*, p, p*, d, d*, f, and f* '
    stop
   end if
   read (12,*) (lssh(issh,ispec), issh = 1, nssh(ispec))
   read (12,*) nsshPP(ispec)
   read (12,*) (lsshPP(issh,ispec), issh = 1, nsshPP(ispec))
   read (12,*) rc_PP(ispec)
   read (12,*) (Qneutral(issh,ispec), issh = 1, nssh(ispec))
   read (12,*) (cutoff(issh,ispec), issh = 1, nssh(ispec))
   do issh = 1, nssh(ispec)
    rcutoff(ispec, issh) = cutoff(issh,ispec)*0.529177d0
   end do
    read (12,103) (wavefxn(issh,ispec), issh = 1, nssh(ispec))
   read (12,103) (napot(issh,ispec), issh = 0, nssh(ispec))
! jel-grid
! adjust the potential file name
   do issh = 0, nssh(ispec)
    nchar = len_trim(napot(issh,ispec))
    do imu = 1,nchar
     charx = napot(issh,ispec)(imu:imu)
     if ( lle(charx,'/') .and. lge(charx,'/') ) then 
      ins = imu
     endif
    enddo ! do imu
    napot(issh,ispec) = trim('/basis/'//napot(issh,ispec)(ins+1:nchar))
   enddo ! do issh
! adjust the wavefunction file name
   do issh = 1, nssh(ispec)
    nchar = len_trim(wavefxn(issh,ispec))
    do imu = 1,nchar
     charx = wavefxn(issh,ispec)(imu:imu)
     if ( lle(charx,'/') .and. lge(charx,'/') ) then 
      ins = imu
     endif
    enddo ! do imu
   wavefxn(issh,ispec) = trim('/basis/'//wavefxn(issh,ispec)(ins+1:nchar))
   enddo ! do issh 
! end je-grid
   read (12,*) etotatom(ispec)
   read (12,*)

! Write out.
   if (verbosity .ge. 3)  then
     write (*,100)
     write (*,301) ispec
     write (*,302) symbolA(ispec)
     write (*,303) nzx(ispec)
     write (*,304) smass(ispec)
     write (*,305) nssh(ispec)
     write (*,306) (lssh(issh,ispec), issh = 1, nssh(ispec))
     write (*,307) nsshPP(ispec)
     write (*,308) (lsshPP(issh,ispec), issh = 1, nsshPP(ispec))
     write (*,314) rc_PP(ispec)
     write (*,309) (Qneutral(issh,ispec), issh = 1, nssh(ispec))
     write (*,310) (cutoff(issh,ispec), issh = 1, nssh(ispec))
     write (*,311) (wavefxn(issh,ispec), issh = 1, nssh(ispec))
     write (*,312) (napot(issh,ispec), issh = 0, nssh(ispec))
     write (*,313) etotatom(ispec)
     write (*,100)
   endif !verbosity 
! Increment ispec, since we read in data
   ispec = ispec + 1
 end do
  

 if(foundspecies /= nspecies) then
   write(*,*)'In Fdata/info.dat is defined only ',foundspecies,' species from ',nspecies
   write(*,*)'Exiting. Please check the info.dat'
   stop
 end if

  deallocate (cutoff)
! end jel-grid
 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
101     format (2x, a25)
102     format (2x, a2)
103     format (9(2x,a25))
301     format (2x, i2, ' - Information for this species ')
302     format (2x, a2, ' - Element ')
303     format (2x, i3, ' - Nuclear Z ')
304     format (2x, f7.3, ' - Atomic Mass ')
305     format (2x, i2, ' - Number of shells ')
306     format (2x, 8(2x,i1), ' - L; quantum number for each shell ')
307     format (2x, i2, ' - Number of shells (Pseudopotential) ')
308     format (2x, 8(2x,i1), ' - L; quantum number for each shell ')
309     format (2x, 8(2x,f5.2), ' - Occupation numbers ')
310     format (2x, 8(2x,f5.2), ' - Radial cutoffs ')
311     format (2x, 9(2x,a25), ' - Wavefunction files ')
312     format (2x, 9(2x,a25), ' - (Non)-neutral atom potentials ')
313     format (2x, f12.4, ' - Atomic energy ')
314     format (2x, f12.4, ' - Radial cutoffs (Pseudopotential) ')
        close (unit = 12)

end



