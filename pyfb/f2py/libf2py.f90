subroutine start()
  use mpi_main
  use interactions
  wrtout = .true.
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
  wrtout = .true.
  call cpu_time (time_begin)
  !call initbasics () = f2py_initbasics(fdatalocation) + fb.f2py_init()
  !en f2py_initbasics cambiamos readinfo por readinfoall para que lea todo el Fdata
  !El problema es que cuando intentamos cargar otra estructura 
  !ya tiene alocateado muchas matrices con natoms y hay que reallocatear...
  !rehacerlo escribiendo en disco duro solo con start() !!!!
  !pensar si esto sirve para algo.... yo creo que no, quedarme solo con start..ufff
  !ultimo commit no push
  call main_loop ()
  call cpu_time (time_end)
  write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'
end

real function f2py_getenergy()
  use energy
  call scf_loop (1)
  call postscf ()
  call getenergy (1)
  f2py_getenergy=etot
  return
end


subroutine f2py_print_charges()
  use dimensions
  use interactions
  use charges
  use configuration
   do iatom = 1, natoms
    in1 = imass(iatom)
    write (*,601) (Qin(issh,iatom), issh = 1, nssh(in1))
   end do
601     format (2x, 10f14.8)
end

subroutine f2py_print_pcharges()
  use dimensions
  use interactions
  use charges
  use configuration
  call load_partial_charges()
  write(*,*) 'natoms = ',natoms
  do iatom = 1, natoms
    write(*,'(2x, i4,2x,f10.6)') iatom,Q_partial(iatom)
  end do
end

character(140) function f2py_pcharge(iaux)
  use dimensions
  use interactions
  use charges
  use configuration
  integer,intent(in)::iaux 
  call load_partial_charges()
  write (f2py_pcharge,'(f14.8)') Q_partial(iaux)
  return
end


character(140) function f2py_charge(iaux)
  use dimensions
  use interactions
  use charges
  use configuration
  integer,intent(in)::iaux 
  in1 = imass(iaux)
  write (f2py_charge,'(2x, 10f14.8)') (Qin(issh,iaux), issh = 1, nssh(in1))
  return
end

subroutine set_icluster(iaux)
  use options
  integer,intent(in)::iaux
  icluster=iaux
end

subroutine set_iqout(iaux)
  use options
  integer,intent(in)::iaux
  iqout=iaux
end

subroutine set_iquench(iaux)
  use options
  integer,intent(in)::iaux
  iquench=iaux
end
subroutine set_dt(aux)
  use options
  real,intent(in)::aux
  dt=aux
end

subroutine set_nstepf(aux)
  use md
  integer,intent(in)::aux
  nstepf=aux
end


subroutine set_idipole(iaux)
  use options
  integer,intent(in)::iaux
  idipole=iaux
end

subroutine set_iwrtcharges(iaux)
  use outputs
  integer,intent(in)::iaux
  iwrtcharges=iaux
end

subroutine set_iwrtdipole(iaux)
  use outputs
  integer,intent(in)::iaux
  iwrtdipole=iaux
end


subroutine set_iks(iaux)
  use options
  integer,intent(in)::iaux
  iks=iaux
end

subroutine set_imcweda(iaux)
  use options
  integer,intent(in)::iaux
  imcweda=iaux
end

subroutine set_idogs(iaux)
  use options
  integer,intent(in)::iaux
  imcweda=iaux
end



subroutine set_iwrtxyz(iaux)
  use outputs
  integer,intent(in)::iaux
  iwrtxyz=iaux
end

subroutine set_verbosity(iaux)
  use options
  integer,intent(in)::iaux
  verbosity=iaux
  end

subroutine f2py_natoms(aux)
  use configuration
  implicit none
  integer, intent(inout) :: aux
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

subroutine f2py_getbas(zaux,nuczaux,raux)
  use configuration
  use interactions
  use charges
  implicit none
  integer, intent(in) :: zaux
!  integer, intent(in) :: naux
  integer, intent(in), dimension(zaux) ::  nuczaux
  real, intent(in), dimension(zaux,3) ::  raux
  integer i,nucz, iatom,j,in1,ispec
  logical zindata
  natoms = zaux
  allocate (ratom (3, natoms))
  allocate (symbol (natoms))
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
  character (len = 200),intent(in) ::  f2py_fdataLocation
  call f2py_initbasics_opt(f2py_fdataLocation,0)
end



subroutine f2py_initbasics_opt(f2py_fdataLocation,idipoleaux)
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
   integer,intent(in)::idipoleaux
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
   !pero fireball.in no existe, ponemos opt
   idipole=idipoleaux

   fdatalocation = f2py_fdatalocation
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
   call init_MPI (iammaster, iammpi, my_proc, nprocs)
! Initialize aux. variable
        if (nstepi .eq. 1) then
         T_average = T_initial
         T_previous = 0.0d0
         time = 0.0d0
        end if
! Allocate more arrays.
!        allocate (imass (natoms))
!        allocate (ratom (3, natoms))
!        allocate (symbol (natoms))
        allocate (degelec (natoms))
        allocate (nowMinusInitialPos (3, natoms))
        allocate (initialPosition (3, natoms))
        allocate (vatom (3, natoms))
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
        write (*,*) ' Initiallizing arrays allocate_neigh '
        call allocate_neigh (nprocs, my_proc, iordern, icluster,     &
     &                       ivdw, ifixneigh, iwrthampiece,  &
     &                       iwrtatom)
        write (*,*) ' Initiallizing arrays allocate_f'
        call allocate_f (natoms, neigh_max, neighPP_max, numorb_max, nsh_max,&
     &                   itheory, itheory_xc, igauss, ivdw, iharmonic, ibias)
        if ( allocated ( dip)) deallocate (dip)
        write (*,*) ' Initiallizing arrays allocate_h'
        call allocate_h (natoms, neigh_max, neighPP_max, itheory, itheory_xc,&
     &                   igauss, iwrtdos, iwrthop, iwrtatom)
! jel-grid
        write (*,*) ' Initiallizing arrays allocate_rho'
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



subroutine f2py_deallocate_all()
 use barrier
 use bias
 use charges
 use classicMD
 use configuration
 use constants_fireball
 use density
 use dftd3_api
 use dftd3_common
 use dftd3_core
 use dftd3_pars
 use dftd3_sizes
 use dimensions
 use dynamo
 use energy
 use forces
 use fragments
 use gaussG
 use grid
 use hartree_fock
 use integrals
 use interactions
 use kpoints
 use matmult
 use md
 use MD
 use module_dos
 use mpi_main
 use neb
 use neighbor_map
 use nonadiabatic
 use noseHoover
 use omp_lib
 use optimization
 use options
 use outputs
 use scf
 use steered
 use tdse
 use transport
 use umbrella
 use vnneutral
 use wavefunction
!if ( allocated ( alfa )) deallocate ( alfa )
if ( allocated ( arhoij_off )) deallocate ( arhoij_off )
if ( allocated ( arhoi_on )) deallocate ( arhoi_on )
if ( allocated ( arho_off )) deallocate ( arho_off )
if ( allocated ( arho_on )) deallocate ( arho_on )
if ( allocated ( arhopij_off )) deallocate ( arhopij_off )
if ( allocated ( arhop_off )) deallocate ( arhop_off )
if ( allocated ( arhop_on )) deallocate ( arhop_on )
!if ( allocated ( atompos )) deallocate ( atompos )
if ( allocated ( bar_density_2c )) deallocate ( bar_density_2c )
if ( allocated ( bar_density_3c )) deallocate ( bar_density_3c )
if ( allocated ( cape )) deallocate ( cape )
if ( allocated ( cape_es )) deallocate ( cape_es )
if ( allocated ( degelec )) deallocate ( degelec )
if ( allocated ( deigen )) deallocate ( deigen )
if ( allocated ( density_2c )) deallocate ( density_2c )
if ( allocated ( density_3c )) deallocate ( density_3c )
if ( allocated ( dewald )) deallocate ( dewald )
!if ( allocated ( dewaldl )) deallocate ( dewaldl )
!if ( allocated ( dip )) deallocate ( dip )
if ( allocated ( dipc )) deallocate ( dipc )
if ( allocated ( dipp )) deallocate ( dipp )
if ( allocated ( dippc )) deallocate ( dippc )
if ( allocated ( dippcm )) deallocate ( dippcm )
if ( allocated ( dipcm )) deallocate ( dipcm )
!if ( allocated ( dqmmm )) deallocate ( dqmmm )
if ( allocated ( dusr )) deallocate ( dusr )
if ( allocated ( dxcdcc )) deallocate ( dxcdcc )
if ( allocated ( dxcdcc_zw )) deallocate ( dxcdcc_zw )
if ( allocated ( dxcv )) deallocate ( dxcv )
if ( allocated ( ewald )) deallocate ( ewald )
!if ( allocated ( ewaldl )) deallocate ( ewaldl )
if ( allocated ( ewaldlr )) deallocate ( ewaldlr )
if ( allocated ( ewaldqmmm )) deallocate ( ewaldqmmm )
if ( allocated ( ewaldsr )) deallocate ( ewaldsr )
if ( allocated ( f0 )) deallocate ( f0 )
if ( allocated ( f3caa )) deallocate ( f3caa )
if ( allocated ( f3cab )) deallocate ( f3cab )
if ( allocated ( f3cac )) deallocate ( f3cac )
if ( allocated ( f3naa )) deallocate ( f3naa )
if ( allocated ( f3nab )) deallocate ( f3nab )
if ( allocated ( f3nac )) deallocate ( f3nac )
if ( allocated ( f3nla )) deallocate ( f3nla )
if ( allocated ( f3nlb )) deallocate ( f3nlb )
if ( allocated ( f3nlc )) deallocate ( f3nlc )
if ( allocated ( f3xca )) deallocate ( f3xca )
if ( allocated ( f3xca_ca )) deallocate ( f3xca_ca )
if ( allocated ( f3xcb )) deallocate ( f3xcb )
if ( allocated ( f3xcb_ca )) deallocate ( f3xcb_ca )
if ( allocated ( f3xcc )) deallocate ( f3xcc )
if ( allocated ( f3xcc_ca )) deallocate ( f3xcc_ca )
if ( allocated ( faca )) deallocate ( faca )
if ( allocated ( fana )) deallocate ( fana )
if ( allocated ( fanl )) deallocate ( fanl )
if ( allocated ( faxc )) deallocate ( faxc )
if ( allocated ( faxc_ca )) deallocate ( faxc_ca )
if ( allocated ( fbias )) deallocate ( fbias )
if ( allocated ( fcoulomb )) deallocate ( fcoulomb )
if ( allocated ( fewald )) deallocate ( fewald )
!if ( allocated ( fewald1 )) deallocate ( fewald1 )
!if ( allocated ( fewald2 )) deallocate ( fewald2 )
if ( allocated ( fharmonic )) deallocate ( fharmonic )
!if ( allocated ( fix )) deallocate ( fix )
if ( allocated ( flrew )) deallocate ( flrew )
!if ( allocated ( flrewl )) deallocate ( flrewl )
if ( allocated ( flrew_qmmm )) deallocate ( flrew_qmmm )
if ( allocated ( Fneb )) deallocate ( Fneb )
if ( allocated ( fotca )) deallocate ( fotca )
if ( allocated ( fotna )) deallocate ( fotna )
if ( allocated ( fotnl )) deallocate ( fotnl )
if ( allocated ( fotxc )) deallocate ( fotxc )
if ( allocated ( fotxc_ca )) deallocate ( fotxc_ca )
if ( allocated ( fragatm )) deallocate ( fragatm )
if ( allocated ( fraggots )) deallocate ( fraggots )
if ( allocated ( fragxyz )) deallocate ( fragxyz )
if ( allocated ( Frec )) deallocate ( Frec )
if ( allocated ( fro )) deallocate ( fro )
if ( allocated ( Fs )) deallocate ( Fs )
if ( allocated ( ft )) deallocate ( ft )
if ( allocated ( ftot )) deallocate ( ftot )
if ( allocated ( ftot1 )) deallocate ( ftot1 )
if ( allocated ( ftot_dftd3 )) deallocate ( ftot_dftd3 )
if ( allocated ( ftot_neb )) deallocate ( ftot_neb )
if ( allocated ( ftotnew )) deallocate ( ftotnew )
if ( allocated ( ftotold )) deallocate ( ftotold )
!if ( allocated ( fumb )) deallocate ( fumb )
if ( allocated ( fvdw )) deallocate ( fvdw )
if ( allocated ( fxcnu )) deallocate ( fxcnu )
if ( allocated ( fxcro )) deallocate ( fxcro )
if ( allocated ( g )) deallocate ( g )
if ( allocated ( g2nu )) deallocate ( g2nu )
if ( allocated ( g2nup )) deallocate ( g2nup )
if ( allocated ( gh_2c )) deallocate ( gh_2c )
if ( allocated ( gh_3c )) deallocate ( gh_3c )
if ( allocated ( gh_atm )) deallocate ( gh_atm )
if ( allocated ( gh_lrew_qmmm )) deallocate ( gh_lrew_qmmm )
if ( allocated ( gh_pp_3c )) deallocate ( gh_pp_3c )
if ( allocated ( gh_pp_atm )) deallocate ( gh_pp_atm )
if ( allocated ( gh_pp_otl )) deallocate ( gh_pp_otl )
if ( allocated ( gh_pp_otr )) deallocate ( gh_pp_otr )
if ( allocated ( gks )) deallocate ( gks )
if ( allocated ( gks_old )) deallocate ( gks_old )
if ( allocated ( gover )) deallocate ( gover )
if ( allocated ( gvhxc )) deallocate ( gvhxc )
if ( allocated ( gvhxcs )) deallocate ( gvhxcs )
if ( allocated ( h )) deallocate ( h )
if ( allocated ( hf_mat )) deallocate ( hf_mat )
if ( allocated ( h_mat )) deallocate ( h_mat )
if ( allocated ( hr_box )) deallocate ( hr_box )
!if ( allocated ( iatomtype )) deallocate ( iatomtype )
if ( allocated ( ideta )) deallocate ( ideta )
if ( allocated ( imass )) deallocate ( imass )
if ( allocated ( initialPosition )) deallocate ( initialPosition )
!if ( allocated ( ipsi22m )) deallocate ( ipsi22m )
if ( allocated ( jatoms_dm )) deallocate ( jatoms_dm )
if ( allocated ( Jialpha )) deallocate ( Jialpha )
if ( allocated ( mask )) deallocate ( mask )
if ( allocated ( neigh_b )) deallocate ( neigh_b )
if ( allocated ( neigh_back )) deallocate ( neigh_back )
!if ( allocated ( neighb_aux )) deallocate ( neighb_aux )
if ( allocated ( neigh_b_classic )) deallocate ( neigh_b_classic )
if ( allocated ( neighb_tot )) deallocate ( neighb_tot )
if ( allocated ( neigh_b_vdw )) deallocate ( neigh_b_vdw )
if ( allocated ( neigh_classic )) deallocate ( neigh_classic )
if ( allocated ( neigh_comb )) deallocate ( neigh_comb )
if ( allocated ( neigh_comj )) deallocate ( neigh_comj )
if ( allocated ( neigh_comm )) deallocate ( neigh_comm )
if ( allocated ( neigh_comn )) deallocate ( neigh_comn )
if ( allocated ( neigh_com_ng )) deallocate ( neigh_com_ng )
if ( allocated ( neigh_j )) deallocate ( neigh_j )
!if ( allocated ( neighj_aux )) deallocate ( neighj_aux )
if ( allocated ( neighj_tot )) deallocate ( neighj_tot )
if ( allocated ( neigh_j_vdw )) deallocate ( neigh_j_vdw )
if ( allocated ( neighn )) deallocate ( neighn )
if ( allocated ( neighn_classic )) deallocate ( neighn_classic )
if ( allocated ( neighn_tot )) deallocate ( neighn_tot )
if ( allocated ( neighn_vdw )) deallocate ( neighn_vdw )
if ( allocated ( neigh_pair_a1 )) deallocate ( neigh_pair_a1 )
if ( allocated ( neigh_pair_a2 )) deallocate ( neigh_pair_a2 )
if ( allocated ( neigh_pair_n1 )) deallocate ( neigh_pair_n1 )
if ( allocated ( neigh_pair_n2 )) deallocate ( neigh_pair_n2 )
if ( allocated ( neighPP_b )) deallocate ( neighPP_b )
if ( allocated ( neighPP_comb )) deallocate ( neighPP_comb )
if ( allocated ( neighPP_comj )) deallocate ( neighPP_comj )
if ( allocated ( neighPP_comm )) deallocate ( neighPP_comm )
if ( allocated ( neighPP_comn )) deallocate ( neighPP_comn )
if ( allocated ( neighPP_j )) deallocate ( neighPP_j )
if ( allocated ( neighPPn )) deallocate ( neighPPn )
if ( allocated ( neighPP_self )) deallocate ( neighPP_self )
if ( allocated ( neigh_self )) deallocate ( neigh_self )
if ( allocated ( neigh_vdw_self )) deallocate ( neigh_vdw_self )
!if ( allocated ( nelectron )) deallocate ( nelectron )
if ( allocated ( nij )) deallocate ( nij )
!if ( allocated ( nitimeshole )) deallocate ( nitimeshole )
if ( allocated ( nowMinusInitialPos )) deallocate ( nowMinusInitialPos )
if ( allocated ( nPP_b )) deallocate ( nPP_b )
if ( allocated ( nPP_j )) deallocate ( nPP_j )
if ( allocated ( nPP_map )) deallocate ( nPP_map )
if ( allocated ( nPPn )) deallocate ( nPPn )
if ( allocated ( nPP_self )) deallocate ( nPP_self )
if ( allocated ( nPPx_b )) deallocate ( nPPx_b )
if ( allocated ( nPPx_j )) deallocate ( nPPx_j )
if ( allocated ( nPPx_map )) deallocate ( nPPx_map )
if ( allocated ( nPPxn )) deallocate ( nPPxn )
if ( allocated ( nPPx_point )) deallocate ( nPPx_point )
if ( allocated ( nPPx_self )) deallocate ( nPPx_self )
!if ( allocated ( npsi22m )) deallocate ( npsi22m )
if ( allocated ( nuxc_total )) deallocate ( nuxc_total )
if ( allocated ( nuxc_3c )) deallocate ( nuxc_3c )
if ( allocated ( hamk )) deallocate ( hamk )
if ( allocated ( getiatom )) deallocate ( getiatom)
if ( allocated ( phidm )) deallocate ( phidm )
!if ( allocated ( psi22m )) deallocate ( psi22m )
!if ( allocated ( Q0 )) deallocate ( Q0 )
!if ( allocated ( Q0_TOT )) deallocate ( Q0_TOT )
!if ( allocated ( Qin )) deallocate ( Qin )
!if ( allocated ( Qin_es )) deallocate ( Qin_es )
!if ( allocated ( Qinmixer )) deallocate ( Qinmixer )
!if ( allocated ( QLowdin_TOT )) deallocate ( QLowdin_TOT )
!if ( allocated ( QLowdin_TOT_es )) deallocate ( QLowdin_TOT_es )
!if ( allocated ( qmmm_struct%Qneutral_TOT )) deallocate ( qmmm_struct%Qneutral_TOT )
!if ( allocated ( qmmm_struct%scf_mchg )) deallocate ( qmmm_struct%scf_mchg )
!if ( allocated ( QMulliken_TOT )) deallocate ( QMulliken_TOT )
!if ( allocated ( Qout )) deallocate ( Qout )
!if ( allocated ( Qoutmixer )) deallocate ( Qoutmixer )
!if ( allocated ( R )) deallocate ( R )
!if ( allocated ( r1_tmp )) deallocate ( r1_tmp )
!if ( allocated ( r2_tmp )) deallocate ( r2_tmp )
!if ( allocated ( r3_tmp )) deallocate ( r3_tmp )
!if ( allocated ( r4_tmp )) deallocate ( r4_tmp )
!if ( allocated ( r5_tmp )) deallocate ( r5_tmp )
if ( allocated ( ratom )) deallocate ( ratom )
if ( allocated ( ratom0 )) deallocate ( ratom0 )
if ( allocated ( ratom2g )) deallocate ( ratom2g )
if ( allocated ( ratom_final )) deallocate ( ratom_final )
if ( allocated ( ratom_frag )) deallocate ( ratom_frag )
if ( allocated ( ratom_frag_save )) deallocate ( ratom_frag_save )
if ( allocated ( ratom_neb )) deallocate ( ratom_neb )
if ( allocated ( ratom_old )) deallocate ( ratom_old )
if ( allocated ( ratom_opt )) deallocate ( ratom_opt )
if ( allocated ( rho )) deallocate ( rho )
if ( allocated ( rhoA )) deallocate ( rhoA )
!if ( allocated ( rhoc )) deallocate ( rhoc )
if ( allocated ( rho_es )) deallocate ( rho_es )
if ( allocated ( rhoij_off )) deallocate ( rhoij_off )
if ( allocated ( rhoi_on )) deallocate ( rhoi_on )
!if ( allocated ( rhom_3c )) deallocate ( rhom_3c )
if ( allocated ( rho_off )) deallocate ( rho_off )
if ( allocated ( rho_old )) deallocate ( rho_old )
if ( allocated ( rho_on )) deallocate ( rho_on )
if ( allocated ( rhopij_off )) deallocate ( rhopij_off )
if ( allocated ( rhop_off )) deallocate ( rhop_off )
if ( allocated ( rhop_on )) deallocate ( rhop_on )
if ( allocated ( rhoPP )) deallocate ( rhoPP )
if ( allocated ( rhoPP_es )) deallocate ( rhoPP_es )
!if ( allocated ( sample1%atom )) deallocate ( sample1%atom )
!if ( allocated ( sample2%atom )) deallocate ( sample2%atom )
if ( allocated ( s_mat )) deallocate ( s_mat )
!if ( allocated ( smatG )) deallocate ( smatG )
if ( allocated ( sm_mat )) deallocate ( sm_mat )
if ( allocated ( sp_mat )) deallocate ( sp_mat )
!if ( allocated ( spmatG )) deallocate ( spmatG )
if ( allocated ( spm_mat )) deallocate ( spm_mat )
if ( allocated ( spVNL )) deallocate ( spVNL )
if ( allocated ( States_total )) deallocate ( States_total )
if ( allocated ( sVNL )) deallocate ( sVNL )
if ( allocated ( symbol )) deallocate ( symbol )
if ( allocated ( tang )) deallocate ( tang )
!if ( allocated ( temp2 )) deallocate ( temp2 )
if ( allocated ( t_mat )) deallocate ( t_mat )
if ( allocated ( tp_mat )) deallocate ( tp_mat )
if ( allocated ( u0vec )) deallocate ( u0vec )
!if ( allocated ( v )) deallocate ( v )
if ( allocated ( vatom )) deallocate ( vatom )
if ( allocated ( vatom_neb )) deallocate ( vatom_neb )
if ( allocated ( vatom_old )) deallocate ( vatom_old )
if ( allocated ( Vbias_mat )) deallocate ( Vbias_mat )
if ( allocated ( Vbiasp_mat )) deallocate ( Vbiasp_mat )
if ( allocated ( vca )) deallocate ( vca )
if ( allocated ( Vcoulomb )) deallocate ( Vcoulomb )
if ( allocated ( Vdip_1c )) deallocate ( Vdip_1c )
if ( allocated ( Vewaldsr )) deallocate ( Vewaldsr )
if ( allocated ( vna )) deallocate ( vna )
if ( allocated ( vnl )) deallocate ( vnl )
if ( allocated ( vxc )) deallocate ( vxc )
if ( allocated ( vxc_1c )) deallocate ( vxc_1c )
if ( allocated ( vxc_3c )) deallocate ( vxc_3c )
if ( allocated ( vxc_ca )) deallocate ( vxc_ca )
if ( allocated ( Vxcnu )) deallocate ( Vxcnu )
if ( allocated ( x0 )) deallocate ( x0 )
if ( allocated ( xdot )) deallocate ( xdot )
if ( allocated ( ximage )) deallocate ( ximage )
if ( allocated ( xmass )) deallocate ( xmass )
if ( allocated ( xl )) deallocate ( xl )
if ( allocated ( ratom )) deallocate ( ratom )
if ( allocated ( symbol )) deallocate ( symbol )
if ( allocated ( imass )) deallocate ( imass )
if ( allocated ( special_k )) deallocate ( special_k)
if ( allocated ( special_k_orig )) deallocate (special_k_orig)
if ( allocated ( scale_k )) deallocate (scale_k)
if ( allocated ( weight_k )) deallocate (weight_k)
if ( allocated ( weight_k_orig )) deallocate (weight_k_orig)
if ( allocated ( special_k )) deallocate (special_k)
if ( allocated ( special_k_orig )) deallocate (special_k_orig)
if ( allocated ( scale_k )) deallocate (scale_k)
if ( allocated ( weight_k )) deallocate (weight_k)
if ( allocated ( weight_k_orig  )) deallocate (weight_k_orig)

if ( allocated ( getiatom )) deallocate ( getiatom)
if ( allocated ( getissh )) deallocate ( getissh)
if ( allocated ( getlssh )) deallocate ( getlssh)
if ( allocated ( getmssh )) deallocate ( getmssh )
if ( allocated (bbnkre_o)) deallocate ( bbnkre_o)
if ( allocated (blowre_o)) deallocate ( blowre_o)

if ( allocated (nelectron )) deallocate (nelectron )
if ( allocated (Qin )) deallocate (Qin )
if ( allocated (Qinmixer)) deallocate (Qinmixer )
if ( allocated (QLowdin_TOT )) deallocate (QLowdin_TOT )
if ( allocated (QMulliken_TOT )) deallocate (QMulliken_TOT )
if ( allocated (Qout)) deallocate (Qout )
if ( allocated (Q_partial)) deallocate (Q_partial)
if ( allocated (Qoutmixer)) deallocate (Qoutmixer )
if ( allocated (dq)) deallocate (dq )
if ( allocated (Q0_TOT)) deallocate (Q0_TOT )
if ( allocated (Qin_es)) deallocate (Qin_es )
if ( allocated (QLowdin_TOT_es )) deallocate (QLowdin_TOT_es )
if ( allocated (Qout_es)) deallocate (Qout_es )
if ( allocated (qaux)) deallocate (qaux )
if ( allocated (Fv )) deallocate (Fv)
if ( allocated (Xv )) deallocate (Xv)
if ( allocated (delX )) deallocate (delX)
if ( allocated (delF )) deallocate (delF)
if ( allocated (r2_sav )) deallocate (r2_sav)
!if ( allocated (atompos)) deallocate (atompos)

if ( allocated (mesh_wf )) deallocate (mesh_wf )
if ( allocated (drr_wf )) deallocate (drr_wf )
if ( allocated (rmax_wf )) deallocate (rmax_wf )
if ( allocated (rr_wf )) deallocate (rr_wf )
if ( allocated (wf_spline )) deallocate (wf_spline )
if ( allocated (wf )) deallocate (wf )

! vneutral part
if ( allocated (mesh_na )) deallocate (mesh_na )
if ( allocated (rmax_na )) deallocate (rmax_na )
if ( allocated (drr_na )) deallocate (drr_na )
if ( allocated (rr_na )) deallocate (rr_na )
if ( allocated (vnna_spline )) deallocate ( vnna_spline )
if ( allocated (vnna )) deallocate (vnna )

if ( allocated (e2r )) deallocate (e2r )
if ( allocated (am2rc)) deallocate (am2rc)
if ( allocated (ram2rc)) deallocate (ram2rc)
if ( allocated (ratom2g)) deallocate (ratom2g)
if ( allocated (e2n )) deallocate (e2n )
if ( allocated (n2e)) deallocate (n2e)
if ( allocated (vnaG )) deallocate (vnaG )
if ( allocated (drhoG )) deallocate (drhoG )
if ( allocated (rhoG0 )) deallocate (rhoG0 )
if ( allocated (vcaG )) deallocate (vcaG )
if ( allocated (vxcG )) deallocate (vxcG )

if ( allocated ( rho_in )) deallocate ( rho_in ) 
if ( allocated ( rho_out )) deallocate ( rho_out ) 
if ( allocated ( mwe )) deallocate ( mwe ) 
if ( allocated ( drwe )) deallocate ( drwe ) 
if ( allocated ( mwe )) deallocate ( mwe ) 
if ( allocated ( drwe )) deallocate ( drwe ) 
if ( allocated ( dngcof )) deallocate ( dngcof ) 
if ( allocated ( E_KS )) deallocate ( E_KS ) 
if ( allocated ( DOS_total )) deallocate ( DOS_total ) 

if ( allocated ( gcoefficientsVNA )) deallocate ( gcoefficientsVNA )
if ( allocated ( alphaVNA )) deallocate ( alphaVNA )
if ( allocated ( nalphaVNA )) deallocate ( nalphaVNA )
if ( allocated ( gcoefficientsN )) deallocate ( gcoefficientsN )
if ( allocated ( alphaN )) deallocate ( alphaN )
if ( allocated ( nalphaN )) deallocate ( nalphaN )
if ( allocated ( gcoefficientsPSI )) deallocate ( gcoefficientsPSI )
if ( allocated ( alphaPSI )) deallocate ( alphaPSI )
if ( allocated ( nalphaPSI )) deallocate ( nalphaPSI )
if ( allocated ( gcoefficientsPSIS )) deallocate ( gcoefficientsPSIS )
if ( allocated ( alphaPSIS )) deallocate ( alphaPSIS )
if ( allocated ( nalphaPSIS )) deallocate ( nalphaPSIS )
if ( allocated ( gcoefficientsVNA_SH )) deallocate ( gcoefficientsVNA_SH )
if ( allocated ( alphaVNA_SH )) deallocate ( alphaVNA_SH )
if ( allocated ( nalphaVNA_SH )) deallocate ( nalphaVNA_SH )
if ( allocated ( R_na )) deallocate ( R_na )
if ( allocated (blowim)) deallocate ( blowim )
if ( allocated (bbnkre)) deallocate ( bbnkre)
if ( allocated ( E_KS)) deallocate ( E_KS)
if ( allocated (Uisigma )) deallocate ( Uisigma)
if ( allocated (Jijsigma )) deallocate ( Jijsigma)

if ( allocated (xintegral_2c )) deallocate (xintegral_2c )
!if ( allocated ( )) deallocate ( )
!if ( allocated ( )) deallocate ( )
end
