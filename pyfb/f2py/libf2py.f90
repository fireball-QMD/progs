subroutine start()
  call cpu_time (time_begin)
  call initbasics ()
  call readdata ()
  call main_loop ()
  call cpu_time (time_end)
  write (*,*) ' FIREBALL RUNTIME : ',time_end-time_begin,'[sec]'
end

subroutine f2py_options(iclusteraux)
  use options
  integer,intent(in)::iclusteraux
  icluster=iclusteraux
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

  print*,'aux',nuczaux,natoms,zaux
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
  integer i,nucz, iatom,j
  print*,'raux',raux,zaux
  print*,raux(2,1)
  allocate (ratom (3, natoms))
  do iatom = 1, natoms
  print*,raux(iatom,:)
  end do
  do iatom = 1, natoms
    ratom(:,iatom) = raux(iatom,:)
  end do
  
  do iatom = 1, natoms
     print*,imass(iatom),ratom(:,iatom)
  end do
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
   print*,'f2py_fdatalocation = ',f2py_fdatalocation
   fdatalocation = f2py_fdatalocation
   !verbosity=0
   !cambiamos readinfo para cargar Fdata completa sin necesidad de leer 
   !las posiciones en el bas
   !call readinfo ()
   call readinfoall()
   if( iclassicMD > 0 ) call readdata_classicMD ()
   if (nstepi .eq. 1) then
     T_average = T_initial
     T_previous = 0.0d0
     time = 0.0d0
   end if
   !call readbasis (nzx, imass)   
   call readquench (iquench, dt, energy_tol, force_tol, iensemble, T_initial, T_want, taurelax)
   if (iendtemp .eq. 1 .and. iquench .eq. 0) then
     if (T_initial .eq. T_final) then
       write (*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       write (*,*) 'T_initial = T_final'
       write (*,*) 'If you wanted this then set iendtemp = 0'
       write (*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       stop
     end if
     T_increment = (T_final - T_initial)/nstepf
   else if (iendtemp .eq. 0) then
     T_increment = 0.00
   end if
   if (iendtemp .eq. 1 .and. iquench .ne. 0) then
     write (*,*) 'STOPPING'
     write (*,*) 'iendtemp = 1  and iquench != 0'
     stop
   end if
   if (ivdw .eq. 1) call readvdw (nspecies, symbolA, ivdw)
   call make_munu (nspecies)
   call make_munuPP (nspecies)
   call make_munuS (nspecies)
   if (idipole .eq. 1) then
     call make_munuDipY (nspecies)
     call make_munuDipX (nspecies)
   end if
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

! check if we need the grid
   igrid = 0
   if (iwrtden .eq. 1) igrid = 1
   if (iwrtewf .eq. 1) igrid = 1
   if (iks .eq. 1) igrid = 1
   if (iwrtdipole .eq. 1) igrid = 1
   if (igrid .eq. 1) call readgrid (iwrtewf)
end subroutine


subroutine f2py_init() !zauxf2py) !,pos,Zin)
   use options
   use configuration
   use interactions
   use outputs
   use optimization
   use md
   use charges
   use barrier
   use energy
   use neighbor_map
   use forces
   use mpi_main
   implicit none
!   integer,intent(in):: N
!   real, intent(in),dimension(N,3) :: pos
!   integer, intent(in),dimension(Nf2py) :: zauxf2py
   integer iatom
   integer in1
   integer issh

   real distance
  
  allocate (degelec (natoms))
!  allocate (imass (natoms))
!  allocate (ratom (3, natoms))
  allocate (nowMinusInitialPos (3, natoms))
  allocate (initialPosition (3, natoms))
  allocate (vatom (3, natoms))
  allocate (symbol (natoms))
  allocate (xmass (natoms))
  allocate (ximage (3, natoms))
  allocate (mask (3,natoms))


 !call readbasis (nzx, imass)  lo hacemos antes con f2py_loadbas

  mask = 1.0d0
  ximage = 0.0d0
  ishiftO = 0
  do iatom = 1, natoms
   distance = ratom(1,iatom)**2 + ratom(2,iatom)**2 + ratom(3,iatom)**2
   distance = sqrt(distance)
   if (distance .lt. 1.0d-4) ishiftO = 1
  end do
  call initmasses (natoms, symbol, smass, xmass)
  if (iwrtdos.ge.1 .or. iwrtatom .ge. 1) call readdos ( )
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

  if (iKS .eq. 1) then
   call initcharges_KS (natoms, nspecies, itheory, ifixcharge, symbol)
  else
   call initcharges (natoms, nspecies, itheory, ifixcharge, symbol)
  endif
  call get_info_orbital (natoms)

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

  call initneb (natoms, nspecies, imass, nzx, ratom)

  allocate (xdot (0:5, 3, natoms)) ! Initialized below
  if (nstepi .ne. 1) then
    call restart (nstepi, dt, acfile, xvfile, T_average, T_previous, time)
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

  call initatomicE (natoms, etotatom, imass, atomic_energy)

  call setgear

  if (iensemble .eq. 2) call initNH(natoms,T_want)

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
  if (nstepi .eq. 1) then
    etotnew = 0.0d0
  if (.not. allocated (ftotnew)) allocate (ftotnew (3, natoms))
    ftotnew = 0.0d0
  end if

  ! jel-grid
  ! Initialize the density matrix
  if (igrid .eq. 1 ) then
    call neighbors (nprocs, my_proc, iordern, icluster, iwrtneigh, ivdw)
    call neighbors_pairs(icluster)
    write (*,*) 'Initialize density matrix'
    call initdenmat (natoms)
  elseif ( iclassicMD > 0 .and. igrid /= 1 )then
    call neighbors (nprocs, my_proc, iordern, icluster, iwrtneigh, ivdw)
    call neighbors_pairs(icluster)
  endif
  if (itdse .eq. 1) then
  write (*,*) ' Read TD parameters'
  call readtdse ()
  write (*,*) ' Allocate TD-matrices'
  call allocate_tdse ()
  endif

end


subroutine loadlvs()

   use kpoints
  ! Read data from the lattice vectors file - XXX.lvs.
  ! Set up the boxes surrounding the central unit cell.
!  call readlvs (lvsfile, a1vec, a2vec, a3vec, icluster, rescal)

  ! Define volume of unit cell
!  call cross (a2vec, a3vec, vector)
!  Vouc=a1vec(1)*vector(1)+a1vec(2)*vector(2)+a1vec(3)*vector(3)
!  call initboxes (1)

  ! Get kpoints for Brillouin zone integration or bandstructure calculation
!  if(iclassicMD /= 1) then !not usefull with empirical potentials
!  call getkpoints(icluster, Vouc, a1vec, a2vec, a3vec, &
!   &       lvsfile, basisfile, iquench, ireducekpts, rescal)
!  end if
  if (iordern .eq. 1 .and. nkpoints .ne. 1) then
    write (*,*) ' Order-N method only works with one k-point! '
    write (*,*) ' Sorry, you must abort your run! '
    stop
  end if
end





subroutine loadfragments()
  call readfragments()
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



