! Copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Institute of Physics, Czech Republic - Pavel Jelinek

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

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


! initbasics.f90
! Program Description
! ===========================================================================
! Initialize basic variables controlling flow of the code
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
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine initbasics ()

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

! Argument Declaration and Description
! ===========================================================================
! Input

! Output

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer in1
   integer icount
   integer isorp
   integer ideriv
   integer issh
   integer numorbPP_max
   integer numorb

   real distance
   real, dimension (3) :: vector

   logical file_exists

! Allocate arrays
! ===========================================================================

! Procedure
! ===========================================================================

! ***************************************************************************
! ===========================================================================
! ---------------------------------------------------------------------------
!                             W E L C O M E !
! ---------------------------------------------------------------------------
! ===========================================================================
        call welcome
! ===========================================================================
! ---------------------------------------------------------------------------
!                  I N I T I A L I Z E    P A R A M E T E R S
!                      A N D   R E A D   D A T A F I L E S
! ---------------------------------------------------------------------------
! ===========================================================================
! Initialize some constants
        call initconstants (sigma, sigmaold, scf_achieved)
! Read in the diagnostics parameters. This allows the user to turn
! particular interactions on or off.
        call diagnostics (ioff2c, ioff3c, itestrange, testrange)

! Read the parameter file - fireball.param
        call readparam ()

! Read the info.dat file.  Allocate lsshPP, etc. inside!
        call readinfo ()
!CHROM
! Cutoff dist. of the classical interaction is used in neighbor routine
		if( iclassicMD > 0 ) call readdata_classicMD ()
!END CHROM


! Read data from the basis file - XXX.bas.
        open (unit = 69, file = basisfile, status = 'old')
        read (69, *) natoms
        close (unit = 69)



! Initialize aux. variable
        if (nstepi .eq. 1) then
         T_average = T_initial
         T_previous = 0.0d0
         time = 0.0d0
        end if

! Allocate more arrays.
        allocate (degelec (natoms))
        allocate (imass (natoms))
        allocate (ratom (3, natoms))
        allocate (nowMinusInitialPos (3, natoms))
        allocate (initialPosition (3, natoms))
        allocate (vatom (3, natoms))
        allocate (symbol (natoms))
        allocate (xmass (natoms))
        allocate (ximage (3, natoms))
        allocate (mask (3,natoms))
        mask = 1.0d0
        ximage = 0.0d0
        call readbasis (nzx, imass)
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



! Call make_munu. This routine determines all of the non-zero matrix elements
! for the two- and three-center matrix elements.  These non-zero matrix
! elements are determined based on selection rules.
        call make_munu (nspecies)

! We now call make_munuPP.
        call make_munuPP (nspecies)

! We now call make_munuS
        call make_munuS (nspecies)

! JIMM: New long-range theory idipole = 1
!  We now call make_munuDipY and make_munuDipX
        if (idipole .eq. 1) then
         call make_munuDipY (nspecies)
         call make_munuDipX (nspecies)
        end if

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

! Count the maximum number of orbital interactions between any given two atoms.
        numorb_max = 0
        do in1 = 1, nspecies
         numorb = 0
         do issh = 1, nssh(in1)
          numorb = numorb + 2*lssh(issh,in1) + 1
         end do
         if (numorb .gt. numorb_max) numorb_max = numorb
        end do

! We originally use the neutral atom charges from the info.dat file.
! If this is a restart run, or the user wishes to start with charges other
! than the neutral atom stuff, then information is read in from a charge file.
        if (iKS .eq. 1) then
         call initcharges_KS (natoms, nspecies, itheory, ifixcharge, symbol)
        else
         call initcharges (natoms, nspecies, itheory, ifixcharge, symbol)
        endif

! Calculate isorpmax and ideriv_max
        isorpmax = 0
        if (itheory .eq. 1) then
         do in1 = 1, nspecies
          isorpmax = max(isorpmax,nssh(in1))
         end do
        end if
! this is only for density-oslxc, in future should be the same as above
! depending on the harris-oslxc
        isorpmax_xc = 0
        do in1 = 1, nspecies
           isorpmax_xc = max(isorpmax_xc,nssh(in1))
        end do

        ideriv_max = 0
        if (itheory .eq. 1) ideriv_max = 6

! Set up the index field ind2c:
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
        if (itheory_xc .eq. 1 .or. itheory_xc .eq. 2) then
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
!dani.JOM.jel-fr-end

        end if
        interactions2c_max = icount

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
 end subroutine initbasics
