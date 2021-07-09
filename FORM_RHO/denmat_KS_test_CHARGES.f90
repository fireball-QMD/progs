! cp denmat_KS_test_CHARGES.f90  denmat_KS.f90 
! and make
! For testing purposes only 
! For testing purposes only 
! For testing purposes only 
! For testing purposes only 
! For testing purposes only 
! copyright info:f
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

 
! denmat.f90
! Program Description
! ===========================================================================
!       This routine calculates the density matrices and the band-structure
! energy, as well as similar density matrices.
!
! ===========================================================================
! Code rewritten by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! (modified by P. Jelinek; May 2005 Utah)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine denmat_KS (ifixcharge, iqout, icluster, iwrtefermi, tempfe, &
  &                           ebs)
        use charges
        use configuration
        use constants_fireball
        use density
        use dimensions
        use interactions
        use neighbor_map
        use kpoints
        use hartree_fock
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icluster
        integer, intent (in) :: ifixcharge
        integer, intent (in) :: iqout
        integer, intent (in) :: iwrtefermi

        real, intent (in) :: tempfe

! Output
        real, intent (out) :: ebs
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iband
        integer jband
        integer ikpoint
        integer imu, inu
        integer ineigh
        integer in1, in2
        integer iorbital
        integer issh, jssh
        integer jatom
        integer jneigh
        integer mqn
        integer mbeta
        integer mmu
        integer noccupy
        integer nnu
        real    Qtot, Qtot1, Qtot2
        real    dip_x, dipQin_x, dipTot_x, dipProy_x, dipIntra_x
        real    dip_y, dipQin_y, dipTot_y, dipProy_y, dipIntra_y
        real    dip_z, dipQin_z, dipTot_z, dipProy_z, dipIntra_z
        real    dip_tot, dip_proy, dipQin_tot, dipProy_tot, dipTot_tot, dipIntra_tot
        real, dimension(3) :: Rbc,u21
        integer, dimension (norbitals) :: ioccupy
        integer, dimension (norbitals, nkpoints) :: ioccupy_k
 
        real aux1, aux2, aux3
        real qmulliken_mix
        real deltae

        real dot
        real gutr
        real pcharge
        real ztest
        real norm 
        real y
        real, dimension (norbitals, nkpoints) :: foccupy
        real, dimension (numorb_max, natoms) :: QMulliken
        real, dimension (3) :: vec, r1, r2, r21
        real, parameter ::  Debye = 0.208194
 
        complex ai
        complex phase, phasex
        complex step1, step2
 
        logical read_occupy

        integer issh1, mu_min, mu_max, l, inumorb
        integer :: info, lwork
        integer, dimension(nssh_tot) :: ipiv
        integer, dimension(100) :: work
        integer :: beta, alpha, ialp, ina, matom
        real,dimension(nssh_tot+1,nssh_tot+1) :: M
        real,dimension(1,nssh_tot+1) :: B
        real auxgS
        real Ntot
        real, dimension(nsh_max,natoms) :: Qout_est
        real, dimension(natoms) :: Qest_TOT

! A bunch of memory to be used in many ways
        integer jnu,jmu
        complex*16, dimension (:, :), allocatable :: xxxx
        complex*16, dimension (:, :), allocatable :: yyyy
        complex*16 a0
        complex*16 a1
        real, dimension (3) :: k_temp

        real,dimension(nssh_tot,nssh_tot) :: A
        real,dimension(nssh_tot) :: c, SQ ! carga
        real,dimension(nssh_tot) :: LB, UB, nalpha
        real :: diff_err,Ep2

! Procedure
! ===========================================================================
! Initialize some things
        ai = cmplx(0.0d0,1.0d0)
! jel: where you set rho to zero ???
        rhoPP = 0.0d0
        rho = 0.0d0

        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' ****************************************************** '
        write (*,*) '  '
        write (*,*) '                 -  Welcome to denmat_G -               '
        write (*,*) '  '
        write (*,*) ' ****************************************************** '


! ****************************************************************************
!
!                     C H A R G E    O C C U P A T I O N S
! ****************************************************************************
! If there exists a file 'OCCUPATION', then make a list of energy eigenvalues
! you want to 'unoccupy', and give the energy shift to put the eigenvalue
! above E-fermi.
        inquire (file = 'OCCUPATION', exist = read_occupy)
        if (read_occupy) then
 
! Open the file and read information.
         open (unit = 22, file = 'OCCUPATION', status = 'old')
         write (*,*) '  '
         write (*,*) ' Reading from the OCCUPATION file! '
         read (22,*) noccupy
         if (noccupy .gt. norbitals) then
          write (*,*) ' noccupy > norbitals: from OCCUPATION file. '
          stop
         end if
         do imu = 1, noccupy
          read (22,*) iband, deltae
          eigen_k(iband,1:nkpoints) = eigen_k(iband,1:nkpoints) + deltae
          ioccupy(imu) = iband
         end do
         close (unit = 22)
 
! If you are monkeying with the occupation it is probable that you would like
! to know about charge localization on the snuffed sites, how much
! eigenvector there is in orbital mu, etc.
         write (*,*) '  '
         do iorbital = 1, noccupy
          iband = ioccupy(iorbital)
          write (*,*) '  '
          write (*,*) ' Band # ', iband,' shifted E = ', eigen_k(iband,1)
          write (*,*) ' **************************************************** '
 
! Loop over the atoms and the neighbors
          do iatom = 1, natoms
           in1 = imass(iatom)
           pcharge = 0.0d0
           do ineigh = 1, neighn(iatom)
            mbeta = neigh_b(ineigh,iatom)
            jatom = neigh_j(ineigh,iatom)
            in2 = imass(jatom)
            vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
 
! Loop over the special k points.
            do ikpoint = 1, nkpoints
             dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2) &
     &                                         + special_k(3,ikpoint)*vec(3)
             phase = cmplx(cos(dot),sin(dot))*weight_k(ikpoint) 
! Loop over all bands
             if (icluster .ne. 1) then
              do imu = 1, num_orb(in1)
               mmu = imu + degelec(iatom)
               step1 = phase*(bbnkre(mmu,iband,ikpoint)                      &
     &                        - ai*bbnkim(mmu,iband,ikpoint))
               do inu = 1, num_orb(in2)
                nnu = inu + degelec(jatom)
                step2 = step1*(bbnkre(nnu,iband,ikpoint)                     &
     &                         + ai*bbnkim(nnu,iband,ikpoint))
                gutr = real(step2)
                pcharge = pcharge + gutr*s_mat(imu,inu,ineigh,iatom)
 
! Finished loop over bands
               end do
              end do
             else
              do imu = 1, num_orb(in1)
               mmu = imu + degelec(iatom)
               step1 = phase*bbnkre(mmu,iband,ikpoint) 
               do inu = 1, num_orb(in2)
                nnu = inu + degelec(jatom)
                step2 = step1*bbnkre(nnu,iband,ikpoint)  
                gutr = real(step2)
                pcharge = pcharge + gutr*s_mat(imu,inu,ineigh,iatom) 
! Finished loop over bands
               end do
              end do
             end if
 
! Finished loop over kpoints
            end do
 
! Finish remaining loops
           end do
           write (*,100) iband, iatom, pcharge
          end do
         end do
        end if
 
! ****************************************************************************
! Get the Fermi energy.
        call fermie (norbitals, ztot, eigen_k, efermi, ioccupy_k, foccupy)
!        write (*,*) ' Fermi Level = ', efermi
 
        if (iwrtefermi .eq. 1) then
         write (*,*) '  '
         write (*,*) ' We write out the occupancies of the levels from the  '
         write (*,*) ' subroutine fermie ----- ioccupy_k! '
         do ikpoint = 1, nkpoints
          write (*,*) '  '
          write (*,*) ' ------ fermi ioccupy_k for k-point = ', ikpoint
          do iband = 1, norbitals_new
           write (*,200) iband, ioccupy_k(iband,ikpoint)
          end do
         end do
 
         write (*,*) '  '
         write (*,*) ' We write out the electron fraction in each level - from '
         write (*,*) ' the subroutine fermie. '
         do ikpoint = 1, nkpoints
          write (*,*) '  '
          write (*,*) ' ------ fermi foccupy for k-point = ', ikpoint
          do iband = 1, norbitals_new
           write (*,201) iband, foccupy(iband,ikpoint)
          end do
         end do
         write (*,*) '  '
        end if
 
 
! ****************************************************************************
!
!                      C O M P U T E    D E N S I T I E S
! ****************************************************************************
! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
        pcharge = 0.0d0
        do iatom = 1, natoms
         in1 = imass(iatom)
         do ineigh = 1, neighn(iatom)
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom)
          vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
 
! Loop over the special k points.
          do ikpoint = 1, nkpoints
           dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2)   &
     &                                       + special_k(3,ikpoint)*vec(3)
           phasex = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)*spin
 
! Loop over all bands
           if (icluster .ne. 1) then
            do iband = 1, norbitals_new
             if (ioccupy_k(iband,ikpoint) .ne. 0) then
              phase = phasex*foccupy(iband,ikpoint)
              do imu = 1, num_orb(in1)
               mmu = imu + degelec(iatom)
               step1 = phase*(bbnkre(mmu,iband,ikpoint)                      &
     &                        - ai*bbnkim(mmu,iband,ikpoint))
               do inu = 1, num_orb(in2)
                nnu = inu + degelec(jatom)
                step2 = step1*(bbnkre(nnu,iband,ikpoint)                     &
     &                         + ai*bbnkim(nnu,iband,ikpoint))
                gutr = real(step2)
                pcharge = pcharge + gutr*s_mat(imu,inu,ineigh,iatom)
! Finally the expressions.........
                rho(imu,inu,ineigh,iatom) = rho(imu,inu,ineigh,iatom) + gutr
                cape(imu,inu,ineigh,iatom) =                                 &
     &           cape(imu,inu,ineigh,iatom) + eigen_k(iband,ikpoint)*gutr
               end do
              end do
             end if
 
! Finish loop over bands.
            end do
           else
            do iband = 1, norbitals_new
             if (ioccupy_k(iband,ikpoint) .ne. 0) then
              phase = phasex*foccupy(iband,ikpoint)
              do imu = 1, num_orb(in1)
               mmu = imu + degelec(iatom)
               step1 = phase*bbnkre(mmu,iband,ikpoint)
               do inu = 1, num_orb(in2)
                nnu = inu + degelec(jatom)
                step2 = step1*bbnkre(nnu,iband,ikpoint)  
                gutr = real(step2)
                pcharge = pcharge + gutr*s_mat(imu,inu,ineigh,iatom)
! Finally the expressions.........
                rho(imu,inu,ineigh,iatom) = rho(imu,inu,ineigh,iatom) + gutr
                cape(imu,inu,ineigh,iatom) =                                 &
     &           cape(imu,inu,ineigh,iatom) + eigen_k(iband,ikpoint)*gutr
               end do
              end do
             end if
 
! Finish loop over bands.
            end do
           end if

! Finish loop over k-points.
          end do


! Finish loop over atoms and neighbors.
         end do
        end do

        write (*,*) 'pcharge =', pcharge


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                   M A T R I X     D E N S I T Y    
!                          PP-neighbors
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Loop over all atoms iatom in the unit cell, and then over all its neighbors.
        do iatom = 1, natoms
         in1 = imass(iatom)
         do ineigh = 1, neighPPn(iatom)
          mbeta = neighPP_b(ineigh,iatom)
          jatom = neighPP_j(ineigh,iatom)
          in2 = imass(jatom)
          vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
! Loop over the special k points.
          do ikpoint = 1, nkpoints
           dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2)   &
     &                                       + special_k(3,ikpoint)*vec(3)
           phasex = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)*spin
 
! Loop over all bands
           if (icluster .ne. 1) then
            do iband = 1, norbitals_new
             if (ioccupy_k(iband,ikpoint) .ne. 0) then
              phase = phasex*foccupy(iband,ikpoint)
              do imu = 1, num_orb(in1)
               mmu = imu + degelec(iatom)
               step1 = phase*(bbnkre(mmu,iband,ikpoint)                      &
     &                        - ai*bbnkim(mmu,iband,ikpoint))
               do inu = 1, num_orb(in2)
                nnu = inu + degelec(jatom)
                step2 = step1*(bbnkre(nnu,iband,ikpoint)                     &
     &                         + ai*bbnkim(nnu,iband,ikpoint))
                gutr = real(step2)
! Finally the expressions.........
                rhoPP(imu,inu,ineigh,iatom) = rhoPP(imu,inu,ineigh,iatom) &
     &                          + gutr

               end do
              end do
             end if
 
! Finish loop over bands.
            end do
           else
            do iband = 1, norbitals_new
             if (ioccupy_k(iband,ikpoint) .ne. 0) then
              phase = phasex*foccupy(iband,ikpoint)
              do imu = 1, num_orb(in1)
               mmu = imu + degelec(iatom)
               step1 = phase*bbnkre(mmu,iband,ikpoint)
               do inu = 1, num_orb(in2)
                nnu = inu + degelec(jatom)
                step2 = step1*bbnkre(nnu,iband,ikpoint)  
                gutr = real(step2)
 
! Finally the expressions.........
                rhoPP(imu,inu,ineigh,iatom) = rhoPP(imu,inu,ineigh,iatom) &
     &                          + gutr
               end do
              end do
             end if
 
! Finish loop over bands.
            end do
           end if

! Finish loop over k-points.
          end do

! Finish loop over atoms and neighbors.
         end do
        end do

        if (ifixcharge .eq. 1) then

          do iatom = 1, natoms
            if (iqout .eq. 1 .or. iqout .eq. 3) QLowdin_TOT(iatom) = 0.0d0
            if (iqout .eq. 2 .or. iqout .eq. 4) QMulliken_TOT(iatom) = 0.0d0
            do issh = 1, nssh(imass(iatom))
              Qout(issh,iatom) = Qin(issh,iatom)
              if (iqout .eq. 1 .or. iqout .eq. 3) then
                QLowdin_TOT(iatom) =  QLowdin_TOT(iatom) +Qin(issh,iatom)
              end if
              if (iqout .eq. 2 .or. iqout .eq. 4) then
                QMulliken_TOT(iatom) = QMulliken_TOT(iatom) +Qin(issh,iatom)
              end if
            end do
          end do

        else !ifixcharge

          ! iqout=1 ! Lowdin

          Qout = 0.0d0
          QLowdin_TOT = 0.0d0
          if(.not. allocated(blowre)) allocate (blowre (norbitals, norbitals, nkpoints))
          do ikpoint = 1, nkpoints
            k_temp(:) = special_k(:,ikpoint)
            !------------------iqout
            call calc_blowreim (1, icluster, ikpoint, nkpoints, k_temp)
          end do
          call LOWDIN_CHARGES(ioccupy_k,foccupy)
          call write_charges_shell('CHARGES_LOWD_KS')
          call write_dipole('_LOWD_KS')
          
          ! iqout=3 !NPA

          Qout = 0.0d0
          QLowdin_TOT = 0.0d0
          do ikpoint = 1, nkpoints
            k_temp(:) = special_k(:,ikpoint)
            !------------------iqout
            call calc_blowreim (3, icluster, ikpoint, nkpoints, k_temp)
          end do
          call LOWDIN_CHARGES(ioccupy_k,foccupy)
          call write_charges_shell('CHARGES_NAP_KS')
          call write_dipole('_NAP_KS')

          ! iqout =  2 ! Mulliken

          Qout = 0.0d0
          QMulliken_TOT = 0.0d0
          call MULLIKEN_CHARGES()
          call write_charges_shell('CHARGES_MULL_KS') 
          call write_dipole('_MULL_KS')

          ! iqout = 4 ! Mulliken-dipole

          Qout = 0.0d0                                                   
          QMulliken_TOT = 0.0d0 
          call MULLIKEN_DIPOLE_CHARGES()          
          call write_charges_shell('CHARGES_MUDI_KS') 
          call write_dipole('_MUDI_KS')

    
          !if (iqout .eq. 5)  call charges_L2(Kscf,0) !DIAGONAL CHARGES           
       
          ! iqout =  6 ! Stationary
!          call STATIONARY_CHARGES()
!          open(unit = 111, file = 'CHARGES_SC_KS', status = 'unknown')
!            do iatom = 1, natoms
!              write(111,'(2x, i4)',advance='no') iatom
!              Ntot=0
!              do imu = 1, num_orb(imass(iatom))
!                Ntot = Ntot + Qout(imu,iatom)
!              end do
!              write(111,'(f10.6)',advance='no') Ntot
!              do imu = 1, num_orb(imass(iatom))
!                write(111,'(f10.6)',advance='no') Qout(imu,iatom)
!              end do
!              write(111,*)
!            end do
!          close(111)

        end if !ifixcharge


! ****************************************************************************
!
!  C O M P U T E    B A N D - S T R U C T U R E    E N E R G Y
! ****************************************************************************
! Compute ebs, the band structure energy.
        ebs = 0.0d0
        ztest = 0.0d0
        do ikpoint = 1, nkpoints
         do iorbital = 1, norbitals_new
          if (ioccupy_k(iorbital,ikpoint) .eq. 1) then
           ebs = ebs + weight_k(ikpoint)*spin*eigen_k(iorbital,ikpoint)      &
     &                 *foccupy(iorbital,ikpoint)
           ztest = ztest + weight_k(ikpoint)*spin*foccupy(iorbital,ikpoint)
          end if
         end do
        end do
 
! Test to make sure we get the proper number of states.
        if (abs(ztest - ztot) .gt. 1.0d-02) then
         write (*,*) ' *************** error *************** '
         write (*,*) ' ztest = ', ztest, ' ztot = ', ztot
         write (*,*) ' In denmat.f - ztest .ne. ztot! '
         stop
        end if
 
! Format Statements
! ===========================================================================
100     format (2x, 2i4, f8.4)
200     format (' Band n = ', i4, ' k-points: ioccupy = ', i2)
201     format (' Band n = ', i4, ' foccupy = ', f6.3)
300     format (2x, ' This is band number: ',2x, i6)
301     format (2x, i4, f10.6)
110     format (2x, 4f10.6)
!120     format (2x, i4, <norbitals>f10.6) 
!121     format (2x, i4,f10.6,<norbitals>f10.6) 
400     format (2x, 'Qmull =',10f10.6)
444     format (a7,4f10.4)
445     format (a2,4f10.4)
        return
      end subroutine denmat_KS
       

! ===========================================================================


       subroutine calc_blowreim (iqout, icluster, ikpoint, nkpoints, sks)
        use density
        use configuration
        use dimensions
        use interactions
        use neighbor_map
        use charges
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iqout
        integer, intent (in) :: icluster
        integer, intent (in) :: ikpoint
        integer, intent (in) :: nkpoints
        real, intent (in), dimension (3) :: sks

! Output


! Local Parameters and Data Declaration
! ===========================================================================
        real*8, parameter :: overtol = 1.0d-4

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer imu
        integer info
        integer inu
        integer in1, in2
        integer ineigh
        integer jatom
        integer jmu
        integer jnu
        integer mbeta
        integer mineig
        integer lm
        integer issh

        real sqlami

        real*8, dimension (norbitals) :: eigen

        real*8, dimension (norbitals) :: slam

        real*8 a0
        real*8 a1
        real*8 magnitude

! A bunch of memory to be used in many ways
        real*8, dimension (:, :), allocatable :: xxxx
        real*8, dimension (:, :), allocatable :: yyyy
        real*8, dimension (:, :), allocatable :: zzzz
        real*8, dimension (:, :), allocatable, save :: sm12_save
!NPA
        real*8, dimension (:, :), allocatable :: ssss
        real*8, dimension (:), allocatable :: ww

! work vector for diaganlization
        real*8, allocatable, dimension (:) :: work
        integer, allocatable, dimension (:) :: iwork
        integer lwork, liwork

! Procedure
! ===========================================================================
! Initialize some things
        magnitude = sqrt(sks(1)**2 + sks(2)**2 + sks(3)**3)
        if (magnitude .gt. 1.0d-3) then
         write (*,*) ' *************** WARNING *************** '
         write (*,*) ' *************** WARNING *************** '
         write (*,*) ' *************** WARNING *************** '
         write (*,*) ' You have chosen to do a calculation with k-points other '
         write (*,*) ' than just the gamma point.  Unfotunately, this version '
         write (*,*) ' of kspace calls the LAPACK subroutine that assumes '
         write (*,*) ' only real matrices (i.e. for kpoints other than the '
         write (*,*) ' gamma point, the complex version of LAPACK is needed.'
         write (*,*) ' We must stop here! '
         stop
        end if
        a0 = 0.0d0
        a1 = 1.0d0

        allocate (xxxx(norbitals,norbitals))
        allocate (yyyy(norbitals,norbitals))
        allocate (zzzz(norbitals,norbitals))
!NPA
        if (iqout .eq. 3) then
         allocate (ssss(norbitals,norbitals))
         allocate (ww(norbitals))
        endif

        liwork = 15*norbitals
        allocate (iwork(liwork))
!        lwork = 100*norbitals + 3*norbitals*norbitals
!        allocate (work(lwork))
        lwork = 1
        allocate(work(lwork))

        if (.not. allocated(sm12_save)) then
         allocate (sm12_save(norbitals,norbitals))
        end if

! We built H/S(K) going through the list of atoms and their neighbors.
! How do we know where the particular (i,j) spot belongs in the grand
! matrix?  This should help you understand the degelec shifting business:
!
!                  atom 1         atom2            atom  3
!
!              nu1 nu2 nu3   nu1 nu2 nu3 nu4     nu1 nu2 nu3
!
!                           _________________
!         mu1               |               |
!  atom 1 mu2               |    H(1,2)     |
!         mu3               |               |
!                           -----------------
!         mu1
!  atom 2 mu2
!         mu3
!         mu4
!
!         mu1
!  atom3  mu2
!         mu3
!
! to the (1,2) portion at the right place we use degelec(iatom), which is
! passed, it remembers how many orbitals one must skip to get to the spot
! reserved for iatom, e.g. in this example degelec(1)=0, degelc(2)=3.

!
! COMPUTE S(K) AND H(K)
! ****************************************************************************
! Find the overlap and Hamiltonian matrices s(mu,nu,i,j) and h(mu,nu,i,j)
! in k space. Here iatom is an atom in the central cell and jatom a
! neighboring atom.

! Initialize to zero first
        zzzz = a0
        yyyy = a0
        xxxx = a0

!NPA
        if (iqout .eq. 3) then
         do inu = 1, norbitals
          imu = getissh(inu)
          iatom = getiatom(inu)
          lm = getlssh(inu)
          in1 = imass(iatom)
          if(Qneutral(getissh(inu),imass(getiatom(inu))).lt.0.01)then
            ww(inu) = 1.0d0
           else
            ww(inu) = 10.0d0
          endif
         enddo
        endif

! We now loop over all neighbors jatom of iatom.
        do iatom = 1, natoms
         in1 = imass(iatom)

! Now loop over all neighbors jatom of iatom.
         do ineigh = 1, neighn(iatom)
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom)

! So this matrix element goes in the i,j slot.
          do inu = 1, num_orb(in2)
           jnu = inu + degelec(jatom)
           do imu = 1, num_orb(in1)
            jmu = imu + degelec(iatom)
            zzzz(jmu,jnu) = zzzz(jmu,jnu) + s_mat(imu,inu,ineigh,iatom)
            yyyy(jmu,jnu) = yyyy(jmu,jnu) + h_mat(imu,inu,ineigh,iatom)
           end do ! do imu
          end do ! do inu
         end do ! do ineigh

! Now loop over all neighbors jatom of iatom VNL
         do ineigh = 1, neighPPn(iatom)
          mbeta = neighPP_b(ineigh,iatom)
          jatom = neighPP_j(ineigh,iatom)
          in2 = imass(jatom)

          do inu = 1, num_orb(in2)
           jnu = inu + degelec(jatom)
           do imu = 1, num_orb(in1)
            jmu = imu + degelec(iatom)
            yyyy(jmu,jnu) = yyyy(jmu,jnu) + vnl(imu,inu,ineigh,iatom)
           end do ! do imu
          end do ! do inu
         end do ! do inegh

        end do ! do iatom

! xxxx = unused
! zzzz = Overlap in AO basis
! yyyy = Hamiltonian in AO basis

!NPA
        if (iqout .eq. 3) then

         do inu = 1, norbitals
          do imu = 1, norbitals
           ssss(inu,imu) = zzzz(inu,imu)
          end do
         end do
         do inu = 1, norbitals
          do imu = 1, norbitals
           zzzz(inu,imu) = zzzz(inu,imu)*ww(inu)*ww(imu)
          end do
         end do

        endif  ! end if (iqout .eq. 3)


! DIAGONALIZE THE OVERLAP MATRIX
! ****************************************************************************
! If you are within the scf loop, you do not have to recalculate the overlap.
! NPA

! Call the diagonalizer
!         if (wrtout) then
!           write (*,*) ' Call diagonalizer for overlap. '
!           write (*,*) '                  The overlap eigenvalues: '
!           write (*,*) ' ******************************************************* '
!         end if
!           call dsyevd('V', 'U', norbitals, zzzz, norbitals, slam, work,    &
!     &                lwork, iwork, liwork, info )
           call dsyev ('V', 'U', norbitals, zzzz, norbitals, slam, work, -1, info)
! resize working space
           lwork = work(1)
           deallocate (work)
           allocate(work(lwork))
! diagonalize the overlap matrix with the new working space
           call dsyev ('V', 'U', norbitals, zzzz, norbitals, slam, work,     &
     &                lwork, info )
        
         if (info .ne. 0) call diag_error (info, 0)

! xxxx = unused
! zzzz = Overlap eigenvectors
! yyyy = Hamiltonian

! CHECK THE LINEAR DEPENDENCE
! ****************************************************************************
! Fix the linear dependence problem. References: Szabo and Ostlund, Modern
! Quantum Chem. McGraw Hill 1989 p. 142; Szabo and Ostlund, Modern Quantum
! Chem. Dover 1996 p. 145. A tolerance for a small overlap eigenvalue is
! set by overtol.

! Determine the smallest active eigenvector
         mineig = 0
         do imu = 1, norbitals
          if (slam(imu) .lt. overtol) mineig = imu
         end do

! You can specify a specific number of orbitals to drop with this
! next line, by uncommenting it.
! mineig = 0  {Don't drop any}
         mineig = mineig + 1
         norbitals_new = norbitals + 1 - mineig
         if (norbitals_new .ne. norbitals) then
          write (*,*) '  '
          write (*,*) ' WARNING. ############################ '
          write (*,*) ' Linear dependence encountered in basis set. '
          write (*,*) ' An overlap eigenvalue is very small. '
          write (*,*) norbitals - norbitals_new, ' vectors removed. '
          write (*,*) ' Spurious orbital energies near zero will '
          write (*,*) ' appear as a result of dropping these orbitals'
          write (*,*) ' You can change this by adjusting overtol in '
          write (*,*) ' kspace.f '
          write (*,*) '  '
           write(*,*) ' '
           write(*,*) ' Eigenvectors that correspond to eigenvalues'
           write(*,*) ' that were eliminated.  These might provide'
           write(*,*) ' insight into what atomic orbitals are causing'
           write(*,*) ' the problem.'
           write(*,*) ' '
           write (*,*) '            The overlap eigenvalues: '
           write (*,*) ' ********************************************** '
           write (*,200) (slam(imu), imu = 1, norbitals)
           write(*,*) ' '
           write(*,*) ' Eigenvectors that correspond to eigenvalues'
           write(*,*) ' that were eliminated.  These might provide'
           write(*,*) ' insight into what atomic orbitals are causing'
           write(*,*) ' the problem.'
           write(*,*) ' '
           do imu = 1, mineig - 1
            write(*,*) ' eigenvector',imu
            do jmu = 1, norbitals
             write(*,*) jmu,' ',zzzz(jmu,imu)
            end do
          end do
          write (*,*) ' '

          do imu = mineig, norbitals
           jmu = imu - mineig + 1
           zzzz(:,jmu) = zzzz(:,imu)
           slam(jmu) = slam(imu)
          end do
       end if
!
!! CALCULATE (S^-1/2) --> sm1
! ****************************************************************************
! In a diagonal reperesentation (Udagger*S*U = s, s is a diagonal matrix)
!! We just take the inverse of the square roots of the eigenvalues to get
! s^-1/2. Then we 'undiagonalize' the s^-1/2 matrix back to get
! S^-1/2 = U*s^-1/2*Udagger.
! Note: We do S^-1/4 here, because the sqlami contribution get squared
! after it is combined with overlap.
         do imu = 1, norbitals_new
          sqlami = slam(imu)**(-0.25d0)
          zzzz(:,imu) = zzzz(:,imu)*sqlami
         end do
         call dgemm ('N', 'C', norbitals, norbitals, norbitals_new, a1, zzzz,&
     &               norbitals, zzzz, norbitals, a0, xxxx, norbitals)

!NPA  now we do X=W(WSW)^-1/2, before X=S^-1/2
         if (iqout .eq. 3) then
          do imu=1, norbitals
           xxxx(imu,:)=xxxx(imu,:)*ww(imu)
          end do
         endif


! Now put S^-1/2 into s(k)^-1/2, this will be remembered for the duration of
! the scf cycle.
         do inu = 1, norbitals
          do imu = 1, norbitals
           sm12_save(imu,inu) = xxxx(imu,inu)
          end do
         end do


! Now if not first iteration
! Restore S^-1/2 from s(k)^-1/2,
         xxxx(:,:) = sm12_save(:,:)
! xxxx = S^-1/2 in AO basis
! zzzz = Unused (used as temporary work area below)
! yyyy = Hamiltonian in AO basis
!
! CALCULATE (S^-1/2)*H*(S^-1/2)
! ****************************************************************************
        if (iqout .ne. 3) then
! Set M=H*(S^-.5)
         call dsymm ( 'R', 'U', norbitals, norbitals, a1, xxxx, &
     &               norbitals, yyyy, norbitals, a0, zzzz, norbitals )

! Set Z=(S^-.5)*M
         call dsymm ( 'L', 'U', norbitals, norbitals, a1, xxxx, &
     &               norbitals, zzzz, norbitals, a0, yyyy, norbitals )

        else
! Set conjg((W(WSW)^-1/2)T)*H
         call dgemm ( 'C', 'N', norbitals, norbitals, norbitals, a1, xxxx, &
     &               norbitals, yyyy, norbitals, a0, zzzz, norbitals )
! Set M*(W(WSW)^-1/2)
         call dgemm ( 'N', 'N', norbitals, norbitals, norbitals, a1, zzzz, &
     &               norbitals, xxxx, norbitals, a0, yyyy, norbitals )
! so we have conjg((W(WSW)^-1/2)T)*H*(W(WSW)^-1/2) now
        endif
! xxxx = S^-1/2 in AO basis
! zzzz = Unused
! yyyy = Hamiltonian in the MO basis set

! DIAGONALIZE THE HAMILTONIAN.
! ****************************************************************************
        if (wrtout) then
          write (*,*) '  '
          write (*,*) ' Call diagonalizer for Hamiltonian. '
          write (*,*) '            The energy eigenvalues: '
          write (*,*) ' *********************************************** '
        end if
! Eigenvectors are needed to calculate the charges and for forces!

! reset working space
          lwork = 1
          deallocate (work)
          allocate (work(lwork))
! first find optimal working space
          call dsyev ('V', 'U', norbitals, yyyy, norbitals, eigen, work, -1, info)
! resize working space
          lwork = work(1)
          deallocate (work)
          allocate(work(lwork))
          call dsyev ('V', 'U', norbitals, yyyy, norbitals, eigen, work,     &
     &               lwork, info )

        if (info .ne. 0) call diag_error (info, 0)

!
! INFORMATION FOR THE LOWDIN CHARGES
! ****************************************************************************
! xxxx = S^-1/2 in AO basis
! zzzz = Unused
! yyyy = Hamiltonian eigenvectors in the MO basis
        if (iqout .ne. 2) blowre(:,:,ikpoint) = real(yyyy(:,:))
        if (iqout .ne. 2 .and. icluster .ne. 1) blowim(:,:,ikpoint) = 0

        call dsymm ( 'L', 'U', norbitals, norbitals, a1, xxxx, &
     &               norbitals, yyyy, norbitals, a0, zzzz, norbitals )

        if (iqout .ne. 3) then
         call dsymm ( 'L', 'U', norbitals, norbitals, a1, xxxx, &
     &               norbitals, yyyy, norbitals, a0, zzzz, norbitals )
!NPA
        else
         call dgemm ( 'N', 'N', norbitals, norbitals, norbitals, a1, xxxx,   &
     &               norbitals, yyyy, norbitals, a0, zzzz, norbitals )
        end if

        bbnkre(:,:,ikpoint) = real(zzzz(:,:))
        if (icluster .ne. 1) bbnkim(:,:,ikpoint) = 0

! Deallocate Arrays'
! ===========================================================================
        deallocate (xxxx)
        deallocate (yyyy)
        deallocate (zzzz)

! NPA
        if (iqout .eq. 3) then
              deallocate (ww)
             deallocate (ssss)
        endif

        deallocate (work)
        deallocate (iwork)

! Format Statements
! ===========================================================================
100     format (2x, ' eigenvalue(1) = ', f10.6, &
     &              ' eigenvalue(norbitals) = ', f10.6)
200     format (4(2x, f12.4))

        return
      end subroutine calc_blowreim







subroutine write_dipole (estring)
        use charges
        use configuration
        use constants_fireball
        use density
        use dimensions
        use interactions
        use neighbor_map
        use kpoints
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Argument Declaration and Description
! ===========================================================================
! Input

        character (len = 25), intent(in) :: estring
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer imu, inu
        integer ineigh
        integer in1, in2
        integer issh, jssh
        integer jatom
        integer jneigh
        integer mqn
        integer mbeta
        integer mmu
        integer noccupy
        integer nnu
        real    Qtot, Qtot1, Qtot2
        real    dip_x, dipQout_x, dipTot_x, dipProy_x, dipIntra_x
        real    dip_y, dipQout_y, dipTot_y, dipProy_y, dipIntra_y
        real    dip_z, dipQout_z, dipTot_z, dipProy_z, dipIntra_z
        real    dip_res_x, dip_res_y, dip_res_z, dip_res_tot
        real    dip_tot, dip_proy, dipQout_tot, dipProy_tot, dipTot_tot, dipIntra_tot
        real, dimension(3) :: Rbc,u21


        real, dimension (3) :: vec, r1, r2, r21
        real, parameter ::  Debye = 0.20819

! ****************************************************************************
! W R I T E   T H E   D I P O L E   (For testing purposes)
! ****************************************************************************
      allocate (Q0_TOT(natoms))
      do iatom = 1, natoms
         Q0_TOT(iatom) = 0
         in1 = imass(iatom)
         do issh = 1, nssh(in1)
            Q0_TOT(iatom) = Q0_TOT(iatom) + Qneutral(issh,in1)
         end do
        end do

      dip_x=0.0d0
      dip_y=0.0d0
      dip_z=0.0d0

      do iatom = 1, natoms
         Qtot=-Q0_TOT(iatom)
         in1 = imass(iatom)
         do imu = 1,num_orb(in1)
             Qtot = Qtot+Qout(imu,iatom)
         end do
            
            !write(*,*) 'ratom is ',ratom(1,iatom),ratom(2,iatom),ratom(3,iatom)
            dip_x = dip_x+Qtot*ratom(1,iatom)
            dip_y = dip_y+Qtot*ratom(2,iatom)
            dip_z = dip_z+Qtot*ratom(3,iatom)

       enddo !end do iatom = 1,natoms


        dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )   

         dipQout_x=dip_x
         dipQout_y=dip_y
         dipQout_z=dip_z
         dipQout_tot=dip_tot

     open( unit = 1729, file = 'dipole_Qout'//trim(estring), status = 'unknown',  &
                                   &     position = 'append')
       
       write(1729,446) -dip_x/Debye, -dip_y/Debye, -dip_z/Debye, dip_tot/Debye

     close(1729)

      dip_x=0.0d0
      dip_y=0.0d0
      dip_z=0.0d0

      do iatom = 1, natoms
         in1 = imass(iatom)
         r1(:) = ratom(:,iatom)
         dip_x = dip_x-Q0_TOT(iatom)*r1(1)
         dip_y = dip_y-Q0_TOT(iatom)*r1(2)         
         dip_z = dip_z-Q0_TOT(iatom)*r1(3)
      end do !end do iatom = 1,natoms

      do iatom = 1, natoms
         in1 = imass(iatom)
         r1(:) = ratom(:,iatom)         
         do ineigh = 1,neighn(iatom)
            mbeta = neigh_b(ineigh,iatom)
            jatom = neigh_j(ineigh,iatom)
            r2(:) = ratom(:,jatom) + xl(:,mbeta)
            in2 = imass(jatom)

            Rbc(:)=(r1(:)+r2(:))/2.0d0
            u21(:)=(r2(:)-r1(:))/(sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
            do imu = 1,num_orb(in1)
               do inu = 1,num_orb(in2)
                 
                !if (iatom .ne. jatom) then !.or. imu .eq. inu) then
                !else 
                  dip_x =dip_x+rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom)+ &
                         & Rbc(1)*s_mat(imu,inu,ineigh,iatom))
                  dip_y =dip_y+rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom)+ &
                         & Rbc(2)*s_mat(imu,inu,ineigh,iatom))
                  dip_z=dip_z+rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom)+ &
                         & Rbc(3)*s_mat(imu,inu,ineigh,iatom))


                 !dip_proy = dipc(1,imu,inu,ineigh,iatom)*u21(1)+ &
                 !           & dipc(2,imu,inu,ineigh,iatom)*u21(2)+ &
                 !           & dipc(3,imu,inu,ineigh,iatom)*u21(3)

                  !dip_x =  dip_x+rho(imu,inu,ineigh,iatom)*(dip_proy*u21(1)+ &
                  !       & Rbc(1)*s_mat(imu,inu,ineigh,iatom))
                  !dip_y =  dip_y+rho(imu,inu,ineigh,iatom)*(dip_proy*u21(2)+ &
                  !       & Rbc(2)*s_mat(imu,inu,ineigh,iatom))
                  !dip_z =  dip_z+rho(imu,inu,ineigh,iatom)*(dip_proy*u21(3)+ &
                  !       & Rbc(3)*s_mat(imu,inu,ineigh,iatom))

                !end if !end if
               end do !end do inu
            end do !end do imu
         end do !end ineigh = 1,natoms
      end do ! end do iatom = 1,natoms
     dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )

         dipTot_x=dip_x
         dipTot_y=dip_y
         dipTot_z=dip_z
         dipTot_tot=dip_tot
    
     open( unit = 666, file = 'dipole_Tot'//trim(estring), status = 'unknown',  &
                                   &     position = 'append')

     write(666,446) -dip_x/Debye, -dip_y/Debye, -dip_z/Debye, dip_tot/Debye
     close(666)
    !W R I T E     T H E     B I G     F I L E

               !write the intraatomic dipole
!DIP PROY

      dip_x=0.0d0
      dip_y=0.0d0
      dip_z=0.0d0
      
      do iatom = 1, natoms
         in1 = imass(iatom)
         r1(:) = ratom(:,iatom)
         dip_x = dip_x-Q0_TOT(iatom)*r1(1)
         dip_y = dip_y-Q0_TOT(iatom)*r1(2)         
         dip_z = dip_z-Q0_TOT(iatom)*r1(3)
      end do !end do iatom = 1,natoms


      do iatom = 1, natoms
         in1 = imass(iatom)
         r1(:) = ratom(:,iatom)         
         do ineigh = 1,neighn(iatom)
            mbeta = neigh_b(ineigh,iatom)
            jatom = neigh_j(ineigh,iatom)
            r2(:) = ratom(:,jatom) + xl(:,mbeta)
            in2 = imass(jatom)

            Rbc(:)=(r1(:)+r2(:))/2.0d0
            do imu = 1,num_orb(in1)
               do inu = 1,num_orb(in2)
                  if (iatom .ne. jatom) then !.or. imu .eq. inu) then
                   u21(:)=(r2(:)-r1(:))/(sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
                 

                    dip_proy = dipc(1,imu,inu,ineigh,iatom)*u21(1)+ &
                            & dipc(2,imu,inu,ineigh,iatom)*u21(2)+ &
                            & dipc(3,imu,inu,ineigh,iatom)*u21(3)

                    dip_x =  dip_x+rho(imu,inu,ineigh,iatom)*(dip_proy*u21(1)+ &
                         & Rbc(1)*s_mat(imu,inu,ineigh,iatom))
                    dip_y =  dip_y+rho(imu,inu,ineigh,iatom)*(dip_proy*u21(2)+ &
                         & Rbc(2)*s_mat(imu,inu,ineigh,iatom))
                    dip_z =  dip_z+rho(imu,inu,ineigh,iatom)*(dip_proy*u21(3)+ &
                         & Rbc(3)*s_mat(imu,inu,ineigh,iatom))
                  else
                    dip_x =dip_x+rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom)+ &
                         & Rbc(1)*s_mat(imu,inu,ineigh,iatom))
                    dip_y =dip_y+rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom)+ &
                         & Rbc(2)*s_mat(imu,inu,ineigh,iatom))
                    dip_z=dip_z+rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom)+ &
                         & Rbc(3)*s_mat(imu,inu,ineigh,iatom))
                end if !end if
               end do !end do inu
            end do !end do imu
         end do !end ineigh = 1,natoms
      end do ! end do iatom = 1,natoms
     dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )

         dipProy_x=dip_x
         dipProy_y=dip_y
         dipProy_z=dip_z
         dipProy_tot=dip_tot
    
     open( unit = 166, file = 'dipole_Tot_proy'//trim(estring), status = 'unknown',  &
                                   &     position = 'append')

     write(166,446) -dip_x/Debye, -dip_y/Debye, -dip_z/Debye, dip_tot/Debye
     close(166)

     !W R I T E     T H E     B I G     F I L E

               !write the intraatomic dipole

      dip_x=0.0d0
      dip_y=0.0d0
      dip_z=0.0d0


      do iatom = 1, natoms
         in1 = imass(iatom)
         r1(:) = ratom(:,iatom)
         do ineigh = 1,neighn(iatom)
            mbeta = neigh_b(ineigh,iatom)
            jatom = neigh_j(ineigh,iatom)
            r2(:) = ratom(:,jatom) + xl(:,mbeta)
            in2 = imass(jatom)

            Rbc(:)=(r1(:)+r2(:))/2.0d0
            u21(:)=(r2(:)-r1(:))/(sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
            do imu = 1,num_orb(in1)
               do inu = 1,num_orb(in2)

                if ((iatom .eq. jatom) .and. (imu .ne. inu)) then
                !else 
                  dip_x=dip_x+rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom)+ &
                         & Rbc(1)*s_mat(imu,inu,ineigh,iatom))
                  dip_y=dip_y+rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom)+ &
                         & Rbc(2)*s_mat(imu,inu,ineigh,iatom))
                  dip_z=dip_z+rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom)+ &
                         & Rbc(3)*s_mat(imu,inu,ineigh,iatom))

                end if !end if
               end do !end do inu
            end do !end do imu
         end do !end ineigh = 1,natoms
      end do ! end do iatom = 1,natoms
           dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )

         dipIntra_x=dip_x
         dipIntra_y=dip_y
         dipIntra_z=dip_z
         dipIntra_tot=dip_tot
     



         open( unit = 667, file = 'Charges_and_Dipoles'//trim(estring), status ='unknown',  &
                                    &     position = 'append')
         write(667,*)   '+++++++++++++++++++ NEW STEP++++++++++++++++++' 
         write(667,444) 'dip_TOT',-dipTot_x/Debye, -dipTot_y/Debye,-dipTot_z/Debye,dipTot_tot/Debye
         write(667,444) 'dip_Pro',-dipProy_x/Debye, -dipProy_y/Debye,-dipProy_z/Debye,dipProy_tot/Debye
         write(667,444) 'dip_Qout',-dipQout_x/Debye, -dipQout_y/Debye, -dipQout_z/Debye, dipQout_tot/Debye
         write(667,444) 'dip_Int',-dipIntra_x/Debye, -dipIntra_y/Debye,-dipIntra_z/Debye,dipIntra_tot/Debye
         !write(667,*)Qout
!DIP INTRA


      do iatom = 1, natoms
         dip_x=0.0d0
         dip_y=0.0d0
         dip_z=0.0d0
         in1 = imass(iatom)
         r1(:) = ratom(:,iatom)
         do ineigh = 1,neighn(iatom)
            mbeta = neigh_b(ineigh,iatom)
            jatom = neigh_j(ineigh,iatom)
            r2(:) = ratom(:,jatom) + xl(:,mbeta)
            in2 = imass(jatom)

            Rbc(:)=(r1(:)+r2(:))/2.0d0
            u21(:)=(r2(:)-r1(:))/(sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))
            do imu = 1,num_orb(in1)
               do inu = 1,num_orb(in2)

                if ((iatom .eq. jatom) .and. (imu .ne. inu)) then
                !else 
                  dip_x=dip_x+rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom)+ &
                         & Rbc(1)*s_mat(imu,inu,ineigh,iatom))
                  dip_y=dip_y+rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom)+ &
                         & Rbc(2)*s_mat(imu,inu,ineigh,iatom))
                  dip_z=dip_z+rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom)+ &
                         & Rbc(3)*s_mat(imu,inu,ineigh,iatom))

                end if !end if
               end do !end do inu
            end do !end do imu
         end do !end ineigh = 1,natoms
      dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )
      write(667,455) symbol(iatom),'dip_int',-(dip_x/Debye),-(dip_y/Debye),-(dip_z/Debye),(dip_tot/Debye)
      end do ! end do iatom = 1,natoms


!DIP RES = DIP_TOT - PROY
      dip_x=0.0d0
      dip_y=0.0d0
      dip_z=0.0d0
      dip_res_x = 0.0d0
      dip_res_y = 0.0d0
      dip_res_z = 0.0d0
      
      do iatom = 1, natoms
         dip_res_x = 0.0d0
         dip_res_y = 0.0d0
         dip_res_z = 0.0d0
         in1 = imass(iatom)
         r1(:) = ratom(:,iatom)         
         do ineigh = 1,neighn(iatom)
            mbeta = neigh_b(ineigh,iatom)
            jatom = neigh_j(ineigh,iatom)
            r2(:) = ratom(:,jatom) + xl(:,mbeta)
            in2 = imass(jatom)

            Rbc(:)=(r1(:)+r2(:))/2.0d0
            do imu = 1,num_orb(in1)
               do inu = 1,num_orb(in2)
                  if (iatom .ne. jatom) then !.or. imu .eq. inu) then
                   u21(:)=(r2(:)-r1(:))/(sqrt((r2(1)-r1(1))**2+(r2(2)-r1(2))**2+(r2(3)-r1(3))**2))

                    dip_proy = dipc(1,imu,inu,ineigh,iatom)*u21(1)+ &
                            & dipc(2,imu,inu,ineigh,iatom)*u21(2)+ &
                            & dipc(3,imu,inu,ineigh,iatom)*u21(3)

                    dip_res_x = dip_res_x + rho(imu,inu,ineigh,iatom)*(dipc(1,imu,inu,ineigh,iatom)- &
                         & dip_proy*u21(1))

                    dip_res_y = dip_res_y + rho(imu,inu,ineigh,iatom)*(dipc(2,imu,inu,ineigh,iatom)- &
                         & dip_proy*u21(2))

                    dip_res_z = dip_res_z + rho(imu,inu,ineigh,iatom)*(dipc(3,imu,inu,ineigh,iatom)- &
                         & dip_proy*u21(3))

                end if !end if
               end do !end do inu
            end do !end do imu
         end do !end ineigh = 1,natoms
         dip_res_tot = sqrt (dip_res_x**2 + dip_res_y**2 + dip_res_z**2 )
         dip_x=dip_x+dip_res_x
         dip_y=dip_y+dip_res_y
         dip_z=dip_z+dip_res_z
         write(667,455) symbol(iatom),'dip_res',-(dip_res_x/Debye),-(dip_res_y/Debye),-(dip_res_z/Debye),(dip_res_tot/Debye)
      end do ! end do iatom = 1,natoms
      dip_tot = sqrt (dip_x**2 + dip_y**2 + dip_z**2 )
      write(667,444)'dip_res_TOT',-(dip_x/Debye),-(dip_y/Debye),-(dip_z/Debye),(dip_tot/Debye) 
     
      deallocate(Q0_TOT)
444     format (a7,4f10.4)
446     format (4f10.4)
445     format (a2,4f10.4)
455     format (a2,a7,4f10.4)
        return
      end subroutine write_dipole


      subroutine write_charges_shell(estring)
        use charges
        use density
        use interactions
        use configuration
        implicit none
        real, dimension(natoms):: Ntot
        real Ntot_all
        integer iatom,in1,imu
        character (len = 25), intent(in) :: estring
        Ntot_all=0.0
        do iatom = 1, natoms
          in1 = imass(iatom)
          Ntot(iatom) = 0.0d0
          do imu = 1, num_orb(in1)
            Ntot(iatom) = Ntot(iatom) + Qout(imu,iatom)
          end do
          Ntot_all = Ntot_all + Ntot(iatom)
        end do
        open(unit = 101, file = estring, status = 'unknown')
          write(101,'(2x, i4,2x,f10.6)') natoms,Ntot_all
          do iatom = 1, natoms
            in1 = imass(iatom)
            write(101,'(2x, i4,2x,f10.6)',advance='no') iatom, Ntot(iatom)
            do imu = 1, nssh(in1)
              write(101,'(f10.6)',advance='no') Qout(imu,iatom)
            end do
            write(101,*)
          end do ! do iatom
        close(101)
      end subroutine write_charges_shell
 


