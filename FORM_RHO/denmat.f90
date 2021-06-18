! copyright info:
!
!                             @@Copyright 2001
!                           Fireball Committee
! Brigham Young University of Utah - James P. Lewis, Chair
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
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine denmat (ifixcharge, iqout, icluster, iwrtefermi, tempfe,  &
        &                  ebs, iwrtpop, bmix, Kscf, igap)
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
        integer, intent (in) :: iwrtpop
        integer, intent (in) :: bmix
        integer, intent (in) :: Kscf
        integer, intent (in) :: igap

        real, intent (in) :: tempfe

! Output
        real, intent (out) :: ebs

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iband
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
        integer issh1, mu_min, mu_max, l, inumorb
        integer :: info, lwork
        !integer, dimension(nssh_tot) :: ipiv
        integer, dimension(100) :: work
        integer :: beta, alpha, ialp, ina, matom
        !real,dimension(nssh_tot,nssh_tot) :: M
        !real,dimension(1,nssh_tot) :: B
        !real,dimension(nssh_tot+1,nssh_tot+1) :: M
        !real,dimension(1,nssh_tot+1) :: B
        !integer, dimension(:), allocatable :: ipiv
        !real,dimension(:,:),allocatable :: M
        !real,dimension(:,:), allocatable :: B
        real auxgS
        real Ntot

        integer, dimension (norbitals) :: ioccupy
        integer, dimension (norbitals, nkpoints) :: ioccupy_k

        real aux1, aux2, aux3
        real deltae
        real dot
        real gutr
        real pcharge
        real ztest
        real checksum
        real Wmu
        real y
        real, dimension (norbitals, nkpoints) :: foccupy
        real, dimension (numorb_max, natoms) :: QMulliken
        real, dimension (3) :: vec, r1, r2, r21
        real, dimension (natoms) :: pqmu

        complex ai
        complex phase, phasex
        complex step1, step2

        logical read_occupy
        real,dimension(nssh_tot,nssh_tot) :: A
        real,dimension(nssh_tot) :: c, SQ ! carga
        real,dimension(nssh_tot) :: LB, UB, nalpha
        real :: diff_err,Ep2

! Procedure
! ===========================================================================
! Initialize some things
        ai = cmplx(0.0d0,1.0d0)
! jel: where you set rho to zero ??? JOM: in build_rho
        rhoPP = 0.0d0

!        write (*,*) '  '
!        write (*,*) '  '
!        write (*,*) ' ****************************************************** '
!        write (*,*) '  '
!        write (*,*) '                   Welcome to denmat --              '
!        write (*,*) '  '
!        write (*,*) ' ****************************************************** '
!
! ****************************************************************************
!
!             if (iqout .eq. 6)then
!               allocate( M(nssh_tot+1,nssh_tot+1) )
!               print*,'allocate M',nssh_tot+1,nssh_tot+1,'=',nssh_tot+1*nssh_tot+1
!               allocate( B(1,nssh_tot+1) )
!               allocate( ipiv(nssh_tot) )
!             else
!               allocate( M(1,1) )
!               allocate( B(1,1) )
!               allocate( ipiv(1) )
!             end if
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
          write (*,*) ' ***************************************** '

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
           !write (*,100) iband, iatom, pcharge
          end do
         end do
        end if

! ****************************************************************************
! Get the Fermi energy.
        call fermie (norbitals, ztot, eigen_k, efermi, ioccupy_k, foccupy)
!        write (*,*) ' Fermi Level = ', efermi

        if (iwrtefermi .eq. 1) then
         write (*,*) '  '
         write (*,*) ' We write out the occupancies of the levels  '
         write (*,*) ' ----- ioccupy_k, from the subroutine fermie '
         do ikpoint = 1, nkpoints
          write (*,*) '  '
          write (*,*) ' ------ fermi ioccupy_k for k-point = ', ikpoint
          do iband = 1, norbitals_new
           write (*,200) iband, ioccupy_k(iband,ikpoint)
          end do
         end do

         write (*,*) '  '
         write (*,*) ' We write out the electron fraction in each level '
         write (*,*) ' ------ from the subroutine fermie. '
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

! ****************************************************************************
!
!  C O M P U T E    L O W D I N    C H A R G E S
! ****************************************************************************
! Initialize
        if (iqout .eq. 1 .or. iqout .eq. 3) then
         Qout = 0.0d0
         QLowdin_TOT = 0.0d0

         if (ifixcharge .eq. 1) then

          do iatom = 1, natoms
           in1 = imass(iatom)
           do issh = 1, nssh(in1)
            Qout(issh,iatom) = Qin(issh,iatom)
            QLowdin_TOT(iatom) = QLowdin_TOT(iatom) + Qin(issh,iatom)
           end do
          end do

         else

          do iatom = 1, natoms
           in1 = imass(iatom)

! Loop over the special k points.
           do ikpoint = 1, nkpoints
            aux1 = weight_k(ikpoint)*spin
            do iorbital = 1, norbitals
             if (ioccupy_k(iorbital,ikpoint) .eq. 1) then
              aux2 = aux1*foccupy(iorbital,ikpoint)

! Finally the imu loop.
              imu = 0
              do issh = 1, nssh(in1)
               do mqn = 1, 2*lssh(issh,in1) + 1
                imu = imu + 1
                mmu = imu + degelec(iatom)
                if (icluster .ne. 1) then
                 aux3 = aux2*(blowre(mmu,iorbital,ikpoint)**2                &
     &                        + blowim(mmu,iorbital,ikpoint)**2)
                else
                 aux3 = aux2*blowre(mmu,iorbital,ikpoint)**2
                end if
                Qout(issh,iatom) = Qout(issh,iatom) + aux3
                QLowdin_TOT(iatom) = QLowdin_TOT(iatom) + aux3
               end do
              end do
             end if

! End loop over orbitals and kpoints
            end do
           end do

! End loop over atoms
          end do

         end if      ! endif of ifixcharges
        end if       ! endif of iqout .ne. 1


! ****************************************************************************
!
!  C O M P U T E    M U L L I K E N    C H A R G E S
! ****************************************************************************
! Compute Mulliken charges.
        if (iqout .eq. 2) then

         Qout = 0.0d0
         QMulliken = 0.0d0
         QMulliken_TOT = 0.0d0

         if (ifixcharge .eq. 1) then

          do iatom = 1, natoms
           in1 = imass(iatom)
           QMulliken_TOT(iatom) = 0.0d0
           do issh = 1, nssh(in1)
            Qout(issh,iatom) = Qin(issh,iatom)
            QMulliken_TOT(iatom) = QMulliken_TOT(iatom) + Qin(issh,iatom)
           end do
          end do

         else

          do iatom = 1, natoms
           in1 = imass(iatom)

! Loop over neighbors
           do ineigh = 1, neighn(iatom)
            jatom = neigh_j(ineigh,iatom)
            in2 = imass(jatom)

     
            jneigh = neigh_back(iatom,ineigh)
            do imu = 1, num_orb(in1)
             do inu = 1, num_orb(in2)
              QMulliken(imu,iatom) = QMulliken(imu,iatom)                    &
     &        + 0.5d0*(rho(imu,inu,ineigh,iatom)*s_mat(imu,inu,ineigh,iatom) &
     &        + rho(inu,imu,jneigh,jatom)*s_mat(inu,imu,jneigh,jatom))
             end do
            end do

! End loop over neighbors
           end do

! Finally the imu loop.
           imu = 0
           do issh = 1, nssh(in1)
            do mqn = 1, 2*lssh(issh,in1) + 1
               imu = imu + 1
               Qout(issh,iatom) = Qout(issh,iatom) + QMulliken(imu,iatom)
            end do
               QMulliken_TOT(iatom) = QMulliken_TOT(iatom) + Qout(issh,iatom)
           end do

! End loop over atoms
          end do
         end if     ! endif of ifixcharges
        end if      ! endif of iqout .eq. 2


! GAP ENRIQUE-FF
       if ((igap.eq.1).or.(igap.eq.2)) then
         call buildnij(ioccupy_k,foccupy, bmix, Kscf)
       end if
       if (igap .eq. 2 ) then
         call koopman(natoms,nkpoints,ratom,ioccupy_k,foccupy)
       end if
! end GAP ENRIQUE-FF





! ****************************************************************************
!
!  C O M P U T E    M U L L I K E N - D I P O L E    C H A R G E S
! ****************************************************************************
! Compute Mulliken-dipole charges.
        if (iqout .eq. 4) then
         Qout = 0.0d0
         QMulliken = 0.0d0
         QMulliken_TOT = 0.0d0

         if (ifixcharge .eq. 1) then

          do iatom = 1, natoms
           in1 = imass(iatom)
           QMulliken_TOT(iatom) = 0.0d0
           do issh = 1, nssh(in1)
            Qout(issh,iatom) = Qin(issh,iatom)
            QMulliken_TOT(iatom) = QMulliken_TOT(iatom) +Qin(issh,iatom)
           end do
          end do

         else

          do iatom = 1, natoms
           in1 = imass(iatom)
           r1(:) = ratom(:,iatom)
! Loop over neighbors
           do ineigh = 1, neighn(iatom)
            jatom = neigh_j(ineigh,iatom)
            in2 = imass(jatom)
            r2(:) = ratom(:,jatom)


! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the
! bondcharge, and the bc distance, y
           r21(:) = r2(:) - r1(:)
           y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))

            jneigh = neigh_back(iatom,ineigh)
            do imu = 1, num_orb(in1)
             do inu = 1, num_orb(in2)
              QMulliken(imu,iatom) = QMulliken(imu,iatom)                   &
     &        +0.5d0*(rho(imu,inu,ineigh,iatom)*s_mat(imu,inu,ineigh,iatom) &
     &        + rho(inu,imu,jneigh,jatom)*s_mat(inu,imu,jneigh,jatom))
             end do
            end do
             
! dipole correction. Only if the two atoms are different
          if (y .gt. 1.0d-05) then
           

            do imu = 1, num_orb(in1)
             do inu = 1, num_orb(in2)

          !    write(*,*) 'Qb(',imu,',',iatom,')= ',QMulliken(imu,iatom)

              QMulliken(imu,iatom) = QMulliken(imu,iatom)+                  &
     &        (-rho(imu,inu,ineigh,iatom)*dip(imu,inu,ineigh,iatom)         &
     &        + rho(inu,imu,jneigh,jatom)*dip(inu,imu,jneigh,jatom))/y
             
 
           !         write(*,*) 'DIPOLE when iatom,jatom,imu,inu = ',  &
     ! &        iatom,jatom,imu,inu,' is',dip(imu,inu,ineigh,iatom),    &
     ! &        dip(inu,imu,jneigh,jatom) 

     !          write(*,*) 'Qa(',imu,',',iatom,')= ',QMulliken(imu,iatom)

             end do
            end do
          end if !end if y .gt. 1.0d-05)

         

! End loop over neighbors
           end do

! Finally the imu loop.
           imu = 0
           do issh = 1, nssh(in1)
            do mqn = 1, 2*lssh(issh,in1) + 1
               imu = imu + 1
               Qout(issh,iatom) = Qout(issh,iatom) + QMulliken(imu,iatom)
            end do
               QMulliken_TOT(iatom) = QMulliken_TOT(iatom) +Qout(issh,iatom)
           end do

!Check whether there are negative charges and correct
!If there's more than one shell whose charge is negative, more work is
!needed, but that'd be quite pathological a situation...
            do issh = 1, nssh(in1)

               if( Qout(issh,iatom) .lt. 0 .and. nssh(in1) .gt. 1 ) then
           
                  do jssh = 1,nssh(in1)

                     if ( jssh .ne. issh ) then

                        Qout(jssh,iatom) = Qout(jssh,iatom)+            &
                 &                         Qout(issh,iatom)/(nssh(in1)-1)

                     end if !end if jssh .ne. issh 

                  end do !end if jssh = 1,nssh(in1)

                  Qout(issh,iatom) = 0.0d0               

               end if !end if  Qout(issh,iatom) .lt. 0

            end do !end do issh = 1, nssh(in1)


! End loop over atoms
          end do
         end if     ! endif of ifixcharges
        end if      ! endif of iqout .eq. 4


! ***************************************************************************
!
! C O M P U T E   L2   DIAGONAL   CHARGES
!
!****************************************************************************

!      if (iqout .eq. 5) then
!
!         call charges_L2(Kscf,0)
!
!      endif

!****************************************************************************

! ****************************************************************************
!
! CHARGES VARIATIONAL
! ****************************************************************************
! integer :: beta, alpha, ialp, ina
! real,dimension(nssh_tot,nssh_tot) :: M
! real,dimension(nssh_tot) :: B
! In MODULES/interactions:
! real, dimension(:,:,:,:,:,:), allocatable :: gvhxc
! real, dimension(:,:,:,:), allocatable :: gvhxcS
! In ALLOCATIONS/allocate_rho:
! allocate (gvhxc(numorb_max,numorb_max,nssh_max,natoms,neigh_max,natoms))
! allocate (gvhxcS(nssh_max,nssh_max,natoms,natoms))
! (remember to include them also in deallocate_rho!!!)
      !B = 0.0d0
      SQ = 0.0d0
      c = 0.0d0
      if (iqout .eq. 6) then
          alpha = 0
          do ialp = 1, natoms
              ina = imass(ialp)
              !gvhxc
              do issh = 1, nssh(ina)
                 alpha = alpha + 1 ! transform to one index
                 !write(*,*) 'alpha indices', ialp, issh, alpha
                  beta = 0
                  do iatom = 1, natoms
                      in1 = imass(iatom)
                      matom = neigh_self(iatom)                    
                      inumorb = 1 ! counter for number of orbitals in atom iatom
                      do issh1 = 1, nssh(in1)
                          beta = beta + 1 ! transform to one index
                          !write(*,*) 'beta indices', iatom, issh1, beta
                          ! Spherical approximation to matrix elements:
                           l = lssh(issh1,in1)
                           auxgS = 0.0d0
                           ! define mu_min and mu_max: the orbitals
                           ! indices associated to the shell issh1
                           mu_min = inumorb
                           mu_max = mu_min+2*l
                           do imu = mu_min, mu_max 
                                auxgS =  auxgS &
                              & +  gvhxc(imu,imu,issh,ialp,matom,iatom)
                           end do ! end do imu = mu_min, mu_max 
                          auxgS = auxgS/(2*l+1)  ! 4*pi??
                          ! Now:
                          !M(alpha,beta) =  auxgS !gvhxcs(issh1,issh,iatom,ialp)                   
                          A(alpha,beta) = auxgS 
                          !write(*,*) alpha, beta, M(alpha,beta) 
                          inumorb = inumorb + 2*l+1
                      end do ! end do issh1
                      do ineigh = 1, neighn(iatom)
                          mbeta = neigh_b(ineigh,iatom)
                          jatom = neigh_j(ineigh,iatom)
                          in2 = imass(jatom) 
                          do imu = 1, num_orb(in1)
                             do inu = 1, num_orb(in2)
                                  !B(1,alpha) = B(1,alpha) + &
                                  c(alpha) = c(alpha) + &
                                  rho(imu,inu,ineigh,iatom)*gvhxc(imu, &
                                      &    inu,issh,ialp,ineigh,iatom)
                              end do ! end do inu
                         end do ! end do imu
                      end do ! end do ineigh
                  end do ! end do iatom
              end do ! end do issh
          end do ! end do ialp
          !M(nssh_tot+1,nssh_tot+1) = 0
          !B(1,nssh_tot+1) = ztot
          !do alpha = 1,nssh_tot
          ! M(nssh_tot+1,alpha) = 1
          ! M(alpha,nssh_tot+1) = 1
          !end do
              
                 !do beta = 1, nssh_tot
                 !write(*,*) 'alpha B ', beta, B(1,beta)
                 !end do
!  +++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  SOLVE SYSTEM Mx = B.  x are the charges
! 
      !write(*,*) gvhxc(1,1,1,1,neigh_self(1),1)
      !write(*,*) M

      ! do alpha = 1,nssh_tot
      !   do beta = alpha+1,nssh_tot
      !     M(alpha,beta) = 0.5*(M(alpha,beta)+M(beta,alpha))
      !     M(beta,alpha) = M(alpha,beta)
      !   end do !beta
      ! end do !alpha 

         !LWMAX = 100
         !call ssysv( 'U', nssh_tot, 1, M, nssh_tot, ipiv, B, &
         !              &      nssh_tot, work, lwork, info )
         !call sgesv(nssh_tot,1,M,nssh_tot,ipiv,B,nssh_tot,info )
         !print*,'=========== cargamos valores de A y c ==============='
         !do i = 1, nssh_tot
         !  do j = 1, nssh_tot
         !    A(i,j)=M(i,j)
         !  end do
         !  c(i)=B(1,i)
         !end do
   
         !call dgesv(nssh_tot+1,1,M,nssh_tot+1,ipiv,B,nssh_tot+1,info)
         !call sgetrs(nssh_tot,1,M,nssh_tot,ipiv,B,nssh_tot,info )
         
!*
!*     Check for the exact singularity.
!*
               !  write(*,*) 'B output ', B
!      IF( info.GT.0 ) THEN
!         WRITE(*,*)'The element of the diagonal factor '
!         WRITE(*,*)'D(',info,',',info,') is zero, so that'
!         WRITE(*,*)'D is singular; the solution could not be computed.'
!         STOP
!      END IF

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
         alpha = 0
         Ntot = 0.0
         do iatom = 1, natoms
             in1 = imass(iatom)
             do issh = 1, nssh(in1)
                 alpha = alpha + 1
!                 Qout(issh,iatom) = B(1,alpha)
!                 carga(alpha)=B(1,alpha)
                 SQ(alpha) = Qin(issh,iatom)
                 LB(alpha) = 0.00
                 if ( lssh(issh,in1) .eq. 0 ) then 
                   UB(alpha) = 2.00
                   nalpha(alpha) = 1.00 
                 end if
                 if ( lssh(issh,in1) .eq. 1 ) then
                    UB(alpha) = 6.00
                   nalpha(alpha) = 3.00 
                 end if
                 if ( lssh(issh,in1) .eq. 2 ) then
                    UB(alpha) = 10.00
                   nalpha(alpha) = 5.00 
                 end if 
                 !UB(alpha) = 100.00
                 !nalpha(alpha) = 1.00
                 Ntot = Ntot + Qin(issh,iatom)                 
             end do ! end do issh
         end do ! end do iatom
!         print*,''
!         print*,"=============================  solución dgesv ============================="
!         print*,'Qtot =',Ntot 
!         write(*,'(A6,<nssh_tot>F7.3,A3)')'Q = (/ ',(carga(alpha),alpha = 1, nssh_tot),' /)'
!         diff_err = Ep2(carga,A,c,nssh_tot)
!         print*,'err dgesv   = ',diff_err
!
!         alpha = 0
!         do iatom = 1, natoms
!           in1 = imass(iatom)
!           do issh = 1,nssh(in1)
!            alpha=alpha+1
!           ! carga(alpha)=Qneutral(issh,in1)
!           end do
!          end do

         !diff_err=Ep2(carga,A,c,nssh_tot)
         !print*,'err cargas neutras =',diff_err

 
!         Ntot=0
!         do i = 1, nssh_tot
!           Ntot = Ntot + carga(i)
!         end do
         !  SOLVE SYSTEM Mx = B.  x are the charges
         print*,''
         print*,"============  step_size  ============================="
         print*,'Qtot =',Ntot
         write(*,'(A6,<nssh_tot>F7.3,A3)')'Qin = (/ ',(SQ(alpha),alpha = 1, nssh_tot),' /)'
         call step_size(nssh_tot,A,c,SQ, LB, UB, nalpha) !,,LB,UB) B(1,alpha
         Ntot=0
         do alpha = 1, nssh_tot
           Ntot = Ntot + SQ(alpha)
         end do
         !write(*,'(A6,<nssh_tot>F7.3,A3)')'LB =(/ ',(LB(alpha),alpha = 1, nssh_tot),' /)'
         !write(*,'(A6,<nssh_tot>F10.3,A3)')'UB =(/ ',(UB(alpha),alpha = 1, nssh_tot),' /)'
         write(*,'(A6,<nssh_tot>F7.3,A3)')'Qout =(/ ',(SQ(alpha),alpha = 1, nssh_tot),' /)'
         print*,'Qtot(out) =',Ntot
         diff_err=Ep2(SQ,A,c,nssh_tot,nalpha)
         print*,'err step_size   =',diff_err
         
         alpha = 0 
         do iatom = 1, natoms
           in1 = imass(iatom)
           do issh = 1, nssh(in1)
             alpha = alpha + 1
             Qout(issh,iatom) = SQ(alpha) 
           end do ! end do issh
         end do ! end do iatom
         print*,"===================================="

!         alpha = 0
!         do iatom = 1, natoms
!           in1 = imass(iatom)
!           do issh = 1,nssh(in1)
!            alpha=alpha+1
!            carga(alpha)=Qneutral(issh,in1)
!           end do
!          end do
!
!         !diff_err=Ep2(carga,A,c,nssh_tot)
!         !print*,'err cargas neutras =',diff_err
!
! 
!         Ntot=0
!         do i = 1, nssh_tot
!           Ntot = Ntot + carga(i)
!         end do
!         !  SOLVE SYSTEM Mx = B.  x are the charges
!         print*,''
!         print*,"=============================  step_size  ============================="
!         print*,'Qtot =',Ntot
!         write(*,'(A6,<nssh_tot>F7.3,A3)')'Qin=(/ ',(carga(alpha),alpha = 1, nssh_tot),' /)'
!         call step_size(nssh_tot,A,c,carga) !,,LB,UB) B(1,alpha
!         Ntot=0
!         do i = 1, nssh_tot
!           Ntot = Ntot + carga(i)
!         end do
!         print*,'Qtot(out) =',Ntot
!         write(*,'(A6,<nssh_tot>F7.3,A3)')'Q=(/ ',(carga(alpha),alpha = 1, nssh_tot),' /)'
!         diff_err=Ep2(carga,A,c,nssh_tot)
!         print*,'err step_size   =',diff_err
!

     
!          write(*,*) 'Now Ntot is: ', Ntot
      ! renorm total charge  (ztot = total charge)
!         do iatom = 1, natoms
!            in1 = imass(iatom)
!             do issh = 1, nssh(in1)
                 !alpha = alpha + 1
                 !write(*,*) 'indices: ', iatom, issh, alpha, B(1,alpha)
                 !Qout(issh,iatom) = Qout(issh,iatom) + (ztot-Ntot)/nssh_tot
!            end do ! end do issh
!         end do ! end do iatom
        !  write(*,*) 'Ntot = ' , Ntot  
!          write(*,*) 'Qoutno corrected = ', Qout 
! ******************************************************************************
! ******************************************************************************


!!Check whether there are negative charges and correct
!!If there's more than one shell whose charge is negative, more work is
!!needed, but that'd be quite pathological a situation...
!           do iatom = 1,natoms
!            in1 = imass(iatom)
!            do issh = 1, nssh(in1)
!
!               if( Qout(issh,iatom) .lt. 0 .and. nssh(in1) .gt. 1 ) then
!           
!                  do jssh = 1, nssh(in1)
!
!                     if ( jssh .ne. issh ) then
!
!                        Qout(jssh,iatom) = Qout(jssh,iatom)+            &
!                &                         Qout(issh,iatom)/(nssh(in1)-1)
!
!                     end if !end if jssh .ne. issh 
!
!                  end do !end if jssh = 1,nssh(in1)
!
!                  Qout(issh,iatom) = 0.0d0               
!
!               end if !end if  Qout(issh,iatom) .lt. 0
!
!            end do !end do issh = 1, nssh(in1)
!          end do ! end do iatom = 1,natoms
!            ! write(*,*) 'Qout corrected = ', Qout 
      end if  !end if (iqout .eq. 6)
!
!                 CORRECT NEGATIVE CHARGES BY FIXING THEM TO ZERO AND
!                 RESOLVING THE SYSTEM
! ****************************************************************************
!
! C O M P U T E    M U L L I K E N    P O P U L A T I O N    F O R   MOs
! 
! ****************************************************************************

        if (iwrtpop .eq. 1) then
! It is probable that you would like to know about charge localization on
! certain sites, how much eigenvector there is in orbital mu, etc.
! The MOs are selected from predefined energy windows
! see J.P. Lewis et al, J. Phys. Chem. B 107, 2581 (2003)
! Note: from principles it should be applied only to G k-point case

         open (unit = 34, file = 'populations.dat', status = 'unknown')
         write (34,*) ' # Charge localizations (per level): '

! Loop over states
         do iband = 1,  norbitals_new

! Reset charge per atom
          pqmu = 0.0d0

! Loop over the special k points.
          do ikpoint = 1, nkpoints

! is MO in the given energy range?
           if ((eigen_k(iband, ikpoint) .gt. Epop_L) .and. (eigen_k(iband, ikpoint) .lt. Epop_U)) then 

! Loop over the atoms and the neighbors
            do iatom = 1, natoms
             in1 = imass(iatom)
             do ineigh = 1, neighn(iatom)
              mbeta = neigh_b(ineigh,iatom)
              jatom = neigh_j(ineigh,iatom)
              in2 = imass(jatom)
              vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)

! calculate phase factor
              dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2)  &
     &                                      + special_k(3,ikpoint)*vec(3)
              phase = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)

! Loop over all bands
              if (icluster .ne. 1) then
               do imu = 1, num_orb(in1)
                mmu = imu + degelec(iatom)
                step1 = phase*(bbnkre(mmu,iband,ikpoint)                       &
     &                      - ai*bbnkim(mmu,iband,ikpoint))
                do inu = 1, num_orb(in2)
                 nnu = inu + degelec(jatom)
                 step2 = step1*(bbnkre(nnu,iband,ikpoint)                      &
     &                       + ai*bbnkim(nnu,iband,ikpoint))
                 gutr = real(step2)
                 pqmu(iatom) = pqmu(iatom) + gutr*s_mat(imu,inu,ineigh,iatom)

! Finished loop over bands
                end do  ! inu
               end do ! imu
              else ! icluster
               do imu = 1, num_orb(in1)
                mmu = imu + degelec(iatom)
                step1 = phase*bbnkre(mmu,iband,ikpoint)
                do inu = 1, num_orb(in2)
                 nnu = inu + degelec(jatom)
                 step2 = step1*bbnkre(nnu,iband,ikpoint)
                 gutr = real(step2)
                 pqmu(iatom) = pqmu(iatom) + gutr*s_mat(imu,inu,ineigh,iatom)

! Finished loop over bands
                end do ! inu
               end do ! imu
              end if ! icluster

! Finished loop over atoms
             end do  ! ineigh           
            end do ! iatom

! Calculate localization factor of MOs W(mu)
!            checksum = 0.0d0
!            do iatom = 1, natoms
!             checksum = checksum + pqmu(iatom) 
!            end do
!            pqmu = pqmu / checksum
!            checksum = 0.0d0
!            do iatom = 1, natoms
!             checksum = checksum + pqmu(iatom) 
!            end do
 
            checksum = 0.0d0
            do iatom = 1, natoms
             checksum = checksum + pqmu(iatom)*log(abs(pqmu(iatom)))
            end do
            Wmu = exp(-1.0d0*checksum)
            write (34,301) iband,eigen_k(iband,ikpoint),Wmu 

           end if ! eigen_k

! Finished loop over special k-points
          end do  ! ikpoint 
! Finished loop over MOs
         end do ! iband

         close (unit = 34)

        end if ! end if of iwrt_pop = 1

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
201     format (' Band n = ', i4, ' foccupy = ', f18.8)
300     format (2x, ' This is band number: ',2x, i6)
301     format (2x, i8, 2f16.8)
800     format (2x,4i3,f12.6)
        return
      end subroutine denmat


!               allocate( M(nssh_tot+1,nssh_tot+1) )
!               allocate( B(1,nssh_tot+1) )
! call step_size(nssh_tot,M,B) !,,LB,UB) B(1,alpha)
!==================================================================
      subroutine step_size(nssh_tot,A,c,Q,LB,UB,nalpha) !,LB,UB,nstep)
        integer, intent(in) :: nssh_tot
        real,intent(in),dimension(nssh_tot,nssh_tot) :: A
        real,intent(in),dimension(nssh_tot) :: c
        real,intent(inout),dimension(nssh_tot) :: Q
        real,intent(in),dimension(nssh_tot) :: LB, UB, nalpha
!        real,dimension(nssh_tot) :: LB, UB
!        integer, intent(inout) :: nstep
        integer :: nstep = 0
        logical, dimension(nssh_tot) :: cero
        integer :: i, j, k,nceros
        real :: err0,err1,dp,g0,dqmax,x,y1,y2,y3,gmod,qmod
        real,dimension(nssh_tot) ::Q0,g,fc
        real,dimension(nssh_tot,nssh_tot) :: F
        real :: dq , diff_err
        logical :: g_err, igualceros
        !========== control ==========
        real :: tol_err=0.001, nstepmax=2000
        real :: tol = 0.001
        !=============================
        Q0=Q
!        write(*,'(A12,<nssh_tot>F7.3,A3)')'Q in =(/ ',(Q(i),i = 1, nssh_tot),' /)'
        F=0
        do i=1,nssh_tot
          do j=1, nssh_tot
            do k=1, nssh_tot
               F(i,j)=F(i,j)+2*A(k,i)*nalpha(k)*A(k,j)
            end do
          end do
        end do

        fc=0
        do i=1,nssh_tot
          do j=1, nssh_tot
            fc(i)=fc(i)-2*c(j)*nalpha(j)*A(j,i)
          end do
        end do

        nstep=0
        diff_err=100.0
        do while (diff_err > tol_err ) ! .and. nstep .lt. nstepmax) 
          nstep=nstep+1
          ! proyectamos en plano Qtot = cte
          g=0.00
          do i=1,nssh_tot
            do j=1, nssh_tot
              g(i)=g(i)-F(i,j)*Q(j)
            end do
            g(i)=g(i)-fc(i)
          end do


          g0=0.00
          do i=1,nssh_tot
              g0=g0+g(i)
          end do

          g0=g0/nssh_tot

          do i=1,nssh_tot
            g(i)=g(i)-g0
          end do

            !proyectamos en las componentes q<0

          do i=1,nssh_tot
            if((Q0(i) < LB(i)+tol .and. g(i) < 0.00 ) .or. (Q0(i) > UB(i)-tol .and. g(i) > 0.00 )) then
              nceros=nceros+1
              g(i)=0.00
              cero(i)=.True.
            else
              cero(i)=.False.
              g0=g0+g(i)
            end if
          end do
          igualceros=.False.
          do while (igualceros .eq. .False.)
            call getceros(nssh_tot,Q0,g,LB,UB,tol,cero,igualceros)
          end do

          nceros=0
          g0=0.00
          do i=1,nssh_tot
            if(cero(i) .eq. .True.) then
              nceros=nceros+1
              g(i)=0.00
            else
              g0=g0+g(i)
            end if
          end do

          g0=g0/(nssh_tot-nceros)

          do i=1,nssh_tot
            if (cero(i) .eq. .False.) then
              g(i)=g(i)-g0
            else
              g(i)=0.00
            end if
          end do

          gmod=0.00
          do i=1,nssh_tot
            gmod=gmod+g(i)**2
          end do
          gmod=gmod**0.5

          qmod=0.00
          do i=1,nssh_tot
            qmod=qmod+Q(i)**2
          end do

          qmod=qmod**0.5
          dq=qmod/gmod
          !obtenemos el min
          x=get_min_parabola(-dq,0,dq,Ep2(Q0-dq*g,A,c,nssh_tot,nalpha),Ep2(Q0,A,c,nssh_tot,nalpha),Ep2(Q0+dq*g,A,c,nssh_tot,nalpha))
          Q=Q0+x*g

          diff_err=Ep2(Q0,A,c,nssh_tot,nalpha)-Ep2(Q,A,c,nssh_tot,nalpha)
         dqmax=1.0E+10
          k=0
          do i=1,nssh_tot
            if(cero(i) .eq. .False.) then  !ojo !!!!
              if(Q(i) < LB(i)) then
                if (dqmax > -Q0(i)/g(i)) then
                  dqmax=-Q0(i)/g(i)
                  k=k+1
                end if
              end if
              if(Q(i) > UB(i)) then
                if (dqmax > (UB(i)-Q0(i))/g(i)) then
                  dqmax=(UB(i)-Q0(i))/g(i)
                  k=k+1
                end if
              end if
            end if
          end do
          if(k .eq.0 ) then
            Q=Q
          else
            Q=Q0+dqmax*g
          end if

      !  aux=0
      !  do i = 1, nssh_tot
      !    aux=aux+g(i)**2
      !  end do
      !  aux=aux**0.5

      !  aux2=0
      !  do i = 1, nssh_tot
      !    aux2=aux2+abs(Q0(i)-Q(i))
      !  end do


      !     diff=diffQ(nssh_tot,Q0,Q)
      !    print*,Ep2(Q,A,c,nssh_tot),aux,aux2
!        write(*,'(A6,<nssh_tot>F7.3,A3)')'Q     = (/ ',(Q(i),i = 1, nssh_tot),' /)'
!        print*,Ep2(Q,A,c,nssh_tot)
          Q0=Q

        !write(*,'(A6,<nssh_tot>F7.3,A3)')'g     = (/ ',(g(i),i = 1, nssh_tot),' /)'
        end do

        aux=0
        do i = 1, nssh_tot
          aux=aux+Q(i)
        end do
      
     print*,'nstep = ',nstep
 
      !  write(*,'(A6,<nssh_tot>L7,A3)')'ceros =(/ ',(cero(i),i = 1, nssh_tot),' /)'
      !  do i = 1, nssh_tot
      !    print*,g(i),cero(i),Q(i),dqmax,k,x
      !  end do
      
      !  print*,'diff err = ',diff_err,';  nstep    =',nstep,'; dqmax =',dqmax
      !  print*,'Qtot     =',aux,' ;  n ceros  =',nceros,'; err   =',Ep2(Q,A,c,nssh_tot)
      !  if ( nstep .eq. nstepmax ) then
      !    print*,'*************** stop nstepmax *************************'
      !    stop
      !  end if
      end
      !=============================================
      real function Ep2(q,A,c,nssh_tot,nalpha)
        integer, intent(in) :: nssh_tot
        real,intent(in),dimension(nssh_tot,nssh_tot) :: A
        real,intent(in),dimension(nssh_tot) :: q,c,nalpha
        integer i,j,k
        real aux
        !Ep=(Aq-c)**2
        Ep2=0.00
        do i = 1, nssh_tot
          aux=0.00
          do j = 1, nssh_tot
            aux=aux+A(i,j)*q(j)
          end do
          aux=aux-c(i)
          aux=aux**2*nalpha(i)
          Ep2=Ep2+aux
        end do
      end function Ep2
 
      subroutine getceros(nssh_tot,Q0,g,LB,UB,tol,cero,igualceros)
      integer, intent(in) :: nssh_tot
      real, intent(in) :: tol
      real,intent(in),dimension(nssh_tot) :: Q0
      real,intent(in),dimension(nssh_tot) :: LB, UB
      real,intent(in),dimension(nssh_tot) :: g
      
      logical,intent(inout) ,dimension(nssh_tot) :: cero
      logical, intent(inout) :: igualceros
      
      real,dimension(nssh_tot) :: gaux
      real :: g0
      integer :: i, j,nceros
      real :: aux
          nceros=0
          g0=0.00
          gaux=g
          !write(*,'(A6,<nssh_tot>F12.3,A3)')'g0=(/ ',(gaux(i),i = 1, nssh_tot),' /)'
          !write(*,'(A6,<nssh_tot>F12.3,A3)')'Q0=(/ ',(Q0(i),i = 1, nssh_tot),' /)'
          do i=1,nssh_tot
            if(cero(i) .eq. .True.) then
              nceros=nceros+1
              gaux(i)=0.00
            else
              g0=g0+g(i)
            end if
          end do
      
          g0=g0/(nssh_tot-nceros)
      
          do i=1,nssh_tot
            if (cero(i) .eq. .False.) then
              gaux(i)=gaux(i)-g0
            else
              gaux(i)=0.00
            end if
          end do

        igualceros=.True.
         do i=1,nssh_tot
           if(cero(i)  .eq. .False.)then
              if((Q0(i) < LB(i)+tol .and. gaux(i) < 0.00 ) .or. (Q0(i) > UB(i)-tol .and. gaux(i) > 0.00 )) then
                nceros=nceros+1
                cero(i)=.True.
                igualceros=.False.
              end if
           end if
         end do


      ! write(*,'(A6,<nssh_tot>L)')'cein =(/ ',(cero(i),i = 1, nssh_tot),' /)'
      ! print*,'igualceros=',igualceros
      end subroutine getceros
 
      real function get_min_parabola(x1,x2,x3,y1,y2,y3)
        real,intent(in) :: x1,x2,x3,y1,y2,y3
        real a,b,c
        b=(y2-y3)-(y1-y2)*(x2**2-x3**2)/(x1**2-x2**2)
        b=b/( (x2-x3)-(x1-x2)*(x2**2-x3**2)/(x1**2-x2**2) )
        a=(y1-y2)/(x1**2-x2**2)-b*(x1-x2)/(x1**2-x2**2)
        c=y1-a*x1**2-b*x1
        get_min_parabola=-b/(2*a)
        ! print*,a,'*x**2+',b,'*x+',c
      end function get_min_parabola

!=================================================
 
      
