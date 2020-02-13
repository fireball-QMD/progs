

! copyright info:
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
 
        integer, dimension (norbitals) :: ioccupy
        integer, dimension (norbitals, nkpoints) :: ioccupy_k
 
        real aux1, aux2, aux3
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
 
        complex ai
        complex phase, phasex
        complex step1, step2
 
        logical read_occupy

! A bunch of memory to be used in many ways
        integer jnu,jmu
        complex*16, dimension (:, :), allocatable :: xxxx
        complex*16, dimension (:, :), allocatable :: yyyy
        complex*16 a0
        complex*16 a1
        real, dimension (3) :: k_temp

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
 
! ****************************************************************************
!
!  C O M P U T E    L O W D I N    C H A R G E S
! ****************************************************************************
! Initialize
!        if (iqout .eq. 1) then
         Qout = 0.0d0
         QLowdin_TOT = 0.0d0 
         if(.not. allocated(blowre)) allocate (blowre (norbitals, norbitals, nkpoints))

         if (ifixcharge .eq. 1) then

          do iatom = 1, natoms
           in1 = imass(iatom)
           do imu = 1, num_orb(in1)
            Qout(imu,iatom) = Qin(imu,iatom)
            QLowdin_TOT(iatom) = QLowdin_TOT(iatom) + Qin(imu,iatom) 
           end do
          end do

         else 
          open(unit = 101, file = 'CHARGES_LOWD_KS', status = 'unknown')
          do iatom = 1, natoms
           in1 = imass(iatom)
! Loop over the special k points.
           do ikpoint = 1, nkpoints
           k_temp(:) = special_k(:,ikpoint)
           !-------------------iqout
           call calc_blowreim (1, icluster, ikpoint, nkpoints, k_temp)
           !deallocate(blowre)
            aux1 = weight_k(ikpoint)*spin
            do iorbital = 1, norbitals
             if (ioccupy_k(iorbital,ikpoint) .eq. 1) then
              aux2 = aux1*foccupy(iorbital,ikpoint)
! Finally the imu loop.
              do imu = 1, num_orb(in1)
                mmu = imu + degelec(iatom)
                if (icluster .ne. 1) then
                 aux3 = aux2*(blowre(mmu,iorbital,ikpoint)**2                &
     &                        + blowim(mmu,iorbital,ikpoint)**2)
                else
                 aux3 = aux2*blowre(mmu,iorbital,ikpoint)**2
                end if 
                Qout(imu,iatom) = Qout(imu,iatom) + aux3
                QLowdin_TOT(iatom) = QLowdin_TOT(iatom) + aux3 
              end do ! do imu
             end if
 
! End loop over orbitals and kpoints
            end do ! do iorbital
           end do ! do ikpoints
! End loop over atoms
          print*,'Write CHARGES_LOWD_KS file'
          !write(101,120) iatom, (Qout(imu,iatom), imu = 1, num_orb(in1))
          write(101,121) iatom,QLowdin_TOT(iatom) , (Qout(imu,iatom), imu = 1, num_orb(in1))
          end do ! do iatom
          close(101)

         Qout = 0.0d0
         QLowdin_TOT = 0.0d0 
          open(unit = 103, file = 'CHARGES_NAP_KS', status = 'unknown')
          do iatom = 1, natoms
           in1 = imass(iatom)
! Loop over the special k points.
           do ikpoint = 1, nkpoints
           k_temp(:) = special_k(:,ikpoint)
           !-------------------iqout
           call calc_blowreim (3, icluster, ikpoint, nkpoints, k_temp)
           !deallocate(blowre)
            aux1 = weight_k(ikpoint)*spin
            do iorbital = 1, norbitals
             if (ioccupy_k(iorbital,ikpoint) .eq. 1) then
              aux2 = aux1*foccupy(iorbital,ikpoint)
! Finally the imu loop.
              do imu = 1, num_orb(in1)
                mmu = imu + degelec(iatom)
                if (icluster .ne. 1) then
                 aux3 = aux2*(blowre(mmu,iorbital,ikpoint)**2                &
     &                        + blowim(mmu,iorbital,ikpoint)**2)
                else
                 aux3 = aux2*blowre(mmu,iorbital,ikpoint)**2
                end if 
                Qout(imu,iatom) = Qout(imu,iatom) + aux3
                QLowdin_TOT(iatom) = QLowdin_TOT(iatom) + aux3 
              end do ! do imu
             end if
 
! End loop over orbitals and kpoints
            end do ! do iorbital
           end do ! do ikpoints
! End loop over atoms
          print*,'Write CHARGES_NAP_KS file'
          !write(103,120) iatom, (Qout(imu,iatom), imu = 1, num_orb(in1))
          write(103,121) iatom,QLowdin_TOT(iatom) , (Qout(imu,iatom), imu = 1, num_orb(in1))
          end do ! do iatom
          close(103)


         end if      ! endif of ifixcharges
!        end if       ! endif of iqout .eq. 1

! ****************************************************************************
!
!  C O M P U T E    M U L L I K E N    C H A R G E S
! ****************************************************************************
! Compute Mulliken charges.
!        if (iqout .eq. 2) then

         Qout = 0.0d0
         QMulliken = 0.0d0
         QMulliken_TOT = 0.0d0

         if (ifixcharge .eq. 1) then

          do iatom = 1, natoms
           in1 = imass(iatom)
           QMulliken_TOT(iatom) = 0.0d0
           do imu = 1, num_orb(in1)
            Qout(imu,iatom) = Qin(imu,iatom) 
            QMulliken_TOT(iatom) = QMulliken_TOT(iatom) + Qin(imu,iatom) 
           end do
          end do

         else
          print*,'Write CHARGES_MULL_KS file'
          open(unit = 102, file = 'CHARGES_MULL_KS', status = 'unknown')
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
           do imu = 1, num_orb(in1)
             Qout(imu,iatom) = Qout(imu,iatom) + QMulliken(imu,iatom)
             QMulliken_TOT(iatom) = QMulliken_TOT(iatom) + Qout(imu,iatom) 
           end do
 
! End loop over atoms
         ! write(102,120) iatom, (Qout(imu,iatom), imu = 1, num_orb(in1))
          write(102,121) iatom, QMulliken_TOT(iatom) , (Qout(imu,iatom), imu = 1, num_orb(in1))
          end do
          close(102)
         end if     ! endif of Gfixcharges
!        end if      ! endif of iqout .eq. 2

! ****************************************************************************
!
!  C O M P U T E    M U L L I K E N - D I P O L E    C H A R G E S
! ****************************************************************************
! Compute Mulliken-dipole charges.

         Qout = 0.0d0
         QMulliken = 0.0d0
         QMulliken_TOT = 0.0d0

         if (ifixcharge .eq. 1) then
          do iatom = 1, natoms
           in1 = imass(iatom)
           QMulliken_TOT(iatom) = 0.0d0
           do imu = 1, num_orb(in1)
            Qout(imu,iatom) = Qin(imu,iatom)
            QMulliken_TOT(iatom) = QMulliken_TOT(iatom) + Qin(imu,iatom)
           end do
          end do

         else
          print*,'Write CHARGES_MUDI_KS file'
          open(unit = 104, file = 'CHARGES_MUDI_KS', status = 'unknown')

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
           do imu = 1, num_orb(in1)
             Qout(imu,iatom) = QMulliken(imu,iatom)
             QMulliken_TOT(iatom) = QMulliken_TOT(iatom) + Qout(imu,iatom)
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
          !write(104,120) iatom, (Qout(imu,iatom), imu = 1, num_orb(in1))
          write(104,121) iatom,QMulliken_TOT(iatom) , (Qout(imu,iatom), imu = 1, num_orb(in1))
          end do
          close(104)
         end if     ! endif of ifixcharges

 

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
120     format (2x, i4, <norbitals>f10.6) 
121     format (2x, i4,f10.6,<norbitals>f10.6) 
400     format (2x, 'Qmull =',10f10.6)

        return
      end subroutine denmat_KS
       




! ===========================================================================
! ===========================================================================
! ===========================================================================
! ===========================================================================
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
print*,ikpoint,nkpoints
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

print*,'CHECK THE LINEAR DEPENDENCE'
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

print*,'CHECK THE LINEAR DEPENDENCE',slam
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
print*,'CALCULATE (S^-1/2) --> sm1'
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

print*,'slam',slam

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




 
