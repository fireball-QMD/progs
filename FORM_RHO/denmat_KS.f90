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
        integer issh
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
        real, dimension (norbitals, nkpoints) :: foccupy
        real, dimension (numorb_max, natoms) :: QMulliken
        real, dimension (3) :: vec
 
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
        write (*,*) ' Fermi Level = ', efermi
 
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
        if (iqout .eq. 1) then
         Qout = 0.0d0
         QLowdin_TOT = 0.0d0 

         if (ifixcharge .eq. 1) then

          do iatom = 1, natoms
           in1 = imass(iatom)
           do imu = 1, num_orb(in1)
            Qout(imu,iatom) = Qin(imu,iatom)
            QLowdin_TOT(iatom) = QLowdin_TOT(iatom) + Qin(imu,iatom) 
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
          end do ! do iatom

         end if      ! endif of ifixcharges
        end if       ! endif of iqout .eq. 1

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
           do imu = 1, num_orb(in1)
            Qout(imu,iatom) = Qin(imu,iatom) 
            QMulliken_TOT(iatom) = QMulliken_TOT(iatom) + Qin(imu,iatom) 
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
           do imu = 1, num_orb(in1)
             Qout(imu,iatom) = Qout(imu,iatom) + QMulliken(imu,iatom)
             QMulliken_TOT(iatom) = QMulliken_TOT(iatom) + Qout(imu,iatom) 
           end do
 
! End loop over atoms
          end do
         end if     ! endif of ifixcharges
        end if      ! endif of iqout .eq. 2
 

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
!120     format (2x, <norbitals>f10.6) 
400     format (2x, 'Qmull =',10f10.6)

        return
      end subroutine denmat_KS
        
