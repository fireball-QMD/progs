! copyright info:
!
!                             @Copyright 2001
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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! mdetdenmat.f90
! Program Description
! ===========================================================================
!       This routine calculates the density matrices and the band-structure
! energy, as well as similar density matrices.
!
! JOM-info : changed for MDET (from denmat)
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
        subroutine mdetdenmat (ifixcharge, iqout, icluster, iwrtefermi, &
     &                  tempfe, ebs, iwrtpop)
        use charges
        use configuration
        use constants_fireball
        use density
        use dimensions
        use interactions
        use neighbor_map
        use kpoints
        use nonadiabatic ! vlada

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icluster
        integer, intent (in) :: ifixcharge
        integer, intent (in) :: iqout
        integer, intent (in) :: iwrtefermi
        integer, intent (in) :: iwrtpop

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
	

        real, dimension (norbitals, nkpoints) :: foccupy
	
	
        real, dimension (numorb_max, natoms) :: QMulliken
        real, dimension (3) :: vec

        complex ai
        complex phase, phasex
        complex step1, step2

        logical read_occupy

! Procedure
! ===========================================================================
! Initialize some things
        ai = cmplx(0.0d0,1.0d0)
! jel: where you set rho to zero ??? JOM: in build_rho
        rhoPP = 0.0d0

        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' ************************************************ '
        write (*,*) '  '
        write (*,*) '                Welcome to denmat --              '
        write (*,*) '  '
        write (*,*) ' ************************************************ '

! ****************************************************************************
!
!                     C H A R G E    O C C U P A T I O N S
! ****************************************************************************
! If there exists a file 'OCCUPATION', then make a list of energy eigenvalues
! you want to 'unoccupy', and give the energy shift to put the eigenvalue
! above E-fermi.
! JOM-info : we will not use this in MDET (I think)
!       inquire (file = 'OCCUPATION', exist = read_occupy)
!       if (read_occupy) then
!
! Open the file and read information.
!        open (unit = 22, file = 'OCCUPATION', status = 'old')
!        write (*,*) '  '
!        write (*,*) ' Reading from the OCCUPATION file! '
!        read (22,*) noccupy
!        if (noccupy .gt. norbitals) then
!         write (*,*) ' noccupy > norbitals: from OCCUPATION file. '
!         stop
!        end if
!        do imu = 1, noccupy
!         read (22,*) iband, deltae
!         eigen_k(iband,1:nkpoints) = eigen_k(iband,1:nkpoints) + deltae
!         ioccupy(imu) = iband
!        end do
!        close (unit = 22)

!       end if

! ****************************************************************************
! JOM-info : we will not use this in MDET (I think)


! Get the Fermi energy.
!       call fermie (norbitals, ztot, eigen_k, efermi, ioccupy_k, foccupy)
!       write (*,*) ' Fermi Level = ', efermi
! JOM
!if (trans .eq. 1) then
!        foccupy = foccupy_na_TS
!else
        foccupy = foccupy_na
!end if
! JOM-info : ioccupy_na can be 0, 1 or 2, but ioccupy_k is only 0 or 1

     ! write (91244,*)  itime_step
 !     write (91244,'(<norbitals>i6.1)') (ioccupy_na(inu,1), inu = 1, norbitals)

     ! write (912444,*)  
   !   write (912444,'(<norbitals>i6.1)') (ioccupy_na_TS(inu,1), inu = 1, norbitals)


        ioccupy_k = 1
        do ikpoint = 1, nkpoints
         do iband = 1, norbitals_new
          if (ioccupy_na(iband,ikpoint) .eq. 0 ) then
     !     if (ioccupy_na_TS(iband,ikpoint) .eq. 0 ) then
                 ioccupy_k(iband,ikpoint) = 0
      !     end if 
          end if
         end do
        end do

        if (iwrtefermi .eq. 1) then
         write (*,*) '  '
         write (*,*) ' We write out the occupancies of the levels  '
         write (*,*) ' ----- ioccupy_k, from ioccupy_na '
         do ikpoint = 1, nkpoints
          write (*,*) '  '
          write (*,*) ' ------ fermi ioccupy_k for k-point = ', ikpoint
          do iband = 1, norbitals_new
           write (*,200) iband, ioccupy_k(iband,ikpoint)
          end do
         end do

         write (*,*) '  '
         write (*,*) 'We write out the electron fraction in each level '
         write (*,*) ' ------ from the foccupy_na. '
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
        if (iqout .ne. 2) then
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


! ****************************************************************************
!
! C O M P U T E    M U L L I K E N    P O P U L A T I O N    F O R    H O M O
! JPL - Generalize this later!
! ****************************************************************************

        if (iwrtpop .eq. 1) then
! It is probable that you would like to know about charge localization on
! certain sites, how much eigenvector there is in orbital mu, etc.
! FOR NOW - ONLY DO HOMO and ONLY FOR DNA!!!
        open (unit = 34, file = 'populations.dat', status = 'unknown')
        write (34,*) '  '
        write (34,*) ' Charge localizations (per level): '
        write (34,*) ' **************************************************** '
        do iband = norbitals_new, 1, -1
         write (34,300) iband

! Loop over the atoms and the neighbors
         do iatom = 1, natoms
          in1 = imass(iatom)
          pcharge = 0.0d0
          do ineigh = 1, neighn(iatom)
           mbeta = neigh_b(ineigh,iatom)
           jatom = neigh_j(ineigh,iatom)
           in2 = imass(jatom)
           vec = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
!
! Loop over the special k points.
           do ikpoint = 1, nkpoints
            dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2)  &
     &                                        + special_k(3,ikpoint)*vec(3)
            phase = cmplx(cos(dot),sin(dot))*weight_k(ikpoint)

! Loop over all bands
            if (icluster .ne. 1) then
             do imu = 1, num_orb(in1)
              mmu = imu + degelec(iatom)
              step1 = phase*(bbnkre(mmu,iband,ikpoint)                       &
     &                       - ai*bbnkim(mmu,iband,ikpoint))
              do inu = 1, num_orb(in2)
               nnu = inu + degelec(jatom)
               step2 = step1*(bbnkre(nnu,iband,ikpoint)                      &
     &                        + ai*bbnkim(nnu,iband,ikpoint))
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
!
! Finished loop over kpoints
           end do

! Finish remaining loops
          end do
          write (34,301) iatom, pcharge
         end do
         write (34,*) ' **************************************************** '
        end do
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
201     format (' Band n = ', i4, ' foccupy = ', f12.8)
300     format (2x, ' This is band number: ',2x, i6)
301     format (2x, i4, f10.6)
800     format (2x,4i3,f12.6)
        return
      end subroutine mdetdenmat

