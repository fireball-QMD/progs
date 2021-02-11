! copyright info:
!
!                             @Copyright 2002
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
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

 
! assemble_ca_3c.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the true three-center interactions.
! These are true three-center in that iatom .ne. jatom .ne. katom.
! The matrix elements look like <Psi(1)|V(3)|Psi(2)>.
!
! The bulk of the work is done in trescentros.f. This program assembles the
! results.
!
! A third party term is when we consider the NA (or etc.) to be at the origin
! and take the matrix element between a pair of neighbors, neither of which is
! the NA (or etc.), and the pair is truly a pair, and not an atom.
!
! The generic form of the matrix element is
! v = <mu(r-(r1+ratm)) ! v(r-ratm) ! nu(r-(r2+ratm))>.
!
! See notes
! "general theory of third party terms" for the definition on p. 2 of
! these three terms.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658 
! Provo, UT 841184602-4658
! FAX 801-422-2265
! Office telephone 801-422-7444
! ===========================================================================
 
! Program Declaration
! ===========================================================================
        subroutine assemble_ca_3c (nprocs, iordern, igauss)
        use charges
        use configuration
        use constants_fireball
        use dimensions
        use interactions
        use neighbor_map
        use gaussG
        use options, only : iqout
        use scf, only : Kscf

        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: igauss
        integer, intent (in) :: iordern
        integer, intent (in) :: nprocs
 
!$ volatile rcutoff
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ialp
        integer iatom
        integer iatomstart
        integer ibeta
        integer icount
        integer icount_sav
        integer ierror
        integer imu
        integer in1
        integer in2
        integer indna
        integer ineigh
        integer interaction
        integer inu
        integer isorp
        integer issh
        integer jatom
        integer jneigh
        integer jbeta
        integer jcount
        integer jcount_sav
        integer jssh
        integer matom
        integer mbeta
        integer mneigh
        integer my_proc
        integer natomsp
 
        real cost
        real distance_13
        real distance_23
        real dq3
        real dstn_temp
        real dterm
        real dxn
        real rcutoff_ialp
        real rend1
        real rend2
        real sterm
        real stn_temp1
        real stn_temp2
        real x
        real y
        real rcutoff_i
        real rcutoff_j
 
        real, dimension (numorb_max, numorb_max) :: bcca
        real, dimension (numorb_max, numorb_max) :: bccax
        real, dimension (numorb_max, numorb_max) :: emnpl
        real, dimension (numorb_max, numorb_max) :: emnpl_noq
        real, dimension (3, 3, 3) :: deps
        real, dimension (3, 3) :: eps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: rhat
        real, dimension (3) :: rna
        real, dimension (3) :: rnabc
        real, dimension (3) :: sighat
        real, dimension (numorb_max, numorb_max) :: stn1
        real, dimension (numorb_max, numorb_max) :: stn2

! jel: dynamical use of arrays listed bellow fixes stack size problem occured on some machines
!        real, dimension (numorb_max, numorb_max) :: smG
!        real, dimension (3, numorb_max, numorb_max) :: spmG
!        real, dimension (numorb_max, numorb_max, neigh_max, natoms) :: smatG
!        real, dimension (numorb_max, numorb_max, neigh_max, natoms) :: spmatG
        real, dimension (:,:), allocatable :: smG
        real, dimension (:,:,:), allocatable :: spmG
        real, dimension (:, :, :, :), allocatable :: smatG
        real, dimension (:, :, :, :), allocatable :: spmatG
      
! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ========================================================================
! The arrays ewaldsr, vna, and vca will not be initialized to zero here.
! Presumably, the two-center interactions have already been calculated.
! As such, at this point and time these arrays should not be zero.
!
! Determine which atoms are assigned to this processor.
         allocate (smG (numorb_max, numorb_max))
         allocate (spmG (3, numorb_max, numorb_max))
         allocate (smatG (numorb_max, numorb_max, neigh_max, natoms))
         allocate (spmatG (numorb_max, numorb_max, neigh_max, natoms))

         smG = 0.0d0
         smatG = 0.0d0
         spmG = 0.0d0
         spmatG = 0.0d0

! Determine which atoms are assigned to this processor.
        if (iordern .eq. 1) then
         call MPI_COMM_RANK (MPI_BTN_WORLD, my_proc, ierror)
         natomsp = natoms/nprocs
         if (my_proc .lt. mod(natoms,nprocs)) then
          natomsp = natomsp + 1
          iatomstart = natomsp*my_proc + 1
         else
          iatomstart = (natomsp + 1)*mod(natoms,nprocs)                 &
                      + natomsp*(my_proc - mod(natoms,nprocs)) + 1
         end if
        else
         iatomstart = 1
         natomsp = natoms
        end if

        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neigh_self(iatom)
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)

! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)          ! <==== loop over i's neighbors
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)

! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
          r21(:) = r2(:) - r1(:)
          y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))

! Find the unit vector in sigma direction.
          if (y .lt. 1.0d-05) then
           sighat(1) = 0.0d0
           sighat(2) = 0.0d0
           sighat(3) = 1.0d0
          else
           sighat(:) = r21(:)/y
          end if

          call epsilon (r2, sighat, eps)
          call deps2cent (r1, r2, eps, deps)

           if (igauss .eq. 1) then
              call doscentrosG_overlap (in1, in2, y, eps, deps, smG,    &
     &                                  spmG, rcutoff)

             do inu = 1, num_orb (in2)
              do imu = 1, num_orb (in1)
               smatG (imu,inu,ineigh,iatom) = smG (imu,inu)
              end do
             end do
            end if

! End loop over iatom and its neighbors - jatom.
         end do
        end do

! Choose atom ialp in the central cell. This is the atom whose position
! we take the derivative, and is the atom who has the the neutral atom
! potential.
! Loop over the atoms in the central cell.
!$omp parallel do private (iatom, ibeta, icount, icount_sav, in1, in2, indna)&
!$omp&            private (interaction, jatom, jbeta, jcount, jcount_sav)    &
!$omp&            private (mneigh, cost, distance_13, distance_23, dq3)      &
!$omp&            private (dstn_temp, dterm, dxn, rcutoff_ialp, rend1, rend2)&
!$omp&            private (sterm, stn_temp1, stn_temp2, x, y, bcca, bccax)   &
!$omp&            private (emnpl, eps, r1, r2, r21, rhat, rna, rnabc, sighat)&
!$omp&            private (stn1, stn2)
        do ialp = iatomstart, iatomstart - 1 + natomsp
         rna(:) = ratom(:,ialp)
         indna = imass(ialp)

! Find the charge of the third-center
         dq3 = 0.0d0
         do issh = 1, nssh(indna)
          dq3 = dq3 + (Qin(issh,ialp) - Qneutral(issh,indna))
         end do
 
         rcutoff_ialp = 0.0d0
         do imu = 1, nssh(indna)
          if (rcutoff(indna,imu) .gt. rcutoff_ialp)                     &
           rcutoff_ialp = rcutoff(indna,imu)
         end do
 
! Loop over the neighbors of each ialp.
! Now look at all common neighbor pairs which were figured out in main.
! The common neighbors were figured out in common_neighbors.f90
         do ineigh = 1, neigh_comn(ialp)
          mneigh = neigh_comm(ineigh,ialp)
 
! The second atom (jatom) is the mneigh'th neighbor of iatom.
          if (mneigh .ne. 0) then
           iatom = neigh_comj(1,ineigh,ialp)
           ibeta = neigh_comb(1,ineigh,ialp)
           r1(:) = ratom(:,iatom) + xl(:,ibeta)
           in1 = imass(iatom)
! ENRIQUE-JOM: calculate rcutoff_i for smoother
          rcutoff_i = 0
          do imu = 1, nssh(in1)
           if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
          end do
! End ENRIQUE-JOM 
 
           jatom = neigh_comj(2,ineigh,ialp)
           jbeta = neigh_comb(2,ineigh,ialp)
           r2(:) = ratom(:,jatom) + xl(:,jbeta)
           in2 = imass(jatom)
           jneigh = neigh_back(iatom,mneigh)
! ENRIQUE-JOM: calculate rcutoff_j for smoother 
          rcutoff_j = 0
          do imu = 1, nssh(in2)
           if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
          end do
! End ENRIQUE-JOM

 
 
! ****************************************************************************
!
! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
           r21(:) = r2(:) - r1(:)
           y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
 
! Find the unit vector in sigma direction.
           if (y .lt. 1.0d-05) then
            sighat(1) = 0.0d0
            sighat(2) = 0.0d0
            sighat(3) = 1.0d0
            write (*,*) ' There is an error here in assemble_3c.f '
            write (*,*) ' r1 = r2!!!! BAD!!!! '
           else
            sighat(:) = r21(:)/y
           end if
 
! Find rnabc = vector pointing from center of bondcharge to the neutral atom.
! The center of the bondcharge is at r1 + r21/2.  This gives us the distance
! dnabc (x value in 2D grid).
           rnabc(:) = rna(:) - (r1(:) + r21(:)/2.0d0)
           x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)
 
! Find the unit vector in rnabc direction.
           if (x .lt. 1.0d-05) then
            rhat(1) = 0.0d0
            rhat(2) = 0.0d0
            rhat(3) = 0.0d0
           else
            rhat(:) = rnabc(:)/x
           end if
           cost = sighat(1)*rhat(1) + sighat(2)*rhat(2) + sighat(3)*rhat(3)
           call epsilon (rhat, sighat, eps)
 
! For now we just do the neutral atom interactions.
! Charged atom interactions are assembled in assemble_ca_3c.f
! So set isorp = 0 within this subroutine.
!
!              interaction    subtypes     index
!
!      bcna         1           0..9(max)   1..10
 
! ****************************************************************************
!
! PERFORM ACTUAL CALCULATIONS
! SET-UP AND ASSEMBLE EWALDSR AND DEMNPL
! ****************************************************************************
! Need direction cosines from atom 1 to ratm (rhatA1),
! and atom 2 to ratm (rhatA2).
           distance_13 = sqrt((rna(1) - r1(1))**2 + (rna(2) - r1(2))**2 &
     &                                           + (rna(3) - r1(3))**2)
           distance_23 = sqrt((rna(1) - r2(1))**2 + (rna(2) - r2(2))**2 &
     &                                           + (rna(3) - r2(3))**2)
 
! Now let's calculate the asymptotic form so that we can match these
! better with the sticky smooters.
! Reminder: <B|A|C>= stn1*exact+(1-stn1)*asympt
! So d/dr = stn1 * dexact + dstn1*exact + (1-stn1)*dasympt + d(1-stn1)*asympt.
!               A              B               C                 D.
! Calculate d(asymptote).
! Note that we are using the effective dipole theory here.
 
! First, calculate the interaction ewaldsr - This is the correction due to what
! is included (more accurately) in the short range terms.
! Second, calculate the long-range effective monopole.  This term is included
! so that we obtain no discontinuities when atoms leave or enter the 2*rc
! range criteria.  Therefore, "close" three-center interactions are exact,
! while more distant three-center integrals go to effective monopoles.  The
! monopoles are effective in the sense that the two atoms in the matrix element,
! each has a different charge.  Since they are separated, this gives a
! monopole + dipole contribution at long range.
! The smoothing function is smoother(r,rbegin,rend).  We define our answer as
! smoother(r)*exact + (1 - smoother(r))*longrange.  The distance r is the
! distance of the third center from the "effective" center of the bondcharge.
! The effective center of the bondcharge is (d + rc1 - rc2)/2 from r1 except in
! wierd cases (see below).
! The distance rbegin is the distance at which we include only exact answers
! and do not smooth. The distance rend is the distance past which smooth(r) is
! zero, so that the result is long-range only.
           icount_sav = 0
           do issh = 1, nssh(in1)
            jcount_sav = 0
            do jssh = 1, nssh(in2)
! ENRIQUE-JOM: use rcutoff_i,rcutoff_j + rcutoff_ialp as the smoother
!             rend1 = rcutoff(in1,issh) + rcutoff_ialp
             rend1 = rcutoff_i + rcutoff_ialp
!             rend2 = rcutoff(in2,jssh) + rcutoff_ialp
             rend2 = rcutoff_j + rcutoff_ialp
             call smoother (distance_13, rend1, smt_elect, stn_temp1, dstn_temp)
             call smoother (distance_23, rend2, smt_elect, stn_temp2, dstn_temp)
! Different centers, multiply smoothers (this way it goes to zero for the 
! the smaller smoother).  Can't just pick smaller rc based one (like 
! two-center case), since distances also vary.
             stn_temp1 = stn_temp1*stn_temp2
             stn_temp2 = 1.0d0 - stn_temp1
             do inu = 1, lssh(issh,in1)*2 + 1
              icount = icount_sav + inu
              do imu = 1, lssh(jssh,in2)*2 + 1
               jcount = jcount_sav + imu
               stn1(icount,jcount) = stn_temp1
               stn2(icount,jcount) = stn_temp2
              end do
             end do
             jcount_sav = jcount
            end do
            icount_sav = icount
           end do

           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             sterm = dq3*s_mat(imu,inu,mneigh,iatom)/2.0d0
             dterm = dq3*dip(imu,inu,mneigh,iatom)/y
             emnpl(imu,inu) = (sterm - dterm)/distance_13               &
     &                       + (sterm + dterm)/distance_23

             emnpl_noq(imu,inu) = ((s_mat(imu,inu,mneigh,iatom)/2.0d0)-&
             & (dip(imu,inu,mneigh,iatom)/y))/distance_13+ &
             &  ((s_mat(imu,inu,mneigh,iatom)/2.0d0)+&
             & (dip(imu,inu,mneigh,iatom)/y))/distance_23

             ewaldsr(imu,inu,mneigh,iatom) =                            &
     &        ewaldsr(imu,inu,mneigh,iatom) + emnpl(imu,inu)*eq2
            !SFIRE
            ewaldsr(inu,imu,jneigh,jatom)=ewaldsr(imu,inu,mneigh,iatom)
            !SFIRE
            
             if (Kscf .eq. 1 .and. iqout .eq. 6) then
            do issh = 1, nssh(indna)
             gvhxc(imu,inu,issh,ialp,mneigh,iatom) = &
           &   gvhxc(imu,inu,issh,ialp,mneigh,iatom) &
              &- emnpl_noq(imu,inu)*eq2
             ! symmetrize
             gvhxc(inu,imu,issh,ialp,jneigh,jatom) = &
           &   gvhxc(imu,inu,issh,ialp,mneigh,iatom)
            end do ! end do issh
           end if ! end if Kscf .eq. 1 .and. iqout .eq. 6
            end do
           end do
 
! ****************************************************************************
!
! PERFORM ACTUAL CALCULATIONS
! CALL TRESCENTROS FOR NEUTRAL ATOM PIECE
! ****************************************************************************
! Initialize bcca
           bcca = 0.0d0
           do isorp = 1, nssh(indna)

           if (igauss .eq. 0) then
           interaction = 1
            call trescentros (interaction, isorp, isorpmax, in1, in2,   &
     &                        indna, x, y, cost, eps, bccax, nspecies)
           end if

          if (igauss .eq. 1) then
             call trescentrosG_VNA_SH(isorp, in1, in2, indna, x, y,     &
     &                                cost, eps, bccax, rcutoff)
          end if

! Find the charge associated with this shell
            dxn = (Qin(isorp,ialp) - Qneutral(isorp,indna))

! Add this piece for iatom, jatom, and ialp into the total - bcca
            do inu = 1, num_orb(in2)
             do imu = 1, num_orb(in1)
              if (igauss .eq. 1) then
              bccax(imu,inu) = bccax(imu,inu) +                         &
     &            smatG(imu,inu,mneigh,iatom)/R_na(isorp,indna)
              end if

              bcca(imu,inu) = bcca(imu,inu) + bccax(imu,inu)*dxn
               if (Kscf .eq. 1) then
             if (iqout .eq. 6) then
             gvhxc(imu,inu,isorp,ialp,mneigh,iatom) = &
             &     gvhxc(imu,inu,isorp,ialp,mneigh,iatom) +    &
             &       (stn1(imu,inu)*bccax(imu,inu) +    &
             &       stn2(imu,inu)*emnpl_noq(imu,inu))*eq2
             ! write(*,*) 'in 3c g = ',
             ! gvhxc(imu,inu,isorp,ialp,mneigh,iatom)
             ! symmetry
             gvhxc(inu,imu,isorp,ialp,jneigh,jatom) = &
             &     gvhxc(imu,inu,isorp,ialp,mneigh,iatom)
             end if ! end if iqout .eq. 6
             end if ! end if Kscf .eq. 1
             end do
            end do

           end do ! end do of isorp loop
 
! Now add bcca to vca, including the smoothing term emnpl.
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
!$omp atomic 
             vca(imu,inu,mneigh,iatom) = vca(imu,inu,mneigh,iatom) +    &
                                      (stn1(imu,inu)*bcca(imu,inu) +    &
     &                                 stn2(imu,inu)*emnpl(imu,inu))*eq2

             !Symmetrize Hamiltonian (April 2018): jneigh is the
            !back_neigh:
              vca(inu,imu,jneigh,jatom) = vca(imu,inu,mneigh,iatom)

            end do
           end do
 
! ****************************************************************************
! End loop over ialp and its common neighbors.
          end if
         end do
        end do
 
        deallocate (smG)
        deallocate (spmG)
        deallocate (smatG)
        deallocate (spmatG)
! Format Statements
! ===========================================================================
 
        return
        end
