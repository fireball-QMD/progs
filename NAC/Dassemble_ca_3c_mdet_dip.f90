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

 
! Dassemble_ca_3c_dip.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the true three-center interactions.
! These are true three-center in that iatom .ne. jatom .ne. katom.
! The matrix elements look like <Psi(1)|V(3)|Psi(2)>.
!
! The bulk of the work is done in Dtrescentros.f. This program assembles the
! results.
!
! A third party term is when we consider the NA (or etc.) to be at the origin
! and take the matrix element between a pair of neighbors, neither of which is
! the NA (or etc.), and the pair is truly a pair, and not an atom.
! the results are stored in: f3ca(ix,i), f3xc(ix,i), etc.
!
! The generic form of the matrix element is
! v = <mu(r-(r1+ratm)) ! v(r-ratm) ! nu(r-(r2+ratm))>.
! We want the derivative with respect to ratm keeping the position of
! orbitals mu, nu fixed (const r1+ratm).
!
! The derivative f3 looks like for the neutral atom case:
! f3ca = - sum (all neighbors of atom alpha at (li,bi) but bi.ne.balph)
!    * sum (all neighbors m of (li,bi), and not having b value balph)
!    * rho(mu,nu,i,m)* deriv wrt balpha <i! v(balph) !j>.
! Note the minus sign to make it "force-like".
!
! This program gets fa, fb, and fc pieces from Dtrescentros.
!
! See notes
! "general theory of third party terms" for the definition on p. 2 of
! these three terms. Breifly fa=-d/dratm, fb=-d/d(bi), fc=-d/d(bj)
! where matrix elements are < psi(r-bi) ! v(r-ratm) ! psi(r-bj)>
!
! ===========================================================================
! Code written by:
! James P. Lewis
!
! code rewritten by Jesús I. Mendieta-Moreno & Diego Soler
! for the dipole long-range theory
!
! see equation (11) in PRB 52, 1618 (1995)
! ===========================================================================
 
! Program Declaration
! ===========================================================================
        subroutine Dassemble_ca_3c_mdet_dip (nprocs, iordern, igauss) 
        use charges
        use configuration
        use constants_fireball
        use density
        use dimensions
        use forces
        use interactions
        use neighbor_map
        use nonadiabatic
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: igauss
        integer, intent (in) :: iordern
        integer, intent (in) :: nprocs
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ialp
        integer iatom
        integer iatomstart
        integer ibeta
        integer icount
        integer ierror
        integer icount_sav
        integer imu
        integer in1
        integer in2
        integer indna
        integer ineigh
        integer interaction
        integer inu
        integer isorp
        integer issh
        integer ix
        integer jatom
        integer jbeta
        integer jcount
        integer jcount_sav
        integer jssh
        integer mneigh
        integer my_proc
        integer natomsp
        integer jneigh
 
        real cost
        real distance13
        real distance23
        real dq1
        real dq2
        real dq3
        real dstn_temp1
        real dstn_temp2
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
        real, dimension (3, numorb_max, numorb_max) :: demnplA
        real, dimension (3, numorb_max, numorb_max) :: demnplB
        real, dimension (3, numorb_max, numorb_max) :: demnplC
        real, dimension (3, 3, 3) :: depsA
        real, dimension (3, 3, 3) :: depsB
        real, dimension (3) :: dpterm
! JIMM_JOM
        real dstn1
        real dstn2
        real, dimension (3) :: dstnA
        real, dimension (3) :: dstnB
        real, dimension (3) :: dstnC
!
        real, dimension (numorb_max, numorb_max) :: emnpl
        real, dimension (3, 3) :: eps
        real, dimension (3, numorb_max, numorb_max) :: f3caXa
        real, dimension (3, numorb_max, numorb_max) :: f3caXb
        real, dimension (3, numorb_max, numorb_max) :: f3caXc
        real, dimension (3, numorb_max, numorb_max) :: f3caXa_sorp
        real, dimension (3, numorb_max, numorb_max) :: f3caXb_sorp
        real, dimension (3, numorb_max, numorb_max) :: f3caXc_sorp
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: rhat
        real, dimension (3) :: rhatA1
        real, dimension (3) :: rhatA2
        real, dimension (3) :: rna
        real, dimension (3) :: rnabc
        real, dimension (3) :: sighat
        real, dimension (3) :: spterm
        real, dimension (3) :: ddterm
        real, dimension (3) :: dptermA
        real, dimension (3) :: dptermB
! JIMM_JOM
        real stn1
        real stn2


! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Procedure
! ============================================================================
! Initialize the force contributions to zero.
        f3caa = 0.0d0
        f3cab = 0.0d0
        f3cac = 0.0d0
 
! Determine which atoms are assigned to this processor.
        if (iordern .eq. 1) then
         call MPI_COMM_RANK (MPI_BTN_WORLD, my_proc, ierror)
         natomsp = natoms/nprocs
         if (my_proc .lt. mod(natoms,nprocs)) then
          natomsp = natomsp + 1
          iatomstart = natomsp*my_proc + 1
         else
          iatomstart = (natomsp + 1)*mod(natoms,nprocs)                      &
                      + natomsp*(my_proc - mod(natoms,nprocs)) + 1
         end if
        else
         iatomstart = 1
         natomsp = natoms
        end if

! Choose atom ialp in the central cell. This is the atom whose position
! we take the derivative, and is the atom who has the the neutral atom
! potential.
! Loop over the atoms in the central cell.
!$omp parallel do private (rna, indna, rcutoff_ialp, imu, ineigh, mneigh)     &
!$omp  private (iatom, ibeta, r1, in1, dq1, issh, jatom, jbeta, r2, in2, dq2) &
!$omp  private (r21, y, sighat, rnabc, x, rhat, cost, eps, depsA, depsB)      &
!$omp  private (distance13, distance23, rhatA1, rhatA2, dq3, icount_sav)      &
!$omp  private (jcount_sav, jssh, rend1, rend2, stn_temp1)                    &
!$omp  private (dstn_temp1, stn_temp2, dstn_temp2, inu, icount, jcount, stn1) &
!$omp  private (stn2, dstn1, dstn2, dstnB, dstnC, dstnA, sterm, dterm, spterm)&
!$omp  private (dpterm, emnpl, demnplA, demnplB, demnplC, bcca, f3caXa)       &
!$omp  private (f3caXb, f3caXc, interaction, isorp, f3caXa_sorp, f3caXb_sorp) &
!$omp  private (f3caXc_sorp, bccax, dxn)
        do ialp = iatomstart, iatomstart - 1 + natomsp
         rna(:) = ratom(:,ialp)
         indna = imass(ialp)
         rcutoff_ialp = 0.0d0
         do imu = 1, nssh(indna)
          if (rcutoff(indna,imu) .gt. rcutoff_ialp)                          & 
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
 
! Find the net charge on iatom
           dq1 = 0.0d0
           do issh = 1, nssh(in1)
            dq1 = dq1 + (Qin(issh,iatom) - Qneutral(issh,in1))
           end do
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
! Find the net charge on iatom
           dq2 = 0.0d0
           do issh = 1, nssh(in2)
            dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
           end do
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
           if (x .lt. 1.0d-03) then
            rhat(1) = 0.0d0
            rhat(2) = 0.0d0
            rhat(3) = 0.0d0
           else
            rhat(:) = rnabc(:)/x
           end if
           cost = sighat(1)*rhat(1) + sighat(2)*rhat(2) + sighat(3)*rhat(3)
 
           call epsilon (rhat, sighat, eps)
 
! dera3 = depsA = deps/dratm in the 3-center system
! der13 = dpesB = deps/dr1 in the 3-center system
           call deps3center (r1, r2, r21, y, rna, rnabc, eps, depsA, depsB)
 
! For now we just do the neutral atom interactions.
! Charged atom interactions are assembled in assemble_ca_3c.f
! So set isorp = 0 within this subroutine.
!
!              interaction    subtypes     index
!
!      bcna         1           0..9(max)   1..10
!      xc3c         2           0..6        11..17
 
! ****************************************************************************
!
! PERFORM ACTUAL CALCULATIONS
! SET-UP AND ASSEMBLE EWALDSR AND DEMNPL
! ****************************************************************************
! Need direction cosines from atom 1 to ratm (rhatA1),
! and atom 2 to ratm (rhatA2).
           distance13 = sqrt((rna(1) - r1(1))**2 + (rna(2) - r1(2))**2       &
     &                                           + (rna(3) - r1(3))**2)
           distance23 = sqrt((rna(1) - r2(1))**2 + (rna(2) - r2(2))**2       &
     &                                           + (rna(3) - r2(3))**2)
 
! Find the unit vector in rna-1 direction.
           if (distance13 .gt. 1.0d-05) then
            rhatA1(:) = (rna(:) - r1(:))/distance13
           else
            write (*,*) ' distance13 is too small in Dassemble_ca_3c_dip.f '
            write (*,*) ' This can not be so!!!! '
           end if
 
! Find the unit vector in rna-2 direction.
           if (distance23 .gt. 1.0d-05) then
            rhatA2(:) = (rna(:) - r2(:))/distance23
           else
            write (*,*) ' distance23 is too small in Dassemble_ca_3c_dip.f '
            write (*,*) ' This can not be so!!!! '
           end if
 
! Now let's calculate the asymptotic form so that we can match these
! better with the sticky smooters.
! Reminder: <B|A|C>= stn1*exact+(1-stn1)*asympt
! So d/dr = stn1 * dexact + dstn1*exact + (1-stn1)*dasympt + d(1-stn1)*asympt.
!               A              B               C                 D.
! Calculate d(asymptote).
! Note that we are using the effective dipole theory here.
!
! First, calculate the interaction ewaldsr - This is the correction due to what
! is included (more accurately) in the short range terms.
! Note that dsrewald is equivalent to demnpl calculated below.
           dq3 = 0.0d0
           do issh = 1, nssh(indna)
            dq3 = dq3 + (Qin(issh,ialp) - Qneutral(issh,indna))
           end do
 
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
!
! JIMM-JOM: we simplify stn-stuff
             rend1 = rcutoff_i + rcutoff_ialp
             rend2 = rcutoff_j + rcutoff_ialp
             call smoother (distance13, rend1, smt_elect, stn_temp1, dstn1)
             call smoother (distance23, rend2, smt_elect, stn_temp2, dstn2)
             stn1 = stn_temp1*stn_temp2
             stn2 = 1.0d0 - stn1
! Different centers, multiply smoothers (this way it goes to zero for the
! the smaller smoother).  Can't just pick smaller rc based one (like
! two-center case), since distances also vary.
 
! d/drB (stn1m*stn2)
!  = dstn1 * stn2 * (-eta(from 1 to Atom)) etc. <B|A|C>=<1|atm|2>
! Note:  d/dr(stn) = dstn1*stn2, and not dstn1*stn2+dstn2*stn1, because
! only stn2 depends on B to C, while stn1 depends only on A to C.
           dstnB = 0.0d0
           dstnC = 0.0d0
             dstnB(:) = - dstn1*stn2*rhatA1(:)
             dstnC(:) = - stn1*dstn2*rhatA2(:)
 
! Get dstnA from Newton's 3rd Law.
           dstnA = - dstnB - dstnC

! New long-range dipole theory
!
         if (x .lt. 1.0d-05) then
                 emnpl = 0.0d0
                 demnplA = 0.0d0
                 demnplB = 0.0d0
                 demnplC = 0.0d0

         else

           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)

             sterm = s_mat(imu,inu,mneigh,iatom)

             dterm = (dipc(1,imu,inu,mneigh,iatom)*rnabc(1)   &
                &   + dipc(2,imu,inu,mneigh,iatom)*rnabc(2)   &
                &   + dipc(3,imu,inu,mneigh,iatom)*rnabc(3))

             ddterm(:) = (dippc(:,1,imu,inu,mneigh,iatom)*rnabc(1)   &
               &        + dippc(:,2,imu,inu,mneigh,iatom)*rnabc(2)   &
               &        + dippc(:,3,imu,inu,mneigh,iatom)*rnabc(3))

             spterm(:) = sp_mat(:,imu,inu,mneigh,iatom)

             dptermA(:)=   dipc(:,imu,inu,mneigh,iatom)/(x*x*x)   &
               &         - 3*dterm*rnabc(:)/(x*x*x*x*x)

             dptermB(:)= - 0.50*dipc(:,imu,inu,mneigh,iatom)/(x*x*x)  & 
               &         + 0.50*3*dterm*rnabc(:)/(x*x*x*x*x)          &
               &         + ddterm(:)/(x*x*x)

             emnpl(imu,inu) =  dq3*sterm/x   &
               &             + dq3*dterm/(x*x*x)

             demnplA(:,imu,inu) = - dq3*sterm*rnabc(:)/(x*x*x)   &
               &                  + dq3*dptermA(:)
           
             demnplB(:,imu,inu) =  dq3*0.50*sterm*rnabc(:)/(x*x*x)  &
               &                 + dq3*spterm(:)/(x)                &
               &                 + dq3*dptermB(:)

             demnplC(:,imu,inu) = - demnplA(:,imu,inu) - demnplB(:,imu,inu)
            end do
           end do
         end if
 
! ****************************************************************************
!
! PERFORM ACTUAL CALCULATIONS
! CALL TRESCENTROS FOR NEUTRAL ATOM PIECE
! ****************************************************************************
! Initialize bcca, f3caXa, f3caXb, and f3caXc
           bcca = 0.0d0
           f3caXa = 0.0d0
           f3caXb = 0.0d0
           f3caXc = 0.0d0

           interaction = 1
           do isorp = 1, nssh(indna)
            if (igauss .eq. 0) then
            call Dtrescentros (interaction, isorp, isorpmax, in1,       &
     &                         in2, indna, x, y, cost, eps, depsA,      &
     &                         depsB, rhat, sighat, bccax, f3caXa_sorp, &
     &                         f3caXb_sorp, f3caXc_sorp, nspecies)
            else
            call DtrescentrosG_VNA_SH (isorp, in1, in2, indna, x, y,    &
     &                                 cost, eps, depsA, depsB, rhat,   &
     &                                 sighat, bccax, f3caXa_sorp,      &
     &                                 f3caXb_sorp, f3caXc_sorp,        &
     &                                 rcutoff)
            end if
 
! Find the charge associated with this shell
            dxn = (Qin(isorp,ialp) - Qneutral(isorp,indna))
 
! Add this piece for iatom, jatom, and ialp into the total - bcca and f3caXa
            do inu = 1, num_orb(in2)
             do imu = 1, num_orb(in1)
              bcca(imu,inu) = bcca(imu,inu) + bccax(imu,inu)*dxn
              f3caXa(:,imu,inu) = f3caXa(:,imu,inu) + f3caXa_sorp(:,imu,inu)*dxn
              f3caXb(:,imu,inu) = f3caXb(:,imu,inu) + f3caXb_sorp(:,imu,inu)*dxn
              f3caXc(:,imu,inu) = f3caXc(:,imu,inu) + f3caXc_sorp(:,imu,inu)*dxn
             end do
            end do
           end do

! Combine the charge contribution with the smoothing pieces - demnpl
! Assemble f3caa, f3cab, and f3cac.
! Note that the - signs below make the contributions force-like.
! f3caXa, f3caXb, and f3caXc are already force-like from Dtrescentros.f
!$omp critical (Dca3)
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             do ix = 1, 3
              f3caa(ix,ialp) = f3caa(ix,ialp)                                &
     &         + 2.0d0*rho(imu,inu,mneigh,iatom)*eq2                               &
     &           *(stn1*f3caXa(ix,imu,inu)                          &
     &             - dstnA(ix)*bcca(imu,inu)                         &
     &             - stn2*demnplA(ix,imu,inu)                       &
     &             + dstnA(ix)*emnpl(imu,inu))
              f3cab(ix,iatom) = f3cab(ix,iatom)                              &
     &         + 2.0d0*rho(imu,inu,mneigh,iatom)*eq2                               &
     &           *(stn1*f3caXb(ix,imu,inu)                          &
     &             - dstnB(ix)*bcca(imu,inu)                         &
     &             - stn2*demnplB(ix,imu,inu)                       &
     &             + dstnB(ix)*emnpl(imu,inu))
              f3cac(ix,jatom) = f3cac(ix,jatom)                              &
     &         + 2.0d0*rho(imu,inu,mneigh,iatom)*eq2                               &
     &           *(stn1*f3caXc(ix,imu,inu)                          &
     &             - dstnC(ix)*bcca(imu,inu)                         &
     &             - stn2*demnplC(ix,imu,inu)                       &
     &             + dstnC(ix)*emnpl(imu,inu))   
             end do
            end do
           end do   
!$omp end critical (Dca3)                 
! Add in the dewaldsr term to f3caa, f3cab, and f3cac. This is
! accomplished by just adding in demnpl*rho again. This compensates the
! addition of the long-range terms in Dassemble_lr.f90, which shouldn't
! be included when the atoms are neighbors.
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             do ix = 1, 3
              f3caa(ix,ialp) = f3caa(ix,ialp)                                &
     &         + 2.0d0*rho(imu,inu,mneigh,iatom)*demnplA(ix,imu,inu)*eq2
              f3cab(ix,iatom) = f3cab(ix,iatom)                              &
     &         + 2.0d0*rho(imu,inu,mneigh,iatom)*demnplB(ix,imu,inu)*eq2
              f3cac(ix,jatom) = f3cac(ix,jatom)                              &
     &         + 2.0d0*rho(imu,inu,mneigh,iatom)*demnplC(ix,imu,inu)*eq2
             end do ! do ix

            end do ! do imu
           end do ! do inu

 
! JOM : gh_3c_ca
! Notice the minus sign, since f3caa etc are force-like
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             do ix = 1, 3
              gh_3c(ix,ialp,imu,inu,mneigh,iatom) =                     &
     &        gh_3c(ix,ialp,imu,inu,mneigh,iatom) - eq2                 &
     &         *(stn1*f3caXa(ix,imu,inu)                         &
     &             - dstnA(ix)*bcca(imu,inu)                         &
     &             - stn2*demnplA(ix,imu,inu)                       &
     &             + dstnA(ix)*emnpl(imu,inu))

              gh_3c(ix,ialp,inu,imu,jneigh,jatom) =                     &
     &        gh_3c(ix,ialp,imu,inu,mneigh,iatom)

              gh_3c(ix,iatom,imu,inu,mneigh,iatom) =                     &
     &        gh_3c(ix,iatom,imu,inu,mneigh,iatom) - eq2                 &
     &           *(stn1*f3caXb(ix,imu,inu)                          &
     &             - dstnB(ix)*bcca(imu,inu)                         &
     &             - stn2*demnplB(ix,imu,inu)                       &
     &             + dstnB(ix)*emnpl(imu,inu))

              gh_3c(ix,iatom,inu,imu,jneigh,jatom) =                     &
     &        gh_3c(ix,iatom,imu,inu,mneigh,iatom) 

              gh_3c(ix,jatom,imu,inu,mneigh,iatom) =                     &
     &        gh_3c(ix,jatom,imu,inu,mneigh,iatom) - eq2                 &
     &           *(stn1*f3caXc(ix,imu,inu)                          &
     &             - dstnC(ix)*bcca(imu,inu)                         &
     &             - stn2*demnplC(ix,imu,inu)                       &
     &             + dstnC(ix)*emnpl(imu,inu))

              gh_3c(ix,jatom,inu,imu,jneigh,jatom) =                     &
     &        gh_3c(ix,jatom,imu,inu,mneigh,iatom)

             end do ! do ix
            end do ! do imu
           end do ! do inu

! Add in the dewaldsr term to gh_3c_ca
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             do ix = 1, 3
              gh_3c(ix,ialp,imu,inu,mneigh,iatom) =                     &
     &        gh_3c(ix,ialp,imu,inu,mneigh,iatom) - eq2                 &
     &         *demnplA(ix,imu,inu)


              gh_3c(ix,ialp,inu,imu,jneigh,jatom) =                     &
     &        gh_3c(ix,ialp,imu,inu,mneigh,iatom) 


              gh_3c(ix,iatom,imu,inu,mneigh,iatom) =                     &
     &        gh_3c(ix,iatom,imu,inu,mneigh,iatom) - eq2                 &
     &         *demnplB(ix,imu,inu)


              gh_3c(ix,iatom,inu,imu,jneigh,jatom) =                     &
     &        gh_3c(ix,iatom,imu,inu,mneigh,iatom)

              gh_3c(ix,jatom,imu,inu,mneigh,iatom) =                     &
     &        gh_3c(ix,jatom,imu,inu,mneigh,iatom) - eq2                 &
     &         *demnplC(ix,imu,inu)


              gh_3c(ix,jatom,inu,imu,jneigh,jatom) =                     &
     &        gh_3c(ix,jatom,imu,inu,mneigh,iatom)

             end do ! do ix
            end do ! do imu
           end do ! do inu
 
! ****************************************************************************
! End loop over ialp and its common neighbors.
          end if
         end do
        end do
 
! Format Statements
! ===========================================================================
 
        return
        end
