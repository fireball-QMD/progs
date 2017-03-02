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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
! Dassemble_ca_2c_mdet.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the two-center and degenerate
! two-center DOGS (charge transfer) interactions. 
!
! JOM : added calculation of non-adiabatic couplings
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
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_ca_2c_mdet (nprocs, iordern)
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
        integer, intent (in) :: iordern
        integer, intent (in) :: nprocs
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iatomstart
        integer icount
        integer icount_sav
        integer ierror
        integer imu
        integer in1
        integer in2
        integer in3
        integer ineigh
        integer interaction
        integer inu
        integer isorp
        integer issh
        integer ix
        integer jatom
        integer jcount
        integer jcount_sav
        integer jssh
        integer kforce
        integer matom
        integer mbeta
        integer my_proc
        integer natomsp
 
        real dq1
        real dq2
        real dstn_temp
        real dterm_1
        real dterm_2
        real dxn
        real rcutoff_j
        real rend
        real rend1
        real rend2
        real sterm_1
        real sterm_2
        real stn_temp1
        real stn_temp2
        real y
        real rcutoff_i
 
        real, dimension (numorb_max, numorb_max) :: bcca
        real, dimension (3, numorb_max, numorb_max) :: bccap
        real, dimension (3, numorb_max, numorb_max) :: bccapx
        real, dimension (numorb_max, numorb_max) :: bccax
        real, dimension (numorb_max, numorb_max) :: demnpl
        real, dimension (3, 3, 3) :: deps
        real, dimension (3, numorb_max, numorb_max) :: dewaldsr
        real, dimension (3) :: dpterm_1
        real, dimension (3) :: dpterm_2
        real, dimension (numorb_max, numorb_max) :: dstn1
        real, dimension (numorb_max, numorb_max) :: dstn2
        real, dimension (numorb_max, numorb_max) :: emnpl
        real, dimension (3, 3) :: eps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: sighat
        real, dimension (3) :: spterm_1
        real, dimension (3) :: spterm_2
        real, dimension (numorb_max, numorb_max) :: stn1
        real, dimension (numorb_max, numorb_max) :: stn2
 
! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ===========================================================================
! Initialize the forces to zero.
        faca = 0.0d0
        fotca = 0.0d0
! JOM Initialize gh_2c_ca and gh_atm_ca
       ! gh_2c_ca = 0.0d0
       ! gh_atm_ca = 0.0d0

 
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

! Loop over the atoms in the central cell.
!!$omp parallel do private (icount, icount_sav, in1, in2, in3, interaction)   &
!!$omp&            private (jatom, jcount, jcount_sav, kforce, matom, mbeta)  &
!!$omp&            private (dq1, dq2, dstn_temp, dterm_1, dterm_2, dxn)       &
!!$omp&            private (rcutoff_j, rend, rend1, rend2, sterm_1, sterm_2)  &
!!$omp&            private (stn_temp1, stn_temp2, y, bcca, bccap, bccapx)     &
!!$omp&            private (bccax, demnpl, deps, dewaldsr, dpterm_1, dpterm_2)&
!!$omp&            private (dstn1, dstn2, emnpl, eps, r1, r2, r21, sighat)    &
!!$omp&            private (stn1, stn2, spterm_1, spterm_2)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neigh_self(iatom)
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)

! ENRIQUE-JOM: calculate rcutoff_i for smoother
          rcutoff_i = 0
          do imu = 1, nssh(in1)
           if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
          end do
! End ENRIQUE-JOM

! Find charge on iatom.
         dq1 = 0.0
         do issh = 1, nssh(in1)
          dq1 = dq1 + (Qin(issh,iatom) - Qneutral(issh,in1))
         end do
 
! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)    ! <==== loop 2 over iatom's neighbors
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)

          rcutoff_j = 0
          do imu = 1, nssh(in2)
           if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
          end do
 
! Find charge on jatom.
          dq2 = 0.0d0
          do issh = 1, nssh(in2)
           dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
          end do
 
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
          else
           sighat(:) = r21(:)/y
          end if
 
          call epsilon (r2, sighat, eps)
          call deps2cent (r1, r2, eps, deps)
 
! ****************************************************************************
!
! ASSEMBLE DEWALDSR AND DEMNPL FOR ATM CASE
! ****************************************************************************
! Set emnpl, demnpl, and dewaldsr to zero
          stn1 = 1.0d0
          dstn1 = 0.0d0
          stn2 = 0.0d0
          dstn2 = 0.0d0
          emnpl = 0.0d0
          demnpl = 0.0d0
          dewaldsr = 0.0d0
 
! Now that we are doing the atm case, let's be clear and call the second loop
! ratm. See dnanlxc (fireball96) for tips.
 
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
! Skip self-interaction terms
          if (y .gt. 1.0d-4) then
           icount_sav = 0
           do issh = 1, nssh(in1)
            jcount_sav = 0
            do jssh = 1, nssh(in1)
! ENRIQUE-JOM: use rcutoff_i + rcutoff_j as the smoother
!             rend1 = rcutoff(in1,issh) + rcutoff_j
!             rend2 = rcutoff(in1,jssh) + rcutoff_j
!             rend = min(rend1,rend2)
             rend = rcutoff_i + rcutoff_j
             call smoother (y, rend, smt_elect, stn_temp1, dstn_temp)
             stn_temp2 = 1.0d0 - stn_temp1
             do inu = 1, lssh(issh,in1)*2 + 1
              icount = icount_sav + inu
              do imu = 1, lssh(jssh,in1)*2 + 1
               jcount = jcount_sav + imu
               stn1(icount,jcount) = stn_temp1
               stn2(icount,jcount) = stn_temp2
               dstn1(icount,jcount) = dstn_temp
               dstn2(icount,jcount) = - dstn_temp
              end do
             end do
             jcount_sav = jcount
            end do
            icount_sav = icount
           end do

! This is special case 1 - <1mu|V(2)|1nu> so the correction is
! (s/2)*dq*(1/d12 + 1/d12) - thus the factor of two is canceled.
! Note: these loops are both over in1 indeed!  This is because the two
! wavefunctions are located on the same site, but the potential is located
! at a different site.
           do inu = 1, num_orb(in1)
            do imu = 1, num_orb(in1)
             emnpl(imu,inu) = (s_mat(imu,inu,matom,iatom)/y)*dq2
             demnpl(imu,inu) = - (s_mat(imu,inu,matom,iatom)/(y*y))*dq2
             dewaldsr(:,imu,inu) = - demnpl(imu,inu)*sighat(:)*eq2
            end do
           end do
 
! end if for y .gt. 1.0d-4.
          end if
 
! ****************************************************************************
!
! ASSEMBLE NEUTRAL ATOM FORCE FOR ATOM CASE
! ****************************************************************************
! The vna 2 centers are: ontop (L), ontop (R), and atm.
! First, do vna_atom case. Here we compute <i | v(j) | i> matrix elements.
! We need jatom because the interaction will include both neutral atom and
! isorp pieces. The isorp pieces will invole xn. Here is a snippet from
! doscenatm:
!       scam(i,j)=scam(i,j)+temp(i,j)*(Qin(isorp,jk)-Qneutral(isorp,jk)
 
! If r1 = r2, then this is a case where the two wavefunctions are at the
! same site, but the potential vna is at a different site (atm case).
! The derivative wrt the "atom r1" position (not the NA position) are
! stored in bccapx.
! Form the "force-like" derivative of the atom terms for NA,
! or -(d/dr1) <phi(mu,r-r1)!h(r-ratm)!phi(nu,r-r1)>.
 
! Initialize bcca and bccap for the charge atom interactions
          bcca = 0.0d0
          bccap = 0.0d0
 
          kforce = 1
          interaction = 4
          in3 = in1
          do isorp = 1, nssh(in2)
           call doscentros (interaction, isorp, kforce, in1, in2, in3, y,  &
     &                      eps, deps, bccax, bccapx)
           dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
 
! Note that the loop below involves num_orb(in1) ONLY. Why? Because the
! potential is somewhere else (or even at iatom), but we are computing the
! vna_atom term, i.e. < phi(i) | v | phi(i) > but V=v(j) ) interactions.
           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
             bcca(imu,inu) = bcca(imu,inu) + dxn*bccax(imu,inu)
             bccap(:,imu,inu) = bccap(:,imu,inu) + dxn*bccapx(:,imu,inu)
            end do
           end do
          end do
 
! Now correct bccapx by doing the stinky correction.
! Note: these loops are both over in1 indeed!  This is because the two
! wavefunctions are located on the same site, but the potential is located
! at a different site.
! Here we must get demnpl as a VECTOR derivative. Thst's easy since
! d/dr =  -eta d/dd
! As long as epsilon is called with sighat in the second "spot" as
! call epsilon(R1,sighat,spe), then eps(ix,3)=eta(ix).
! 2 type of terms. d/dd or stn, and d/dd of matrix.
          do inu = 1, num_orb(in3)
           do imu = 1, num_orb(in1)
            bccap(:,imu,inu) =                                               &
     &       stn1(imu,inu)*bccap(:,imu,inu)                                  &
     &        - dstn1(imu,inu)*bcca(imu,inu)*sighat(:)                       &
     &        - (stn2(imu,inu)*demnpl(imu,inu)                               &
     &           + dstn2(imu,inu)*emnpl(imu,inu))*sighat(:)
           end do
          end do
 
! Now write bccap to the charged atom force piece.
! Add dewaldsr to the charged force term.  This is added in, since the energy
! is subtracted.
! Notice the explicit - sign, this makes it force like.
          do inu = 1, num_orb(in3)
           do imu = 1, num_orb(in1)
            do ix = 1, 3
!!$omp atomic
             faca(ix,ineigh,iatom) = faca(ix,ineigh,iatom)                   &
     &        - rho(imu,inu,matom,iatom)*bccap(ix,imu,inu)*eq2               &
     &        + rho(imu,inu,matom,iatom)*dewaldsr(ix,imu,inu)
            end do
           end do
          end do
 
! JOM ; add gradient of vna_atom to gh_atm_ca
          do inu = 1, num_orb(in3)
           do imu = 1, num_orb(in1)
!            gh_atm_ca (:,imu,inu,ineigh,iatom) =                             &
!     &      gh_atm_ca (:,imu,inu,ineigh,iatom) +                             &
!     &      bccap(:,imu,inu)*eq2 - dewaldsr(:,imu,inu)               
            gh_atm (:,imu,inu,ineigh,iatom) =                             &
     &      gh_atm (:,imu,inu,ineigh,iatom) +                             &
     &      bccap(:,imu,inu)*eq2 - dewaldsr(:,imu,inu)               
           end do
          end do
!           gh_atm_ca (:,:,:,ineigh,iatom) =                             &
!    &      gh_atm_ca (:,:,:,ineigh,iatom) +                             &
!    &      bccap(:,:,:)*eq2 - dewaldsr(:,:,:)               
 
! ****************************************************************************
!
! ASSEMBLE NEUTRAL ATOM FORCE FOR ONTOP CASE
! ****************************************************************************
! Ontop case.
! We have Ontop Left, and Ontop Right.
! Left is <1|V(1)|2> and Right is <1|V(2)|2>
! Check to make sure we are not doing <1|V(1)|1>. This is done in the atm case.
          if (iatom .eq. jatom .and. mbeta .eq. 0) then
 
! Do nothing here - special case.
 
          else
 
! First the dewaldsr.
! Note that dewaldsr and dewaldlr are used as matrix elements
! in other programs and they DO HAVE THE eq2 factor!
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
 
! Note: S is off-site. Get sp and dipp from main.f
! Subtraction minus "effect" is because we subtract sr-term from total energy.
             sterm_1 = dq1*eq2*s_mat(imu,inu,ineigh,iatom)/(2.0d0*y)
             sterm_2 = dq2*eq2*s_mat(imu,inu,ineigh,iatom)/(2.0d0*y)
 
             dterm_1 = dq1*eq2*dip(imu,inu,ineigh,iatom)/(y*y)
             dterm_2 = dq2*eq2*dip(imu,inu,ineigh,iatom)/(y*y)
 
             spterm_1(:) = dq1*eq2*sp_mat(:,imu,inu,ineigh,iatom)/(2.0d0*y)
             spterm_2(:) = dq2*eq2*sp_mat(:,imu,inu,ineigh,iatom)/(2.0d0*y)
 
             dpterm_1(:) = dq1*eq2*dipp(:,imu,inu,ineigh,iatom)/(y*y)
             dpterm_2(:) = dq2*eq2*dipp(:,imu,inu,ineigh,iatom)/(y*y)
 
! There are three minus signs in all - one belonging to the sighat term,
! one belonging to the 1/y term and a minus subtraction error.
! No minus signs are "force-like". This is not a force like derivative,
! but a regular derivative d/dr1.  The net effect is a minus sign.
! The two is from (dip/dbc)/dip = - 2 (dip/dbc**3) = 2(dip/dbc)/dbc*(-1/dbc)
! Form the ontop force. Use the derivatives, since the derivative is
! with respect to d/d(ratom) when the atom is ontop atom 1.
             dewaldsr(:,imu,inu) =                                           &
     &        + (spterm_1(:) + dpterm_1(:))                                  &
     &        + (sterm_1 + 2.0d0*dterm_1)*sighat(:)/y                        &
     &        + (spterm_2(:) - dpterm_2(:))                                  &
     &        + (sterm_2 - 2.0d0*dterm_2)*sighat(:)/y
            end do
           end do
 
! Now ontop Left.
! Note that we only compute ontop left. That is because we do cross terms.
! If we were to do both ontop left and ontop right, then we would get
! double counting in the forces.
 
! Initialize bccap for the charge atom interactions
           bccap = 0.0d0
 
! Charged atom piece - Left ontop
           kforce = 1
           interaction = 2
           in3 = in2
           do isorp = 1, nssh(in1)
            call doscentros (interaction, isorp, kforce, in1, in1, in3, y,   &
     &                       eps, deps, bccax, bccapx)
            dxn = (Qin(isorp,iatom) - Qneutral(isorp,in1))
 
! Now add d/dr1 (ewaldsr) to bccapx.
            do imu = 1, num_orb(in1)
             do inu = 1, num_orb(in3)
              bccap(:,imu,inu) = bccap(:,imu,inu) + dxn*bccapx(:,imu,inu)
             end do
            end do
           end do
 
! Notice the explicit - sign which makes f force like.
! Add dewaldsr to the charged force term.  The factor of 0.5d0 is due to
! the fact that we have a factor of 2.0d0 later, since we are combining
! left and right contributions.
           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
             do ix = 1, 3
!!$omp atomic
             fotca(ix,ineigh,iatom) = fotca(ix,ineigh,iatom)                 &
     &        - rho(imu,inu,ineigh,iatom)*bccap(ix,imu,inu)*eq2              &
     &        + 0.5d0*rho(imu,inu,ineigh,iatom)*dewaldsr(ix,imu,inu)
             end do
            end do
           end do

! JOM : add gradient of vna_ontop_left to gh_2c_ca
 
           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
!             gh_2c_ca(:,imu,inu,ineigh,iatom) =                             &
!     &       gh_2c_ca(:,imu,inu,ineigh,iatom) + bccap(:,imu,inu)*eq2   
             gh_2c(:,imu,inu,ineigh,iatom) =                             &
     &       gh_2c(:,imu,inu,ineigh,iatom) + bccap(:,imu,inu)*eq2   
            end do
           end do

! JOM : add now gradient of vna_ontop_right to gh_2c_ca
! JOM-info : Now, add gradient of vna_ontop_right to gh_2c
! JOM-info : Here we need the ontop_right case, but for
! forces (see
! above) the contribution from ontop_right is added later in
! assemble_F.f90 by multiplying by 2.0d0

! Now ontop Right.
 
! Initialize bccap for the charge atom interactions
           bccap = 0.0d0
 
! Charged atom piece - Right ontop
           kforce = 1
           interaction = 3
           in3 = in2
           do isorp = 1, nssh(in2)
            call doscentros (interaction, isorp, kforce, in1, in2, in3, y,   &
     &                       eps, deps, bccax, bccapx)
            dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
 
            do imu = 1, num_orb(in1)
             do inu = 1, num_orb(in3)
              bccap(:,imu,inu) = bccap(:,imu,inu) + dxn*bccapx(:,imu,inu)
             end do
            end do
           end do

           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
! merge large arrays VLADA
!             gh_2c_ca(:,imu,inu,ineigh,iatom) =                             &
!     &       gh_2c_ca(:,imu,inu,ineigh,iatom) + bccap(:,imu,inu)*eq2   
             gh_2c(:,imu,inu,ineigh,iatom) =                             &
     &       gh_2c(:,imu,inu,ineigh,iatom) + bccap(:,imu,inu)*eq2   

            end do
           end do

! Finally,  SUBSTRACT d/dr1 (ewaldsr) to gh_2c_ca
           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
! merge large arrays VLADA
!             gh_2c_ca(:,imu,inu,ineigh,iatom) =                             &
!     &       gh_2c_ca(:,imu,inu,ineigh,iatom) - dewaldsr(:,imu,inu)   
             gh_2c(:,imu,inu,ineigh,iatom) =                             &
     &       gh_2c(:,imu,inu,ineigh,iatom) - dewaldsr(:,imu,inu)   

            end do
           end do

          end if

! ****************************************************************************
! End loop over iatom and its neighbors - jatom.
         end do
        end do

! Format Statements
! ===========================================================================
        return
        end
