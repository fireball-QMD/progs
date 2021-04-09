! copyright info:
!
!                             @Copyright 2006
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

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

 
! assemble_ca_2c.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the two-center and degenerate
! two-center interactions for DOGS.
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
!
! Program Declaration
! ===========================================================================
        subroutine assemble_ca_2c (nprocs, iforce, iordern)
        use charges
        use configuration
        use constants_fireball
        use dimensions
        use forces
        use interactions
        use neighbor_map
        use options, only : iqout
        use scf, only : Kscf
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iforce
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
        real dterm_1
        real dterm_2
        real dstn_temp
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
        real, dimension (3, numorb_max, numorb_max) :: bccapx
        real, dimension (numorb_max, numorb_max) :: bccax
        real, dimension (3, 3, 3) :: deps
        real, dimension (numorb_max, numorb_max) :: dipx
        real, dimension (3, numorb_max, numorb_max) :: dippx
        real, dimension (numorb_max, numorb_max) :: emnpl
        real, dimension (numorb_max, numorb_max) :: emnpl_noq
        real, dimension (3, 3) :: eps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: sighat
        real, dimension (numorb_max, numorb_max) :: stn1
        real, dimension (numorb_max, numorb_max) :: stn2

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ===========================================================================
! Initialize interactions to zero.
        vca = 0.0d0
        ewaldsr = 0.0d0
        dip = 0.0d0
        dipp = 0.0d0
        if (Kscf .eq. 1 .and. iqout .eq. 6) then
          gvhxc = 0.0d0
        end if ! end if Kscf .eq. 1

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
!!$omp&            private (isorp, jatom, jcount, jcount_sav, kforce, matom)  &
!!$omp&            private (mbeta, dq1, dq2, dterm_1, dterm_2, dstn_temp)     &
!!$omp&            private (dxn, rcutoff_j, rend, rend1, rend2, sterm_1)      &
!!$omp&            private (sterm_2, stn_temp1, stn_temp2, y, bcca, bccapx)   &
!!$omp&            private (bccax, deps, eps, dipx, dippx, emnpl, r1, r2, r21)&
!!$omp&            private (sighat, stn1, stn2)
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

! Find charge on iatom
         dq1 = 0.0d0
         do issh = 1, nssh(in1)
          dq1 = dq1 + (Qin(issh,iatom) - Qneutral(issh,in1))
         end do
 
! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)         ! <==== loop 2 over i's neighbors
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)

          rcutoff_j = 0
          do imu = 1, nssh(in2)
           if (rcutoff(in2,imu) .gt. rcutoff_j) rcutoff_j = rcutoff(in2,imu)
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
! CALL DOSCENTROS AND GET DIP
! ****************************************************************************
          isorp = 0
          interaction = 9
          in3 = in2
          call doscentros (interaction, isorp, iforce, in1, in2, in3, y,     &
     &                     eps, deps, dipx, dippx)
 
          do inu = 1, num_orb(in2)
           do imu = 1, num_orb(in1)
            dip(imu,inu,ineigh,iatom) = dipx(imu,inu)
            if (iforce .eq. 1) dipp(:,imu,inu,ineigh,iatom) = dippx(:,imu,inu)
           end do
          end do
 
! ****************************************************************************
!
! ASSEMBLE EWALDSR AND EMNPL FOR ATM CASE
! ****************************************************************************
! First, calculate the interaction ewaldsr - This is the correction due to what
! is included (more accurately) in the short range terms.
          dq2 = 0.0d0
          do issh = 1, nssh(in2)
           dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
          end do
 
! This is special case 1 - <1mu|V(2)|1nu> so the correction is
! (s/2)*delta_q*(1/d12 + 1/d12) - thus the factor of two is canceled.
! Note: these loops are both over in1 indeed!  This is because the two
! wavefunctions are located on the same site, but the potential is located
! at a different site.
! Initialize
          stn1 = 1.0d0
          stn2 = 0.0d0
          emnpl = 0.0d0
 
! Second, calculate the long-range effective monopole.  This term is included
! so that we obtain no discontinuities when atoms leave or enter the 2*rc
! range criteria.  Therefore, "close" three-center interactions are exact,
! while more distant three-center integrals go to effective monopoles.  The
! monopoles are effective in the sense that the two atoms in the matrix
! element, each has a different charge.  Since they are separated, this gives a
! monopole + dipole contribution at long range.
! The smoothing function is smoother(r,rbegin,rend).  We define our answer as
! smoother(r)*exact + (1 - smoother(r))*longrange.  The distance r is the
! distance of the third center from the "effective" center of the bondcharge.
! The effective center of the bondcharge is (d + rc1 - rc2)/2 from r1 except in
! weird cases (see below).
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
! Same center, so use tighter smoother.
             stn_temp2 = 1.0d0 - stn_temp1
             do inu = 1, lssh(issh,in1)*2 + 1
              icount = icount_sav + inu
              do imu = 1, lssh(jssh,in1)*2 + 1
               jcount = jcount_sav + imu
               stn1(icount,jcount) = stn_temp1
               stn2(icount,jcount) = stn_temp2
              end do
             end do
             jcount_sav = jcount
            end do
            icount_sav = icount
           end do

           do inu = 1, num_orb(in1)
            do imu = 1, num_orb(in1)
             emnpl(imu,inu) = (s_mat(imu,inu,matom,iatom)/y)*dq2
             emnpl_noq(imu,inu) = (s_mat(imu,inu,matom,iatom)/y)
!$omp atomic
             ewaldsr(imu,inu,matom,iatom) =                                  &
     &        ewaldsr(imu,inu,matom,iatom) + emnpl(imu,inu)*eq2
             ! Add the ewaldsr to the gvhxc
        if (Kscf .eq. 1 .and. iqout .eq. 6) then
          do issh = 1, nssh(in2)
             gvhxc(imu,imu,issh,jatom,matom,iatom) = &
      &    gvhxc(imu,imu,issh,jatom,matom,iatom) - emnpl_noq(imu,inu)*eq2
          end do ! end- do issh
        end if ! end if Kscf .eq. 1 .and. iqout .eq. 6               


            end do
           end do
 
! End if y .gt. 1.0d-4
          end if
 
! ****************************************************************************
!
! CALL DOSCENTROS AND GET VNA FOR ATOM CASE
! ****************************************************************************
! The vna 2 centers are: ontop (L), ontop (R), and atm.
! First, do vna_atom case. Here we compute <i | v(j) | i> matrix elements.
! We need jatom because the interaction will include both neutral atom and
! isorp pieces. The isorp pieces will invole Qin. Here is a snippet from
! doscenatm:
!         scam(i,j)=scam(i,j)+temp(i,j)*(Qin(isorp,jk)-Qneutral(isorp,jk)
 
! If r1 = r2, then this is a case where the two wavefunctions are at the
! same site, but the potential vna is at a different site (atm case).
 
! Initialize bcca to zero for the charge atom interactions.
          bcca = 0.0d0
 
          kforce = 0
          interaction = 4
          in3 = in1
          do isorp = 1, nssh(in2)
           call doscentros (interaction, isorp, kforce, in1, in2, in3, y,    &
     &                      eps, deps, bccax, bccapx)
           dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
 
! Now correct bccax by doing the stinky correction.
! Note: these loops are both over in1 indeed!  This is because the two
! wavefunctions are located on the same site, but the potential is located
! at a different site.
           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
             bcca(imu,inu) = bcca(imu,inu) + bccax(imu,inu)*dxn
              if (Kscf .eq. 1) then
              if (iqout .eq. 6) then
              gvhxc(imu,inu,isorp,jatom,matom,iatom) = &
            & gvhxc(imu,inu,isorp,jatom,matom,iatom) + &
            & (stn1(imu,inu)*bccax(imu,inu) + stn2(imu,inu)*emnpl_noq(imu,inu))*eq2
              !write(*,*) 'Hello there'
              end if ! end if iqout .eq. 6
              end if ! end if Kscf .eq. 1
            end do
           end do
          end do
 
! Add bcca to vca, including the smoothing.
! Note that the loop below involves num_orb(in1) ONLY. Why?
! Because the potential is somewhere else (or even at iatom), but we are
! computing the vna_atom term, i.e. < phi(i) | v | phi(i) > but V=v(j) )
! interactions.
          do inu = 1, num_orb(in3)
           do imu = 1, num_orb(in1)
!$omp atomic
            vca(imu,inu,matom,iatom) = vca(imu,inu,matom,iatom)              &
     &       + (stn1(imu,inu)*bcca(imu,inu) + stn2(imu,inu)*emnpl(imu,inu))*eq2
           end do
          end do
 
! ****************************************************************************
!
! ASSEMBLE EWALDSR FOR ONTOP CASE
! ****************************************************************************
! This is special case 1 - <1mu|V(1)|2nu> so the correction is
! (s/2)*delta_q*(1/d12).
! Skip self-interaction terms
          if (y .gt. 1.0d-4) then
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             sterm_1 = eq2*s_mat(imu,inu,ineigh,iatom)/(2.0d0*y)
             sterm_2 = eq2*s_mat(imu,inu,ineigh,iatom)/(2.0d0*y)
             dterm_1 = eq2*dip(imu,inu,ineigh,iatom)/(y*y)
             dterm_2 = eq2*dip(imu,inu,ineigh,iatom)/(y*y)
!$omp atomic

             ewaldsr(imu,inu,ineigh,iatom) = ewaldsr(imu,inu,ineigh,iatom)   &
     &        + (dq1*sterm_1 + dq1*dterm_1) + (dq2*sterm_2 - dq2*dterm_2)
           ! ???? ewaldsr on top case??
            if (Kscf .eq. 1 .and. iqout .eq. 6) then
             do issh = 1, nssh(in2)
             emnpl_noq(imu,inu) = sterm_1-dterm_1  ! on top right
               gvhxc(imu,inu,issh,jatom,ineigh,iatom) = &
        &             gvhxc(imu,inu,issh,jatom,ineigh,iatom) - emnpl_noq(imu,inu)
             end do ! end do issh
             do issh = 1, nssh(in1)
             emnpl_noq(imu,inu) = sterm_1+dterm_1   ! on top left
             gvhxc(imu,inu,issh,iatom,ineigh,iatom) = &
      &             gvhxc(imu,inu,issh,iatom,ineigh,iatom) - emnpl_noq(imu,inu)
            end do ! end do issh
            end if ! end Kscf .eq. 1 .and. iqout .eq. 6)
            end do
           end do
          end if
 
! ****************************************************************************
!
! CALL DOSCENTROS AND GET VNA FOR ONTOP CASE
! ****************************************************************************
! if r1 .ne. r2, then this is a case where the potential is located
! at one of the sites of a wavefunction (ontop case).
          if (iatom .eq. jatom .and. mbeta .eq. 0) then
 
! Do nothing here - special case. Interaction already calculated in atm case.
 
          else
 
! Initialize bcca for charged atom interactions.
           bcca = 0.0d0
 
! For the vna_ontopl case, the potential is on the first atom (j):
! Charged atom piece
           interaction = 2
           in3 = in2
           do isorp = 1, nssh(in1)
            call doscentros (interaction, isorp, kforce, in1, in1, in3, y,   &
     &                       eps, deps, bccax, bccapx)
 
            dxn = (Qin(isorp,iatom) - Qneutral(isorp,in1))
            do inu = 1, num_orb(in3)
             do imu = 1, num_orb(in1)
              bcca(imu,inu) = bcca(imu,inu) + dxn*bccax(imu,inu)
               if (Kscf .eq. 1) then
              if (iqout .eq. 6) then
              gvhxc(imu,inu,isorp,iatom,ineigh,iatom) = &
           & gvhxc(imu,inu,isorp,iatom,ineigh,iatom) + &
           &     bccax(imu,inu)*eq2
              end if ! end if iqout .eq. 6
              end if ! end if Kscf .eq. 1
             end do
            end do
           end do
 
! For the vna_ontopr case, the potential is on the second atom (j):
! Charged atom piece
           interaction = 3
           in3 = in2
           do isorp = 1, nssh(in2)
            call doscentros (interaction, isorp, kforce, in1, in2, in3, y,   &
     &                       eps, deps, bccax, bccapx)
 
            dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
            do inu = 1, num_orb(in3)
             do imu = 1, num_orb(in1)
              bcca(imu,inu) = bcca(imu,inu) + dxn*bccax(imu,inu)
              if (Kscf .eq. 1) then
              if (iqout .eq. 6) then
              gvhxc(imu,inu,isorp,jatom,ineigh,iatom) = &
           & gvhxc(imu,inu,isorp,jatom,ineigh,iatom) + &
           &     bccax(imu,inu)*eq2
              end if ! end if iqout .eq. 6
              end if ! end if Kscf .eq. 1
             end do
            end do
           end do
 
! Now put into vca.
           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
!$omp atomic
             vca(imu,inu,ineigh,iatom) =                                     &
     &        vca(imu,inu,ineigh,iatom) + bcca(imu,inu)*eq2
            end do
           end do
 
! End if for r1 .ne. r2 case
          end if
 
! ****************************************************************************
! End loop over iatom and its neighbors - jatom.
         end do ! do ineigh
        end do ! do iatom

! Format Statements
! ===========================================================================
 
        return
        end subroutine assemble_ca_2c
