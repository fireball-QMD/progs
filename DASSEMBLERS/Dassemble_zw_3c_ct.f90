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
 
! Program Declaration
! ===========================================================================
        subroutine Dassemble_zw_3c_ct (nprocs, iordern, igauss) 
        use charges
        use configuration
        use constants_fireball
        use density
        use dimensions
        use forces
        use interactions
        use neighbor_map
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
        integer :: issh1,issh2, ineigh1, ineigh2
 
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
        real :: A,B

        real, dimension (numorb_max, numorb_max) :: bcca
        real, dimension (numorb_max, numorb_max) :: bccax
        real, dimension (3) :: dAb,dBb,dAc,dBc

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
! JIMM_JOM
        real stn1
        real stn2


! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Procedure
! ============================================================================
! Initialize the force contributions to zero.
        f3xca_ca = 0.0d0
        f3xcb_ca = 0.0d0
        f3xcc_ca = 0.0d0
 
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
 
! Loop over the neighbors of each ialp.
! Now look at all common neighbor pairs which were figured out in main.
! The common neighbors were figured out in common_neighbors.f90
         do ineigh = 1, neigh_comn(ialp)
          mneigh = neigh_comm(ineigh,ialp)
 
! The second atom (jatom) is the mneigh'th neighbor of iatom.
          if (mneigh .ne. 0) then
           iatom = neigh_comj(1,ineigh,ialp)
           ibeta = neigh_comb(1,ineigh,ialp)
           ineigh1 = neigh_com_ng(1,ineigh,ialp)
           r1(:) = ratom(:,iatom) + xl(:,ibeta)
           in1 = imass(iatom)
 
! Find the net charge on iatom
           dq1 = 0.0d0
           do issh = 1, nssh(in1)
            dq1 = dq1 + (Qin(issh,iatom) - Qneutral(issh,in1))
           end do
 
           jatom = neigh_comj(2,ineigh,ialp)
           jbeta = neigh_comb(2,ineigh,ialp)
           ineigh2 = neigh_com_ng(2,ineigh,ialp)
           r2(:) = ratom(:,jatom) + xl(:,jbeta)
           in2 = imass(jatom)

 
! Find the net charge on iatom
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
!           call deps3center (r1, r2, r21, y, rna, rnabc, eps, depsA, depsB)
 
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
     !           call Dtrescentros (interaction, isorp, isorpmax, in1,       &
     ! &                         in2, indna, x, y, cost, eps, depsA,      &
     ! &                         depsB, rhat, sighat, bccax, f3caXa_sorp, &
     ! &                         f3caXb_sorp, f3caXc_sorp, nspecies)
 
! Find the charge associated with this shell
            dxn = (Qin(isorp,ialp) - Qneutral(isorp,indna))

!   <B | A | C>, 
! Add this piece for iatom, jatom, and ialp into the total - bcca and f3caXa
            do inu = 1, num_orb(in2)
             do imu = 1, num_orb(in1)
               issh1=orb2shell(imu,in1)
               issh2=orb2shell(inu,in2)
               A = 0.5*s_mat(imu,inu,mneigh,iatom)-dip(imu,inu,mneigh,iatom)/y   
               B = 0.5*s_mat(imu,inu,mneigh,iatom)+dip(imu,inu,mneigh,iatom)/y  
               dAb(:) = 0.5*sp_mat(:,imu,inu,mneigh,iatom)-dipp(:,imu,inu,mneigh,iatom)/y- &
                     & dip(imu,inu,mneigh,iatom)*sighat(:)/(y*y)
               dBb(:) = 0.5*sp_mat(:,imu,inu,mneigh,iatom)+dipp(:,imu,inu,mneigh,iatom)/y+ &
                     & dip(imu,inu,mneigh,iatom)*sighat(:)/(y*y)
               dAc(:) = -dAb(:)   !!?
               dBc(:) = -dBb(:)  
               !force:   !Careful with the signs!!!
                f3caXc(:,imu,inu) = f3caXc(:,imu,inu)+ &
              & dxn*(-B*g2nup(:,isorp,issh2,ineigh2,ialp)+ &
              &      dAc(:)*g2nu(isorp,issh1,ineigh1,ialp)+dBc(:)*g2nu(isorp,issh2,ineigh2,ialp)   )

                f3caXb(:,imu,inu) = f3caXb(:,imu,inu)+&
              & dxn*(-A*g2nup(:,isorp,issh1,ineigh1,ialp)+ &
              &      dAb(:)*g2nu(isorp,issh1,ineigh1,ialp)+dBb(:)*g2nu(isorp,issh2,ineigh2,ialp))
 
                f3caXa(:,imu,inu) = f3caXa(:,imu,inu)+&
              & dxn*(A*g2nup(:,isorp,issh1,ineigh1,ialp)+B*g2nup(:,isorp,issh2,ineigh2,ialp))



               !OJO..!! One variable per atom..!!

              ! f3caXa(:,imu,inu) = f3caXa(:,imu,inu) + f3caXa_sorp(:,imu,inu)*dxn
              ! f3caXb(:,imu,inu) = f3caXb(:,imu,inu) + f3caXb_sorp(:,imu,inu)*dxn
              ! f3caXc(:,imu,inu) = f3caXc(:,imu,inu) + f3caXc_sorp(:,imu,inu)*dxn

             end do !end do imu
            end do !end do inu
           end do !end do isorp
  
             do inu = 1,num_orb(in2)
              do imu = 1,num_orb(in1)
                 do ix = 1, 3
                   f3xca_ca(ix,ialp) = f3xca_ca(ix,ialp) &
                   &         - 2*rho(imu,inu,mneigh,iatom) &
                   &           *f3caXa(ix,imu,inu)    
                  f3xcb_ca(ix,iatom) = f3xcb_ca(ix,iatom) &
                   &         - 2*rho(imu,inu,mneigh,iatom) &
                   &           *f3caXb(ix,imu,inu)                   
                  f3xcc_ca(ix,jatom) = f3xcc_ca(ix,jatom) &
                   &         - 2*rho(imu,inu,mneigh,iatom) &
                   &           *f3caXc(ix,imu,inu)                  
                 end do !end do ix
               end do !end do imu
             end do !end do inu

 
! ****************************************************************************
! End loop over ialp and its common neighbors.
          end if !end if mneigh = 0
         end do !end do ineigh
        end do !end do ialp
        !f3xca_ca = 0.0d0
        !f3xcb_ca = 0.0d0
        !f3xcc_ca = 0.0d0
 
! Format Statements
! ===========================================================================
 
        return
        end
