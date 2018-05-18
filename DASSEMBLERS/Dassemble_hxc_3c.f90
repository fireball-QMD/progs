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
! Brigham Young University - Eduardo Mendez 
! Arizona State University - Otto F. Sankey
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

 
! Dassemble_hxc_3c.f90
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
! the results are stored in: f3na(ix,i), f3xc(ix,i), etc.
!
! The generic form of the matrix element is
! v = <mu(r-(r1+ratm)) ! v(r-ratm) ! nu(r-(r2+ratm))>.
! We want the derivative with respect to ratm keeping the position of
! orbitals mu, nu fixed (const r1+ratm).
!
! The derivative f3 looks like for the neutral atom case:
! f3na = - sum (all neighbors of atom alpha at (li,bi) but bi.ne.balph)
!    * sum (all neighbors m of (li,bi), and not having b value balph)
!    * rho(mu,nu,i,m)* deriv wrt balpha <i! v(balph) !j>.
! Note the minus sign to make it "force-like".
!
! This program gets fa, fb, and fc pieces from Dtrescentros.
!
! See notes "general theory of third party terms" for the definition on p. 2
! of these three terms. Breifly fa=-d/dratm, fb=-d/d(bi), fc=-d/d(bj)
! where matrix elements are < psi(r-bi) ! v(r-ratm) ! psi(r-bj)>
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
!
! Non-local pseudopotential piece by:
! Otto F. Sankey
! Campus Box 1504
! Department of Physics
! Arizona State University
! Tempe, AZ 85287-1504
! (602) 965-4334 (office)      email: otto.sankey@asu.edu
! ===========================================================================
 
! Program Declaration
! ===========================================================================
        subroutine Dassemble_hxc_3c (nprocs, iordern, itheory, igauss) 
        use charges
        use configuration
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
        integer, intent (in) :: itheory
        integer, intent (in) :: nprocs
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ialp
        integer iatom
        integer iatomstart
        integer ibeta
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
        integer ix
        integer jatom
        integer jbeta
        integer mneigh
        integer my_proc
        integer natomsp
 
        real cost
        real dq1
        real dq2
        real dq3
        real x
        real y
 
        real, dimension (numorb_max, numorb_max) :: bcxcx
        real, dimension (numorb_max, numorb_max) :: bcxcx_ca
        real, dimension (6) :: dqfact          ! there are six different isorps
        real, dimension (3, 3, 3) :: depsA
        real, dimension (3, 3, 3) :: depsB
        real, dimension (3, 3) :: eps
        real, dimension (3, numorb_max, numorb_max) :: f3xcXa
        real, dimension (3, numorb_max, numorb_max) :: f3xcXb
        real, dimension (3, numorb_max, numorb_max) :: f3xcXc
        real, dimension (3, numorb_max, numorb_max) :: f3xcXa_ca
        real, dimension (3, numorb_max, numorb_max) :: f3xcXb_ca
        real, dimension (3, numorb_max, numorb_max) :: f3xcXc_ca
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: rhat
        real, dimension (3) :: rna
        real, dimension (3) :: rnabc
        real, dimension (3) :: sighat

        integer itest
 
! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ========================================================================
! Initialize the force contributions to zero.
        f3xca = 0.0d0
        f3xcb = 0.0d0
        f3xcc = 0.0d0
        if (itheory .eq. 1) then
         f3xca_ca = 0.0d0
         f3xcb_ca = 0.0d0
         f3xcc_ca = 0.0d0
        end if

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

! Choose atom ialp in the central cell. This is the atom whose position
! we take the derivative, and is the atom who has the the neutral atom
! potential.
! Loop over the atoms in the central cell.
!$omp parallel do private (iatom, ibeta, in1, in2, indna, interaction, isorp)&
!$omp&            private (jatom, jbeta, mneigh, cost, dq1, dq2, dq3, x, y)  &
!$omp&            private (bcxcx, bcxcx_ca, depsA, depsB, eps, f3xcXa)       &
!$omp&            private (f3xcXb, f3xcXc, f3xcXa_ca, f3xcXb_ca, f3xcXc_ca)  &
!$omp&            private (r1, r2, r21, rhat, rna, rnabc, sighat)
        do ialp = iatomstart, iatomstart - 1 + natomsp
         rna(:) = ratom(:,ialp)
         indna = imass(ialp)

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

           jatom = neigh_comj(2,ineigh,ialp)
           jbeta = neigh_comb(2,ineigh,ialp)
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
           cost = sighat(1)*rhat(1) + sighat(2)*rhat(2) +               &
     &            sighat(3)*rhat(3)
           call epsilon (rhat, sighat, eps)
 
! dera3 = depsA = deps/dratm in the 3-center system
! der13 = depsB = deps/dr1 in the 3-center system
           call deps3center (r1, r2, r21, y, rna, rnabc,eps,depsA,depsB)
 
! Find the charge for the third neighbor.
           dq3 = 0.0d0
           do issh = 1, nssh(indna)
            dq3 = dq3 + (Qin(issh, ialp) - Qneutral(issh, indna))
           end do

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
! CALL TRESCENTROS FOR EXCHANGE-CORRELATION PIECE - NA CASE
! ****************************************************************************
           isorp = 0
           interaction = 2
           if (igauss .eq. 0) then
            call Dtrescentros (interaction, isorp, ideriv_max, in1, in2,&
     &                         indna, x, y, cost, eps, depsA, depsB,    &
     &                         rhat, sighat, bcxcx, f3xcXa, f3xcXb,     &
     &                         f3xcXc, nspecies)

! The terms f3naXa, f3naXb, and f3naXc are already force-like.
!$omp critical (D3hxc_1)
            do inu = 1, num_orb(in2)
             do imu = 1, num_orb(in1)
              do ix = 1, 3
!$omp atomic
               f3xca(ix,ialp) = f3xca(ix,ialp)                          &
     &          + 2*rho(imu,inu,mneigh,iatom)*f3xcXa(ix,imu,inu)
!$omp atomic
               f3xcb(ix,iatom) = f3xcb(ix,iatom)                        &
     &          + 2*rho(imu,inu,mneigh,iatom)*f3xcXb(ix,imu,inu)
!$omp atomic
               f3xcc(ix,jatom) = f3xcc(ix,jatom)                        &
     &          + 2*rho(imu,inu,mneigh,iatom)*f3xcXc(ix,imu,inu)
              end do
             end do
            end do
!$omp end critical (D3hxc_1)
           else 
            call DtrescentrosGHXC_VXC (in1,in2, indna, x, y, cost, eps, &
     &                               depsA, depsB, rhat, sighat, bcxcx, &
     &                               f3xcXa, f3xcXb, f3xcXc, rcutoff)
!$omp critical (D3hxc_2)
            do inu = 1, num_orb(in2)
             do imu = 1, num_orb(in1)
              do ix = 1, 3
!$omp atomic
               f3xca(ix,ialp)=f3xca(ix,ialp)+rho(imu,inu,mneigh,iatom)  &
     &            *2*nuxc_total(imu,inu,mneigh,iatom)*f3xcXa(ix,imu,inu)
!$omp atomic
               f3xcb(ix,iatom)=f3xcb(ix,iatom)+rho(imu,inu,mneigh,iatom)&
     &            *2*nuxc_total(imu,inu,mneigh,iatom)*f3xcXb(ix,imu,inu)
!$omp atomic
               f3xcc(ix,jatom)=f3xcc(ix,jatom)+rho(imu,inu,mneigh,iatom)&
     &            *2*nuxc_total(imu,inu,mneigh,iatom)*f3xcXc(ix,imu,inu)
              end do
             end do
            end do  
!$omp end critical (D3hxc_2)
           end if

! ****************************************************************************
!
! PERFORM ACTUAL CALCULATIONS
! CALL TRESCENTROS FOR EXCHANGE-CORRELATION PIECE - CA CASE
! ****************************************************************************
! Initialize dqfact; dq1, dq2, and dq3 are weighted by dq(in1) and dq(in2)
           if (itheory .eq. 1) then 
            dqfact(1) = -dq1/(2.0d0*dq(in1))
            dqfact(2) =  dq1/(2.0d0*dq(in1))
            dqfact(3) = -dq2/(2.0d0*dq(in2)) 
            dqfact(4) =  dq2/(2.0d0*dq(in2)) 
            dqfact(5) = -dq3/(2.0d0*dq(indna))
            dqfact(6) =  dq3/(2.0d0*dq(indna)) 
            
            interaction = 2 
            do isorp = 1, ideriv_max 
             call Dtrescentros (interaction,isorp,ideriv_max, in1, in2, &
     &                          indna, x, y, cost, eps, depsA, depsB,   &
     &                          rhat, sighat, bcxcx_ca, f3xcXa_ca,      &
     &                          f3xcXb_ca, f3xcXc_ca, nspecies) 
     
! The forces f3xcXa_ca, f3xcXb_ca, and f3xcXc_ca are already force-like.  
             do inu = 1, num_orb(in2) 
              do imu = 1, num_orb(in1) 
               do ix = 1, 3
!$omp atomic
                f3xca_ca(ix,ialp) = f3xca_ca(ix,ialp) +                 &
     &                              2*rho(imu,inu,mneigh,iatom)*          &
     &                              f3xcXa_ca(ix,imu,inu)*dqfact(isorp)
!$omp atomic
                f3xcb_ca(ix,iatom) = f3xcb_ca(ix,iatom) +               &
     &                               2*rho(imu,inu,mneigh,iatom)*         &
     &                               f3xcXb_ca(ix,imu,inu)*dqfact(isorp)
!$omp atomic
                f3xcc_ca(ix,jatom) = f3xcc_ca(ix,jatom) +               &
     &                               2*rho(imu,inu,mneigh,iatom)*         &
     &                               f3xcXc_ca(ix,imu,inu)*dqfact(isorp)
               end do
              end do 
             end do 
            end do 
           end if 

! ****************************************************************************
! End loop over ialp and its common neighbors.
          end if
         end do
        end do

! Format Statements
! ===========================================================================
 
        return
        end
