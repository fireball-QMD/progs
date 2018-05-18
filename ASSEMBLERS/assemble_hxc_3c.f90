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

 
! assemble_hxc_3c.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the true three-center interactions
! for the Horsfield exchange-correlation interactions.  These are true 
! three-center in that iatom .ne. jatom .ne. katom.  The matrix elements look 
! like <Psi(1)|V(3)|Psi(2)>.
!
! The bulk of the work is done in trescentros.f. This program assembles the
! results.
!
! A third party term is when we consider the NA (or etc.) to be at the origin
! and take the matrix element between a pair of neighbors, neither of which is
! the NA (or etc.), and the pair is truly a pair, and not an atom.
! the results are stored in: f3na(ix,i), f3xc(ix,i), etc.
!
! The generic form of the matrix element is
! v = <mu(r-(r1+ratm)) ! v(r-ratm) ! nu(r-(r2+ratm))>.
!
! See notes: "general theory of third party terms" for the definition on p. 2
! of these three terms.
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
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_hxc_3c (nprocs, Kscf, iordern, &
     &                              itheory, igauss)
        use charges
        use configuration
        use constants_fireball
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
        integer, intent (in) :: Kscf
        integer, intent (in) :: nprocs
 
 
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
        integer jatom
        integer jneigh
        integer jbeta
        integer mneigh
        integer my_proc
        integer natomsp

        integer matom
        integer mbeta
        integer in3 
 
        real cost
        real epsx
        real epsxc
        real dq1
        real dq2
        real dq3
        real dpotxc_2c
        real dpotxc_3c
        real dpotxc_total
        real drvexc
        real potx
        real potxc_2c
        real potxc_3c
        real potxc_total
        real rho_2c
        real rho_3c
        real rho_total
        real x
        real y

        real, dimension (numorb_max, numorb_max) :: bcxcx
        real, dimension (3, numorb_max, numorb_max) :: bcxcpx

        real, dimension (6) :: dqfact          ! there are six different isorps
        real, dimension (3, 3) :: eps
        real, dimension (3, 3, 3) :: deps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: rhat
        real, dimension (3) :: rna
        real, dimension (3) :: rnabc
        real, dimension (3) :: sighat

        real, dimension (numorb_max, numorb_max) :: smG
        real, dimension (3, numorb_max, numorb_max) :: spmG

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ========================================================================
! The array vxc will not be initialized to zero here.
! Presumably, the two-center interactions have already been calculated.
! As such, at this point and time these arrays should not be zero.
        if (igauss .eq. 1) then
         density_2c = 0.0d0
         density_3c = 0.0d0
         bar_density_2c = 0.0d0
         bar_density_3c = 0.0d0
         vxc_3c = 0.0d0
          smG = 0.d0
         spmG = 0.d0
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
!$omp&            private (bcxcx, dqfact, eps, r1, r2, r21, rhat, rna, rnabc)&
!$omp&            private (sighat)
        do ialp = iatomstart, iatomstart - 1 + natomsp
         rna(:) = ratom(:,ialp)
         indna = imass(ialp)

! Find the charge of the third-center
         dq3 = 0.0d0
         do issh = 1, nssh(indna)
          dq3 = dq3 + (Qin(issh,ialp) - Qneutral(issh,indna))
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
 

! ****************************************************************************
!
! PERFORM ACTUAL CALCULATIONS
! CALL TRESCENTROS FOR EXCHANGE-CORRELATION PIECE - NA CASE
! ****************************************************************************
           if (Kscf .eq. 1) then
            isorp = 0
            interaction = 2
            if (igauss .eq. 0) then
             call trescentros (interaction, isorp, ideriv_max, in1, in2,     &
     &                         indna, x, y, cost, eps, bcxcx, nspecies)

! Add this piece for iatom, jatom, and katom into the total - bcxc
!$omp critical (vxc3)
             do inu = 1, num_orb(in2)
              do imu = 1, num_orb(in1)
               vxc(imu,inu,mneigh,iatom) =                              &
     &          vxc(imu,inu,mneigh,iatom) + bcxcx(imu,inu)
             vxc(inu,imu,jneigh,jatom) = vxc(imu,inu,mneigh,iatom)

              end do
             end do
!$omp end critical (vxc3)
            else 
             call trescentrosGHXC_VXC (in1, in2, indna, x, y, cost,     &
     &                                 eps, bcxcx, rcutoff)
!$omp critical (den3)
             do inu = 1, num_orb(in2)
              do imu = 1, num_orb(in1)
               density_3c(imu,inu,mneigh,iatom) =                       &
     &          density_3c(imu,inu,mneigh,iatom) + bcxcx(imu,inu)
         density_3c(inu,imu,jneigh,jatom) = density_3c(imu,inu,mneigh,iatom) 
              end do
             end do
!$omp end critical (den3)
            end if
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
             call trescentros (interaction, isorp, ideriv_max, in1, in2,     &
     &                         indna, x, y, cost, eps, bcxcx, nspecies) 
     
! Now add bcxc to vxc_ca.  
             do inu = 1, num_orb(in2) 
              do imu = 1, num_orb(in1) 
!$omp atomic
               vxc_ca(imu,inu,mneigh,iatom) =                                &
     &          vxc_ca(imu,inu,mneigh,iatom) + bcxcx(imu,inu)*dqfact(isorp) 
              end do 
             end do 
            end do
           end if

 
! ****************************************************************************
! End loop over ialp and its common neighbors.
          end if
         end do
        end do

! Approximation for three-center interactions using gaussian fits to 
! wavefunctions.
! ****************************************************************************
! Calculate three-center exchange-correlation
! vxc_3c = xc(nbar_1 + nbar_2 + nbar_3c) - xc(nbar_1 + nbar_2)
! Here, nbar_i = density/s

! XC --- LDA - Ceperley - Alder exchange-correlation potential and energy
! as parameterized by Perdew and Zunger, Phys. Rev. B23, 5048 (1981)
! ****************************************************************************
        if (igauss .eq. 1) then
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

           if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in unocentros.f

           else

! We need to find the two-center density and subtract it off the three-center density.

! For the atom case: these loops are both over in1 indeed!  This is because
! the two wavefunctions are located on the same site, but the potential is
! located at a different site.
             cost = 0.0d0
             in3 = in1
             call dosgaussians (in1, in2, in3, y, cost, eps, deps,      &
     &                          bcxcx, bcxcpx, rcutoff)

             do inu = 1, num_orb(in3)
              do imu = 1, num_orb(in1)
               density_2c(imu,inu,matom,iatom) = bcxcx(imu,inu)
              end do
             end do

! For the ontop left case, the density is in the first atom (iatom):
!             cost = -1.0d0
!             in3 = in2
!             call dosgaussians (in1, in1, in3, y, cost, eps, deps,      &
!     &                          bcxcx, bcxcpx, rcutoff)
!
!             do inu = 1, num_orb(in3)
!              do imu = 1, num_orb(in1)
!              density_2c(imu,inu,ineigh,iatom) =                       &
!     &         density_2c(imu,inu,ineigh,iatom) + bcxcx(imu,inu)
!              end do
!             end do

! For the ontop right case, the density is in the second atom (iatom):
             cost = 1.0d0
             in3 = in2
             call dosgaussians (in1, in2, in3, y, cost, eps, deps,      &
     &                          bcxcx, bcxcpx, rcutoff)
             do inu = 1, num_orb(in3)
              do imu = 1, num_orb(in1)
               density_2c(imu,inu,ineigh,iatom) =                       &
     &         density_2c(imu,inu,ineigh,iatom) + bcxcx(imu,inu)
              end do
             end do

! overlap matrix for gaussian
           call doscentrosG_overlap (in1, in2, y, eps, deps, smG, spmG, &
     &                               rcutoff)


           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
 
! Calculate nbar
             if (abs(smG(imu,inu)) .lt. xc_overtolG) then
              rho_3c = 0.0d0
              rho_2c = 0.0d0

! bar_density_2c means <phi(1)|(N1_2c + N2_2c)|phi(2)>/<phi(1)|phi(2)>
! bar_density_3c means <phi(1)| N_3c |phi(2)>/<phi(1)|phi(2)>
             else
             rho_3c = density_3c(imu,inu,ineigh,iatom)/smG(imu,inu)
             rho_2c = density_2c(imu,inu,ineigh,iatom)/smG(imu,inu)
             end if

! Get the average electron density.
! Setup matries bardens4_3c and bardens4n12 for force calculation
             if (rho_3c .lt. 0.0d0) rho_3c =0.0d0
             if (rho_2c .lt. 0.0d0) rho_2c = 0.0d0
             bar_density_3c(imu,inu,ineigh,iatom) = rho_3c
             bar_density_2c(imu,inu,ineigh,iatom) = rho_2c

! CALL subroutine xc-capal for total density and the two-center density.
! LDA - Ceperley - Alder exchange-correlation potential and energy
! as parameterized by Perdew and Zunger, Phys. Rev. B23, 5048 (1981)
             rho_total = rho_3c + rho_2c
             call ceperley_alder (rho_total, epsxc, potxc_total, dpotxc_total)
             call ceperley_alder (rho_2c, epsxc, potxc_2c, dpotxc_2c)
 
! three-center correction
             potxc_3c = potxc_total - potxc_2c

             vxc_3c(imu,inu,ineigh,iatom) = potxc_3c*smG(imu,inu)
!             vxc_3c(imu,inu,ineigh,iatom) = 0.d0
! The three-center exchange-correlation correction matrix is:
! vxc_3c = vxc_3c * S 
! Finally, inculding the two-center exchange-correlation contribution
             vxc(imu,inu,ineigh,iatom) = vxc(imu,inu,ineigh,iatom) +    &
     &                                vxc_3c(imu,inu,ineigh,iatom)

! The derivative of exchange-correlation energy matrix respect to density are
! in atomic units. We should convert it to hartree (eV) units.
! Finally, the units of exchange-correlation potential matrix dpotxc* is eV*A**2
             dpotxc_3c = dpotxc_total - dpotxc_2c
 
             nuxc_3c(imu,inu,ineigh,iatom) = dpotxc_3c
             nuxc_total(imu,inu,ineigh,iatom) = dpotxc_total
            end do
           end do

          end if ! end if of iatom .eq. jatom .and. mbeta .eq. 0

          end do ! end do of ineigh
         end do  ! end do of iatom
        end if   ! end if of igauss

! Format Statements
! ===========================================================================
118      format(4(2x,i4), 2(2x,f16.8))
 
        return
        end
