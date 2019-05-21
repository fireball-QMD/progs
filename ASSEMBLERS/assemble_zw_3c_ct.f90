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

 
! assemble_zw_3c_ct.f90
! Program Description
! ===========================================================================
!
! ===========================================================================
 
! Program Declaration
! ===========================================================================
        subroutine assemble_zw_3c_ct (nprocs, iordern, igauss)
        use charges
        use configuration
        use constants_fireball
        use dimensions
        use interactions
        use neighbor_map
        use gaussG

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
        integer kneigh
        integer ineigh1
        integer ineigh2
        integer jbeta
        integer jcount
        integer jcount_sav
        integer jssh
        integer issh1
        integer issh2
        integer matom
        integer mbeta
        integer mneigh
        integer my_proc
        integer natomsp
        integer ix
        integer j
 
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
        real dot_product_dipc_x
        real :: A,B
 
        real, dimension (numorb_max, numorb_max) :: bcca
        real, dimension (numorb_max, numorb_max) :: bccax
        real, dimension (numorb_max, numorb_max) :: emnpl
        real, dimension (3, 3, 3) :: deps
        real, dimension (3, 3) :: eps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: rhat
        real, dimension (3) :: rna
        real, dimension (3) :: rnabc
        real, dimension (3) :: sighat
! JIMM-JOM
!       real, dimension (numorb_max, numorb_max) :: stn1
!       real, dimension (numorb_max, numorb_max) :: stn2
        real stn1
        real stn2

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
 
     !     rcutoff_ialp = 0.0d0
     !    do imu = 1, nssh(indna)
     !     if (rcutoff(indna,imu) .gt. rcutoff_ialp)                     &
     !      rcutoff_ialp = rcutoff(indna,imu)
     !    end do
 
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
! ENRIQUE-JOM: calculate rcutoff_i for smoother
     !     rcutoff_i = 0
     !     do imu = 1, nssh(in1)
     !      if (rcutoff(in1,imu) .gt. rcutoff_i) rcutoff_i = rcutoff(in1,imu)
     !     end do
! End ENRIQUE-JOM 
 
           jatom = neigh_comj(2,ineigh,ialp)
           jbeta = neigh_comb(2,ineigh,ialp)
           ineigh2 = neigh_com_ng(2,ineigh,ialp)
           r2(:) = ratom(:,jatom) + xl(:,jbeta)
           in2 = imass(jatom)
           jneigh = neigh_back(iatom,mneigh)
 
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
 
! ****************************************************************************
!
! PERFORM ACTUAL CALCULATIONS
! CALL TRESCENTROS FOR NEUTRAL ATOM PIECE
! ****************************************************************************
! Initialize bcca
           bcca = 0.0d0
           bccax=0.0d0
           do isorp = 1, nssh(indna)

     !      interaction =   !???
     !       call trescentros (interaction, isorp, isorpmax, in1, in2,   &
     !&                        indna, x, y, cost, eps, bccax, nspecies)

               !We can use the two center integrals here computed in
               !zw_2c_ct !!   gcab=A*gcaa+B*gcbb, with A,B the
               !Mulliken-Dipole projection coefficients
                  
               !A*g2nu(imu,inu,ineigh1,ialp)+B*g2nu(imu,inu,ineigh2,ialp)
               !ineigh1 is the neighbor index of iatom with respect to
               !ialp and ineigh2 is the neighbor index of jatom with
               !respect to ialp.
               !Other possibility: Store g2nu in the form
               !g2nu(imu,inu,jatom,iatom), NOT
               !g2nu(imu,inu,ineigh,iatom) ..!?


! Find the charge associated with this shell
            dxn = (Qin(isorp,ialp) - Qneutral(isorp,indna))

! Add this piece for iatom, jatom, and ialp into the total - bcca
!            do inu = 1, num_orb(in2)
!             do imu = 1, num_orb(in1)
!              bcca(imu,inu) = bcca(imu,inu) + bccax(imu,inu)*dxn
!             end do
!            end do

 
        !  end do ! end do of isorp loop.  I've moved this loop to the
        !  end

! Now add bcca to vca, including the smoothing term emnpl.
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
!$omp atomic
              issh1=orb2shell(imu,in1)
              issh2=orb2shell(inu,in2)
              A=0.5*s_mat(imu,inu,mneigh,iatom)-dip(imu,inu,mneigh,iatom)/y
              B=0.5*s_mat(imu,inu,mneigh,iatom)+dip(imu,inu,mneigh,iatom)/y  
           !Here is the split of the INTEGRAL gcab...--->     
                 bccax(imu,inu) = bccax(imu,inu)+ &
      & (A*g2nu(isorp,issh1,ineigh1,ialp)+B*g2nu(isorp,issh2,ineigh2,ialp))*dxn
             end do !end do imu
           end do !end do inu
          end do ! end do of isorp loop
         
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             vxc_ca(imu,inu,mneigh,iatom) = vxc_ca(imu,inu,mneigh,iatom) +    &
                        &              bccax(imu,inu)
            !Symmetrize Hamiltonian (April 2018): jneigh is the
            !back_neigh:
              vxc_ca(inu,imu,jneigh,jatom) = vxc_ca(imu,inu,mneigh,iatom)
            end do  !end do inu
           end do  !end do imu
 
! ****************************************************************************
! End loop over ialp and its common neighbors.

          end if
         end do
        end do
! Format Statements
! ===========================================================================
601     format (9(f8.4))
602     format (9(f8.4))  

        return 
        end
