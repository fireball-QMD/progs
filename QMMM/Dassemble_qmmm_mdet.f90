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

 
! Dassemble_lr.f90
! Program Description
! ===========================================================================
!       This routine calculated the long long-range ewald forces -
! flrew(3, natoms).
!
! ===========================================================================
! Code originally written by:
! Alex A. Demkov
! Predictive Engineering Laboratory
! Motorola, Inc.  M360
! Mesa, AZ 85202 USA
! (602) 655-2454, r37591@email.sps.mot.com
 
! and Otto F. Sankey
! Campus Box 1504
! Department of Physics
! Arizona State University
! Tempe, AZ 85287-1504
! (602) 965-4334 (office)      email: otto.sankey@asu.edu
 
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
        subroutine Dassemble_qmmm_mdet (nprocs, iordern)
        use charges
        use configuration
        use constants_fireball
        use density
        use dimensions
        use forces
        use interactions
        use neighbor_map
        use nonadiabatic
        use qmmm_module, only : qmmm_struct, qmmm_nml
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: iordern
        integer, intent(in) :: nprocs
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iatomstart
        integer ierror
        integer imu
        integer in1
        integer in2
        integer in3
        integer ineigh
        integer inu
        integer issh
        integer jatom
        integer gatom
        integer jmu
        integer katom
        integer mbeta
        integer my_proc
        integer natomsp

        real distance12
        real dij
        real dq1
        real dq2
        real dq3
        real dterm
        real sterm
!        real out_charge
 
        real, dimension (3) :: dewaldlr_i_qmmm
        real, dimension (3) :: dewaldlr_j_qmmm
        real, dimension (3) :: dewaldlr_k
        real, dimension (3) :: dpterm
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: rhat12
        real, dimension (3) :: spterm
        real, dimension (3) :: vij
        real, dimension (qmmm_struct%qm_mm_pairs) :: mm_charges
!        real, dimension (3,natoms,qmmm_struct%qm_mm_pairs) :: dqmmm
        real, dimension (:,:,:), allocatable :: dqmmm 
        real, dimension (natoms) :: sub_ewaldqmmm
        real, dimension (3, natoms) :: sub_dewaldqmmm

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
       

! Procedure
! ===========================================================================
! Initialize the flrew array.
        flrew_qmmm = 0.0d0

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

        allocate(dqmmm(3,natoms,qmmm_struct%qm_mm_pairs))

        sub_ewaldqmmm = 0.0d0
        sub_dewaldqmmm = 0.0d0
        dqmmm = 0.0d0        
        qmmm_struct%dxyzcl = 0.0d0

!        out_charge = 0
!        do katom = 1, qmmm_struct%qm_mm_pairs
!          out_charge = out_charge + qmmm_struct%qm_xcrd(4,katom)
!        end do
!        write(*,*) 'Total mm charge = ', out_charge
!        out_charge = out_charge / qmmm_struct%qm_mm_pairs
!        write(*,*) 'Delta charge = ', out_charge
!        do katom = 1, qmmm_struct%qm_mm_pairs
!          mm_charges(katom) = qmmm_struct%qm_xcrd(4,katom) - out_charge
!        end do

!        do gatom = 1, qmmm_struct%qm_mm_pairs
!          if ( sqrt(qmmm_struct%qm_xcrd_dist(gatom)) > 0.8*qmmm_nml%qmcut ) then
!            qmmm_struct%qm_xcrd(4,gatom) = qmmm_struct%qm_xcrd(4,gatom)*(((qmmm_nml%qmcut - sqrt(qmmm_struct%qm_xcrd_dist(gatom)))/(0.2*qmmm_nml%qmcut)))
!          end if
!        end do



        do iatom = iatomstart, iatomstart - 1 + natomsp
          do katom = 1, qmmm_struct%qm_mm_pairs

            dij = sqrt ( (qmmm_struct%qm_xcrd(1,katom)-ratom(1,iatom))*(qmmm_struct%qm_xcrd(1,katom)-ratom(1,iatom)) + &
                         (qmmm_struct%qm_xcrd(2,katom)-ratom(2,iatom))*(qmmm_struct%qm_xcrd(2,katom)-ratom(2,iatom)) + & 
                         (qmmm_struct%qm_xcrd(3,katom)-ratom(3,iatom))*(qmmm_struct%qm_xcrd(3,katom)-ratom(3,iatom)) )

              sub_ewaldqmmm(iatom) = sub_ewaldqmmm(iatom) - (qmmm_struct%qm_xcrd(4,katom)/dij)
            do jmu = 1, 3
              vij(jmu) = qmmm_struct%qm_xcrd(jmu,katom)-ratom(jmu,iatom)
              sub_dewaldqmmm(jmu,iatom) = sub_dewaldqmmm(jmu,iatom) - qmmm_struct%qm_xcrd(4,katom)*(vij(jmu) / dij**3)
              dqmmm(jmu,iatom,katom) = -(vij(jmu) / dij**3)
            end do
!            write(*,*) 'distance',iatom,katom
!            write(*,*)  vij
          end do
        end do

! Now the meat of the calculation.  Construct ewaldlr(mu,nu,i,m) ===>
! the matrix elements of the long-range parts of the Hamiltonian.
! We make matrix elements for the Long Range Ewald according to our theory:
! ewaldlr(mu,nu,iatom,ineigh) =
! {s(mu,nu,iatom,ineigh)/2}*SUM(j_basis)(Qin(jatom) - Qneutral(jatom))
!                                        *(ewald(iatom,jatom)
!                                          + ewald(ineigh,jatom))*eq2
! eq2 makes it into the units of eV.
! Loop over the atoms in the central cell.
!$omp parallel do private (in1, in2, in3, jatom, mbeta, distance12, dq1, dq2)&
!$omp&            private (dq3, dterm, sterm, dewaldlr_i_qmmm, dewaldlr_j_qmmm)        &
!$omp&            private (dewaldlr_k, dpterm, r1, r2, rhat12, spterm)

        do iatom = iatomstart, iatomstart - 1 + natomsp
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)
! Loop over the neighbors of the atom i.
         do ineigh = 1, neighn(iatom)
          jatom = neigh_j(ineigh,iatom)
          mbeta = neigh_b(ineigh,iatom)
 
! Position of jatom.
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)
 
          distance12 = sqrt((r2(1) - r1(1))**2 + (r2(2) - r1(2))**2  &
     &                                         + (r2(3) - r1(3))**2)
          if (distance12 .gt. 1.0d-4) then
           rhat12(:) = (r2(:) - r1(:))/distance12
          end if
 
! "Charge" on each atom of the bondcharge. We split the charge S and
! dipole p to be S/2-p/d on atom 1 and S/2+p/d on atom 2. If atom 1 is
! equal to atom 2, then drop the p term.
          do inu = 1, num_orb(in2)
           do imu = 1, num_orb(in1)
            sterm = s_mat(imu,inu,ineigh,iatom)/2.0d0
            if (distance12 .gt. 1.0d-4) then
             dterm = dip(imu,inu,ineigh,iatom)/distance12
             dpterm(:) = dipp(:,imu,inu,ineigh,iatom)/distance12             &
     &                   + dip(imu,inu,ineigh,iatom)*rhat12(:)/distance12**2
             spterm(:) = sp_mat(:,imu,inu,ineigh,iatom)/2.0d0

             dewaldlr_i_qmmm(:) = (sterm - dterm)*sub_dewaldqmmm(:,iatom)             &
     &        + (spterm(:) - dpterm(:))*sub_ewaldqmmm(iatom)                     &
     &        + (spterm(:) + dpterm(:))*sub_ewaldqmmm(jatom)
             dewaldlr_j_qmmm(:) = -(spterm(:) - dpterm(:))*sub_ewaldqmmm(iatom)                     &
     &        + (sterm + dterm)*sub_dewaldqmmm(:,jatom)                          &
     &        - (spterm(:) + dpterm(:))*sub_ewaldqmmm(jatom)
            else
             dterm = 0.0d0
             dpterm(:) = 0.0d0
             spterm(:) = 0.0d0

             dewaldlr_i_qmmm(:) = sterm*sub_dewaldqmmm(:,iatom)
             dewaldlr_j_qmmm(:) = sterm*sub_dewaldqmmm(:,jatom)
            end if
            
!            write(*,*) 'iatom',iatom
!            write(*,*) 'dewaldlr_i_qmmm',dewaldlr_i_qmmm
!            write(*,*) 'jatom',jatom
!            write(*,*) 'dewaldlr_j_qmmm',dewaldlr_j_qmmm        

! Combine all of the contributions together before multiplying by rho.
! The - sign makes flrew force-like.
            do jmu = 1, 3
!$omp atomic
             flrew_qmmm(jmu,iatom) = flrew_qmmm(jmu,iatom)                             &
     &        - rho(imu,inu,ineigh,iatom)*dewaldlr_i_qmmm(jmu)*eq2
!$omp atomic
             flrew_qmmm(jmu,jatom) = flrew_qmmm(jmu,jatom)                             &
     &        - rho(imu,inu,ineigh,iatom)*dewaldlr_j_qmmm(jmu)*eq2
            end do
 
! add to gh_lrew_qmmm           
            do jmu = 1, 3
            gh_lrew_qmmm(jmu,iatom,imu,inu,ineigh,iatom) =                     &
     &      gh_lrew_qmmm(jmu,iatom,imu,inu,ineigh,iatom) +                     &
     &        dewaldlr_i_qmmm(jmu)*eq2
            gh_lrew_qmmm(jmu,jatom,imu,inu,ineigh,iatom) =                     &
     &      gh_lrew_qmmm(jmu,jatom,imu,inu,ineigh,iatom) +                     &
     &        dewaldlr_j_qmmm(jmu)*eq2
            end do


! ****************************************************************************
! Now find the force on the rho of katoms (mm atoms)
! ****************************************************************************

            do katom = 1, qmmm_struct%qm_mm_pairs
             do jmu = 1, 3
              qmmm_struct%dxyzcl(jmu,katom) = qmmm_struct%dxyzcl(jmu,katom) - &
              &    ((sterm-dterm)*dqmmm(jmu,iatom,katom)+(sterm+dterm)*dqmmm(jmu,jatom,katom))*  &
              &    rho(imu,inu,ineigh,iatom)*qmmm_struct%qm_xcrd(4,katom)*eq2*23.061d0
!              gh_mm_qmmm(jmu,katom,imu,inu,ineigh,iatom) = gh_mm_qmmm(jmu,katom,imu,inu,ineigh,iatom) - &
!              &    ((sterm-dterm)*dqmmm(jmu,iatom,katom)+(sterm+dterm)*dqmmm(jmu,jatom,katom))* &
!              &    qmmm_struct%qm_xcrd(4,katom)*eq2
             end do
            end do


! End loop over imu, inu
           end do
          end do
 
! End loop over iatom and it neighbors
         end do
        end do

! ****************************************************************************
! Now find the force on the core of katoms (mm atoms)
! ****************************************************************************

        do iatom = iatomstart, iatomstart - 1 + natomsp
          in3 = imass(iatom)
          dq3 = 0.0d0
          do issh = 1, nssh(in3)
               dq3 = dq3  + Qneutral(issh,in3)
          end do
          do katom = 1, qmmm_struct%qm_mm_pairs
            do jmu = 1, 3
              qmmm_struct%dxyzcl(jmu,katom) = qmmm_struct%dxyzcl(jmu,katom) + &
              &        dq3*qmmm_struct%qm_xcrd(4,katom)*dqmmm(jmu,iatom,katom)*eq2*23.061d0
              flrew_qmmm(jmu,iatom) = flrew_qmmm(jmu,iatom) + &
              &        dq3*qmmm_struct%qm_xcrd(4,katom)*dqmmm(jmu,iatom,katom)*eq2
 
!              gh_lrew_qmmm(jmu,iatom,imu,inu,ineigh,iatom) = gh_lrew_qmmm(jmu,iatom,imu,inu,ineigh,iatom) - &
!              &        qmmm_struct%qm_xcrd(4,katom)*dqmmm(jmu,iatom,katom)*eq2
!              gh_mm_qmmm(jmu,katom,imu,inu,ineigh,iatom) = gh_mm_qmmm(jmu,katom,imu,inu,ineigh,iatom) + &
!              &        qmmm_struct%qm_xcrd(4,katom)*dqmmm(jmu,iatom,katom)*eq2
            end do
          end do          
        end do

        deallocate(dqmmm)

! Format Statements
! ===========================================================================
        return
        end
