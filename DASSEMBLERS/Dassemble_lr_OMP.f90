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
! 
! subroutine rewritten by P. Jelinek (openMP optimization) 
! email: jelinekp@fzu.cz
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_lr (natoms, nprocs, impi)
        use charges
        use configuration
        use constants_fireball
        use density
        use dimensions
        use forces
        use interactions
        use neighbor_map
        use omp_lib
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: impi
        integer, intent(in) :: natoms
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
        integer jmu
        integer katom
        integer mbeta
        integer my_proc
        integer natomsp

        real distance12
        real dq1
        real dq2
        real dq3
        real dterm
        real sterm
 
        real, dimension (3) :: dewaldlr_i
        real, dimension (3) :: dewaldlr_j
        real, dimension (3) :: dewaldlr_k
        real, dimension (3) :: dpterm
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: rhat12
        real, dimension (3) :: spterm
        real, dimension (natoms) :: sub_ewald
        real, dimension (3, natoms) :: sub_dewald
        real, allocatable, dimension (:,:,:) :: flrewl
        integer nth
        integer ith

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Procedure
! ===========================================================================
! Initialize the flrew array.
        flrew = 0.0d0

! get number of threads
        nth = omp_get_max_threads()
! allocate aux array
        allocate ( flrewl(3,natoms,nth) )
        flrewl = 0.0d0

! Determine which atoms are assigned to this processor.
        if (impi .eq. 1) then
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
 
! Initialize some arrays to zero.
        sub_ewald = 0.0d0
        sub_dewald = 0.0d0
 
! First a preliminary quantity. ewald(i,j) = sum(L) 1/(| bi-bj +/- L|
! We need to calculate SUMS of ewald sums in which the charge is included.
! These will also be used later in Pnanl3pP.f.  The details are on page 10 of
! the "Derivative of the long range matrix element". In the nutshell we are
! calculating:
!       Sum_(i,L) q(i)/|b(i)-b(alpha)+L|
 
!$omp parallel do private (jatom, in2, dq2, issh)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         do jatom = 1, natoms
          in2 = imass(jatom)
 
! Calculate the charge on jatom
          dq2 = 0.0d0
          do issh = 1, nssh(in2)
           dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
          end do
 
          sub_ewald(iatom) = sub_ewald(iatom) + dq2*ewald(iatom,jatom)
          sub_dewald(:,iatom) = sub_dewald(:,iatom) + dq2*dewald(:,iatom,jatom)
          !write (*,*) 'atoms :',iatom,jatom,' sub_ewald = ',sub_ewald(iatom),sub_dewald(:,iatom)
         end do ! do jatom
        end do ! do iatom
 
! sum sub_ewald, sub_dewald over procs
        if (impi .eq. 1)                                                 &
     &   call assemble_ordern_sub_dewald (natoms, sub_ewald, sub_dewald)

! Now the meat of the calculation.  Construct ewaldlr(mu,nu,i,m) ===>
! the matrix elements of the long-range parts of the Hamiltonian.
! We make matrix elements for the Long Range Ewald according to our theory:
! ewaldlr(mu,nu,iatom,ineigh) =
! {s(mu,nu,iatom,ineigh)/2}*SUM(j_basis)(Qin(jatom) - Qneutral(jatom))
!                                        *(ewald(iatom,jatom)
!                                          + ewald(ineigh,jatom))*eq2
! eq2 makes it into the units of eV.
! Loop over the atoms in the central cell.

!$omp parallel do private (r1, in1, dq1, issh, ineigh, jatom, mbeta, r2, in2) &
!$omp&  private (dq2, distance12, rhat12, inu, imu, sterm, dterm, dpterm)     &
!$omp&  private (spterm, dewaldlr_i, dewaldlr_j, jmu, katom, in3, dq3)        &
!$omp&  private (dewaldlr_k, ith) 
        do iatom = iatomstart, iatomstart - 1 + natomsp

! get thread id    
         ith = omp_get_thread_num () + 1
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)
 
! Calculate the charge on iatom
         dq1 = 0.0d0
         do issh = 1, nssh(in1)
          dq1 = dq1 + (Qin(issh,iatom) - Qneutral(issh,in1))
         end do
 
         !write (*,*) 'iatom = ',iatom, dq1
! Loop over the neighbors of the atom i.
         do ineigh = 1, neighn(iatom)
          jatom = neigh_j(ineigh,iatom)
          mbeta = neigh_b(ineigh,iatom)
 
! Position of jatom.
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)
 
! Find charge on jatom
          dq2 = 0.0d0
          do issh = 1, nssh(in2)
           dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
          end do
          !write (*,*) '   ineigh = ',ineigh,jatom, dq2
          distance12 = sqrt((r2(1) - r1(1))**2 + (r2(2) - r1(2))**2  &
     &                                         + (r2(3) - r1(3))**2)
          if (distance12 .gt. 1.0d-4) then
           rhat12(:) = (r2(:) - r1(:))/distance12
          end if
 
! "Charge" on each atom of the bondcharge. We split the charge S and
! dipole p to be S/2-p/d on atom 1 and S/2+p/d on atom 2. If atom 1 is
! equal to atom 2, then drop the p term.
          if (distance12 .gt. 1.0d-4) then
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             sterm = s_mat(imu,inu,ineigh,iatom)/2.0d0
!            if (distance12 .gt. 1.0d-4) then
             dterm = dip(imu,inu,ineigh,iatom)/distance12
             dpterm(:) = dipp(:,imu,inu,ineigh,iatom)/distance12             &
     &                   + dip(imu,inu,ineigh,iatom)*rhat12(:)/distance12**2
             spterm(:) = sp_mat(:,imu,inu,ineigh,iatom)/2.0d0

             dewaldlr_i(:) = (sterm - dterm)*sub_dewald(:,iatom)             &
     &        + (spterm(:) - dpterm(:))*sub_ewald(iatom)                     &
     &        + (sterm + dterm)*dewald(:,iatom,jatom)*dq1                    &
     &        + (spterm(:) + dpterm(:))*sub_ewald(jatom )
             dewaldlr_j(:) = (sterm - dterm)*dewald(:,jatom,iatom)*dq2       &
     &        - (spterm(:) - dpterm(:))*sub_ewald(iatom)                     &
     &        + (sterm + dterm)*sub_dewald(:,jatom)                          &
     &        - (spterm(:) + dpterm(:))*sub_ewald(jatom)
!            else
!             dterm = 0.0d0
!             dpterm(:) = 0.0d0
!             spterm(:) = 0.0d0

!             dewaldlr_i(:) = sterm*sub_dewald(:,iatom)
!             dewaldlr_j(:) = sterm*sub_dewald(:,jatom)
!            end if
 
! Combine all of the contributions together before multiplying by rho.
! The - sign makes flrew force-like.
             do jmu = 1, 3
              flrewl(jmu,iatom,ith) = flrewl(jmu,iatom,ith)                  &
     &         - rho(imu,inu,ineigh,iatom)*dewaldlr_i(jmu)*eq2
              flrewl(jmu,jatom,ith) = flrewl(jmu,jatom,ith)                  &
     &         - rho(imu,inu,ineigh,iatom)*dewaldlr_j(jmu)*eq2
             end do

! ****************************************************************************
! Now find the force on the other atoms katom, which are implicit in the
! ewald sums. Here katom is neither iatom nor jatom.
! ****************************************************************************
             do katom = 1, natoms
              if (katom .ne. iatom .and. katom .ne. jatom) then
               in3 = imass(katom)
 
! Calculate the charge on katom
               dq3 = 0.0d0
               do issh = 1, nssh(in3)
                dq3 = dq3 + (Qin(issh,katom) - Qneutral(issh,in3))
               end do
 
               dewaldlr_k(:) = ((sterm-dterm)*dewald(:,katom,iatom)          &
     &            + (sterm+dterm)*dewald(:,katom,jatom))*dq3

               do jmu = 1, 3
                flrewl(jmu,katom,ith) = flrewl(jmu,katom,ith)                &
     &           - rho(imu,inu,ineigh,iatom)*dewaldlr_k(jmu)*eq2
               end do
              end if
             end do ! do katom
 
! ****************************************************************************
! End loop over imu, inu
            end do ! do imu
           end do ! do inu
          else 
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             sterm = s_mat(imu,inu,ineigh,iatom)/2.0d0
             dterm = 0.0d0
             dpterm(:) = 0.0d0
             spterm(:) = 0.0d0

             dewaldlr_i(:) = sterm*sub_dewald(:,iatom)
             dewaldlr_j(:) = sterm*sub_dewald(:,jatom)

! Combine all of the contributions together before multiplying by rho.
! The - sign makes flrew force-like.
             do jmu = 1, 3
              flrewl(jmu,iatom,ith) = flrewl(jmu,iatom,ith)                  &
     &         - rho(imu,inu,ineigh,iatom)*dewaldlr_i(jmu)*eq2
              flrewl(jmu,jatom,ith) = flrewl(jmu,jatom,ith)                  &
     &         - rho(imu,inu,ineigh,iatom)*dewaldlr_j(jmu)*eq2
             end do

! ****************************************************************************
! Now find the force on the other atoms katom, which are implicit in the
! ewald sums. Here katom is neither iatom nor jatom.
! ****************************************************************************
             do katom = 1, natoms
              if (katom .ne. iatom .and. katom .ne. jatom) then
               in3 = imass(katom)

! Calculate the charge on katom
               dq3 = 0.0d0
               do issh = 1, nssh(in3)
                dq3 = dq3 + (Qin(issh,katom) - Qneutral(issh,in3))
               end do

               dewaldlr_k(:) = ((sterm-dterm)*dewald(:,katom,iatom)          &
     &            + (sterm+dterm)*dewald(:,katom,jatom))*dq3

               do jmu = 1, 3
                flrewl(jmu,katom,ith) = flrewl(jmu,katom,ith)                &
     &           - rho(imu,inu,ineigh,iatom)*dewaldlr_k(jmu)*eq2
               end do
              end if
             end do ! do katom
! ****************************************************************************
! End loop over imu, inu
            end do ! do imu
           end do ! do inu
          end if ! if(distance12)
! End loop over iatom and it neighbors
         end do ! do ineigh
        end do ! do iatom

! sum results of all threads
        do ith = 1,nth
         do iatom = 1,natoms
          flrew(:,iatom) = flrew(:,iatom) + flrewl(:,iatom,ith)
         end do ! do iatom
        end do ! do ith

! deallocate local arrays
        deallocate ( flrewl )


! Format Statements
! ===========================================================================
        return
        end
