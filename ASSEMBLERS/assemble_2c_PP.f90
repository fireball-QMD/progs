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

 
! assemble_2c_PP.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the two-center and degenerate
! two-center interactions.
!
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
!
! Program Declaration
! ===========================================================================
        subroutine assemble_2c_PP (nprocs, iforce, iordern)
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
        integer, intent (in) :: iforce
        integer, intent (in) :: iordern
        integer, intent (in) :: nprocs
 
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
        integer ineigh
        integer inu
        integer isorp
        integer jatom
        integer jneigh
        integer kneigh
        integer matom
        integer mbeta
        integer mneigh_self
        integer my_proc
        integer natomsp
        integer ncc
 
        real, dimension (numorb_max) :: cl
        real, dimension (numorb_max, numorb_max) :: PPx

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Initialize interactions
        vnl = 0.0d0

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
 
! ****************************************************************************
!
! ASSEMBLE VNL ATM CASE  <phi_i|Psi_j><Psi_j|phi_i>
! ****************************************************************************
!$omp parallel do private (in1, matom, ineigh, jatom, mbeta, in2, cl)        &
!$omp&            private (imu, inu, ncc, PPx ) 
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neighPP_self(iatom)
         in1 = imass(iatom)
 
! Loop over the neighbors of each iatom.
         do ineigh = 1, nPPn(iatom)        !  <==== loop over i's neighbors
          mbeta = nPP_b(ineigh,iatom)
          jatom = nPP_j(ineigh,iatom)
          in2 = imass(jatom)
 

! Get the coeffiecients
          call cl_value (in2, cl)
 
! Now we combine and sum:
! in1 twice because it is an atom case.
          do inu = 1, num_orb(in1)
           do imu = 1, num_orb(in1)
            PPx(imu,inu) = 0.0d0
            do ncc = 1, num_orbPP(in2)
             PPx(imu,inu) = PPx(imu,inu)                                     &
     &        + cl(ncc)*sVNL(imu,ncc,ineigh,iatom)*sVNL(inu,ncc,ineigh,iatom)
            end do
           end do
          end do
 
! Final assembly of vnl - the energy piece.
          do inu = 1, num_orb(in1)
           do imu = 1, num_orb(in1)
!$omp atomic
            vnl(imu,inu,matom,iatom) = vnl(imu,inu,matom,iatom)              &
     &                                    + PPx(imu,inu)
           end do
          end do

         enddo ! do ineigh
        enddo ! do iatom

! ****************************************************************************
!
! ASSEMBLE VNL ONTOP LEFT CASE   <phi_i|Psi_i><Psi_i|phi_j>
! ****************************************************************************
! Loop over iatom
!$omp parallel do private (in1, matom, ineigh, jatom, mbeta, in2, cl)        &
!$omp&            private (mneigh_self, imu, inu, ncc, PPx, kneigh)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neighPP_self(iatom)
         in1 = imass(iatom)
 
! Loop over the neighbors of each iatom.
         do ineigh = 1, nPPxn(iatom)        !  <==== loop over i's neighbors
          mbeta = nPPx_b(ineigh,iatom)
          jatom = nPPx_j(ineigh,iatom)
          in2 = imass(jatom)



! if r1 .ne. r2, then this is a case where the potential is located
! at one of the sites of a wavefunction (ontop case).
          if (iatom .eq. jatom .and. mbeta .eq. 0) then
 
! Do nothing here - special case. Interaction already calculated in atm case.
! sanity check
           if (nPPx_self(iatom) .ne. ineigh) then
            write (*,*) ' Something really wrong in assemble_2c_PP.f90 '
            write (*,*) ' iatom, jatom, mbeta = ', iatom, jatom, mbeta
            write (*,*) ' neigh_self(iatom), ineigh = ',                     &
     &                    nPPx_self(iatom), ineigh  
! FIXME stop is not allowed in an OpenMP parallelized loop
!           stop
           end if ! if(neighPP_self)
          else ! if(iatom .eq. jatom)

! Case 1. PP is iatom.  <i | VNL(i) |j>.
           call cl_value (in1, cl)
           jneigh = nPPx_point(ineigh,iatom)
! <phi_i|Psi_i>  ->  nPP(mneigh_self,iatom)
           mneigh_self = nPP_self(iatom)

! Now we combine and sum:
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             PPx(imu,inu) = 0.0d0
             do ncc = 1, num_orbPP(in1)
              PPx(imu,inu) = PPx(imu,inu)                                    &
     &         + cl(ncc)*sVNL(imu,ncc,mneigh_self,iatom)                     &
     &                  *sVNL(inu,ncc,jneigh,jatom)
             end do ! do ncc
            end do ! do imu
           end do ! do inu

! Mapping to the global matrix
           kneigh = nPPx_map(ineigh,iatom)
! Assemble the global matrix 
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
!$omp atomic
             vnl(imu,inu,kneigh,iatom) =                                    &
     &        vnl(imu,inu,kneigh,iatom) + PPx(imu,inu)
            end do ! do imu
           end do ! do inu

          endif ! if(iatom .eq. jatom)
         enddo ! do ineigh
        enddo ! do iatom


! ****************************************************************************
!
! ASSEMBLE VNL ONTOP RIGHT CASE   <phi_i|Psi_j><Psi_j|phi_j>
! ****************************************************************************
! Loop over iatom
!$omp parallel do private (in1, matom, ineigh, jatom, mbeta, in2, cl)        &
!$omp&            private (mneigh_self, imu, inu, ncc, PPx, kneigh) 
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neighPP_self(iatom)
         in1 = imass(iatom)
 
! Loop over the neighbors of each iatom.
         do ineigh = 1, nPPn(iatom)   !  <==== loop over i's neighbors
          mbeta = nPP_b(ineigh,iatom)
          jatom = nPP_j(ineigh,iatom)
          in2 = imass(jatom)

! if r1 .ne. r2, then this is a case where the potential is located
! at one of the sites of a wavefunction (ontop case).
          if (iatom .eq. jatom .and. mbeta .eq. 0) then
 
! Do nothing here - special case. Interaction already calculated in atm case.
! sanity check
           if (nPP_self(iatom) .ne. ineigh) then
            write (*,*) ' Something really wrong in assemble_2c_PP.f90 '
            write (*,*) ' iatom, jatom, mbeta = ', iatom, jatom, mbeta
            write (*,*) ' neigh_self(iatom), ineigh = ',                     &
     &                    nPP_self(iatom), ineigh  
! FIXME stop is not allowed in an OpenMP parallelized loop
           stop
           end if ! if(neighPP_self)
          else ! if(iatom .eq. jatom)
 
! Now the second case. <i | V(j) | j>.
           call cl_value (in2, cl)
! Looking for <phi_j|Psi_j>, what is jneigh of jatom itself in the nPPx list 
           mneigh_self = nPP_self(jatom)
 
! Now we combine and sum:
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             PPx(imu,inu) = 0.0d0
             do ncc = 1, num_orbPP(in2)
              PPx(imu,inu) = PPx(imu,inu)                                    &
     &         + cl(ncc)*sVNL(imu,ncc,ineigh,iatom)                          &
     &                  *sVNL(inu,ncc,mneigh_self,jatom)
             end do
            end do
           end do

! Mapping to the global matrix
           kneigh = nPP_map(ineigh,iatom) 
! Assemble the global Matrix
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
!$omp atomic
             vnl(imu,inu,kneigh,iatom) =                                    &
     &        vnl(imu,inu,kneigh,iatom) + PPx(imu,inu)
            end do
           end do

          end if ! if(iatom .eq. jatom)

! ****************************************************************************
! END VNL ASSEMBLY.
! End loop over iatom and its neighbors - jatom.
         end do ! do ineigh
        end do ! do iatom

! Format Statements
! ===========================================================================

        return
        end subroutine assemble_2c_PP
