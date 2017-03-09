! copyright info:
!
!                             @Copyright 2004
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Universidad de Madrid - Pavel Jelinek
! Brigham Young University - Hao Wang

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio State University - Dave Drabold
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

 
! assemble_olsxc_2c.f90
! Program Description
! ===========================================================================
!       This routine assembles the two & three-center exchange-correlation
! for the average density approximation. 
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
       subroutine assemble_snxc_off (natoms, nprocs, my_proc, iordern, itheory)
        use charges
        use density
        use dimensions
        use interactions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iordern
        integer, intent (in) :: itheory
        integer, intent (in) :: my_proc
        integer, intent (in) :: natoms
        integer, intent (in) :: nprocs

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iatomstart
        integer imu
        integer in1
        integer in2
        integer in3
        integer ineigh
        integer inu
        integer jatom
        integer matom
        integer jbeta
        integer natomsp

        real, dimension (numorb_max, numorb_max) :: bcxcx
        real, dimension (numorb_max, numorb_max) :: denmx
        real, dimension (numorb_max, numorb_max) :: sx
      
! Procedure
! ===========================================================================
! Determine which atoms are assigned to this processor.
        if (iordern .eq. 1) then
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
!!$omp parallel do private (matom, in1, ineigh, jbeta, jatom, in2, inu, imu)  &
!!$omp&  private (denmx, sx, bcxcx)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neigh_self(iatom)
         in1 = imass(iatom)

! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)       ! <==== loop over i's neighbors
          jbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom)
 
 
! ****************************************************************************
!
! CALL DOSCENTROS AND GET VXC FOR ATM CASE - AVERAGE DENSITY APPROXIMATION 
! ****************************************************************************

!   denmx    .... <i|n|j> = <i|n_i|j> + <i|n_j|j> + S_{i.ne.j.ne.k}<i|n_k|j>
!   sx       ..... <i|j>
!   
          if (iatom .eq. jatom .and. jbeta .eq. 0) then

! skip atom case i-site is identical to j-site
                 
          else ! i.ne.j
           ! restore density and overlap matrices
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             denmx(imu,inu) = rho_off(imu,inu,ineigh,iatom)
             sx(imu,inu) = s_mat(imu,inu,ineigh,iatom)
            end do
           end do
           call build_snxc_off (in1, in2, denmx, sx, ineigh, iatom, bcxcx)

! Now complete 'non-diagonal' terms <i|V(n)|j>
           if (itheory .eq. 1) then
            do inu = 1, num_orb(in2)
             do imu = 1, num_orb(in1)
              vxc_ca(imu,inu,ineigh,iatom) =                                    &
     &         vxc_ca(imu,inu,ineigh,iatom) + bcxcx(imu,inu)
             end do
            end do
           else
            do inu = 1, num_orb(in2)
             do imu = 1, num_orb(in1)
              vxc(imu,inu,ineigh,iatom) =                                 &
     &         vxc(imu,inu,ineigh,iatom) + bcxcx(imu,inu)
             end do
            end do
           end if

! End if for r1 .ne. r2 case
          end if
! ****************************************************************************
! End loop over iatom and its neighbors - jatom.
         end do
        end do

! Format Statements
! ===========================================================================
        return
        end subroutine assemble_snxc_off
