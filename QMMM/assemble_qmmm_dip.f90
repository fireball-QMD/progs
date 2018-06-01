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

 
! assemble_qmmm_dip.f90
! Program Description
! ===========================================================================
!       This routine calculated the long range matrix elements for
!       the qm/mm embedding with the true dipole term
!
! ===========================================================================
! Code originally written by:
! Jesús I. Mendieta-Moreno
!
! Program Declaration
! ===========================================================================
        subroutine assemble_qmmm_dip (nprocs, iordern) 
        use charges
        use configuration
        use constants_fireball
        use dimensions
        use interactions
        use neighbor_map
        use energy
        !use qmmm_module, only : qmmm_struct, qmmm_nml
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
        integer inu
        integer in1
        integer in2
        integer in3
        integer ineigh
        integer issh
        integer jatom
        integer katom
        integer mbeta
        integer my_proc
        integer natomsp
 
        real dij
        real dterm
        real sterm
        real dq3
        real dq4
 
        real, dimension (3) :: rna
        real, dimension (3) :: rnabc
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real :: x
        real :: emnpl

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE


! Procedure
! ===========================================================================
! Initialize interactions to zero.

        ewaldqmmm = 0.0d0
        emnpl = 0.0d0

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

        do iatom = iatomstart, iatomstart - 1 + natomsp
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)

         do ineigh = 1, neighn(iatom)
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          r2(:) = ratom(:,jatom) !+ xl(:,mbeta)
          in2 = imass(jatom)

          do katom = 1, qmmm_struct%qm_mm_pairs
           rna(1) = qmmm_struct%qm_xcrd(1,katom)
           rna(2) = qmmm_struct%qm_xcrd(2,katom)
           rna(3) = qmmm_struct%qm_xcrd(3,katom)
           dq3 = - qmmm_struct%qm_xcrd(4,katom) ! charge in amber have opposite sign
           r21(:) = r2(:) - r1(:)
           rnabc(:) = rna(:) - (r1(:) + r21(:)/2.0d0)
           x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)

             sterm = s_mat(imu,inu,ineigh,iatom)
            
             dterm = (dipc(1,imu,inu,ineigh,iatom)*rnabc(1)    &
              &     + dipc(2,imu,inu,ineigh,iatom)*rnabc(2)    &
              &     + dipc(3,imu,inu,ineigh,iatom)*rnabc(3))

             emnpl = dq3*sterm/x + dq3*dterm/(x*x*x)

             ewaldqmmm(imu,inu,ineigh,iatom) = ewaldqmmm(imu,inu,ineigh,iatom)  &
             &                               + emnpl*eq2

            end do !end do imu = 1, num_orb(in1)
           end do  ! end do inu = 1, num_orb(in2)

          end do    ! end do katom: mm atom

         end do  ! end do ineigh = 1, neighn(iatom)
        end do   ! end do iatom = 1, natoms

        eqmmm = 0.0d0
        do iatom = iatomstart, iatomstart - 1 + natomsp
         in3 = imass(iatom)
         dq4 = 0.0d0
         do issh = 1, nssh(in3)
          dq4 = dq4  + Qneutral(issh,in3)
         end do
         do katom = 1, qmmm_struct%qm_mm_pairs
          dij = sqrt ((qmmm_struct%qm_xcrd(1,katom)-ratom(1,iatom))*(qmmm_struct%qm_xcrd(1,katom)-ratom(1,iatom)) + &
                      (qmmm_struct%qm_xcrd(2,katom)-ratom(2,iatom))*(qmmm_struct%qm_xcrd(2,katom)-ratom(2,iatom)) + & 
                      (qmmm_struct%qm_xcrd(3,katom)-ratom(3,iatom))*(qmmm_struct%qm_xcrd(3,katom)-ratom(3,iatom)) )

          eqmmm = eqmmm + (dq4*qmmm_struct%qm_xcrd(4,katom) / dij)*eq2
         end do
        end do


! Format Statements
! ===========================================================================
        return
        end
