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

 
! assemble_lr_dip.f90
! Program Description
! ===========================================================================
!       This routine calculated the long long-range ewald matrix elements -
! ewaldlr(mu,nu,ineigh,iatom).
!
! ===========================================================================
! Code originally written by:
! Otto F. Sankey
! Campus Box 1504
! Department of Physics
! Arizona State University
! Tempe, AZ 85287-1504
! (602) 965-4334 (office)      email: otto.sankey@asu.edu
 
! Code rewritten by:
! James P. Lewis
!
! code rewritten by Jesús I. Mendieta-Moreno & Diego Soler
! for the dipole long-range theory
!
! see equation (11) in PRB 52, 1618 (1995)
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_lr_dip (nprocs, iordern) 
        use charges
        use configuration
        use constants_fireball
        use dimensions
        use interactions
        use neighbor_map
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
        integer ineigh
        integer issh
        integer jatom
        integer mbeta
        integer my_proc
        integer natomsp
        integer ialp
        integer inalp

        real dist13
        real dist23
        real dq3
        real dterm
        real sterm
        real x
    
        
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: rna
        real, dimension (3) :: r13
        real, dimension (3) :: r23
        real, dimension (3) :: r21
        real, dimension (3) :: rnabc
        real, dimension (natoms) :: sub_ewald
        real, dimension (numorb_max, numorb_max) :: emnpl

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ===========================================================================
! Initialize interactions to zero.
        ewaldlr = 0.0d0
!          open( unit = 788, file = 'ewaldlr', status = 'unknown')

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

!  =========  NEW DIP_XY ============
!Here we compute long-range terms. The 2-center atm case is included,
!because the loop includes the case ineigh=matom, in which the neighbor
!ineigh is the atom iatom itself. Therefore, we need to compute ewaldsr
!in the atm case in assemble_ca_2c.f90. Also notice that this triple loop
!includes the case in which the three atoms are pairwise neighbors, so
!we'll need to compute the ewaldsr in assemble_ca_3c.f90
!
         do iatom = 1, natoms
          r1(:) = ratom(:,iatom)
          in1 = imass(iatom)
          do ineigh = 1, neighn(iatom)
           mbeta = neigh_b(ineigh,iatom)
           jatom = neigh_j(ineigh,iatom)
           r2(:) = ratom(:,jatom) + xl(:,mbeta)
           do ialp = 1, natoms   !the ialp is the "distant" atom.
            rna(:) = ratom(:,ialp)
            inalp = imass(ialp)
            dq3 = 0.0d0
            do issh = 1, nssh(inalp)
             dq3 = dq3 + (Qin(issh,ialp) - Qneutral(issh,inalp))
            end do ! end do issh = 1m nssh(inalp)

            in2 = imass(jatom)
            r13=rna(:)-r1(:)
            r23=rna(:)-r2(:)
            dist13=sqrt(r13(1)*r13(1) + r13(2)*r13(2) + r13(3)*r13(3))
            dist23=sqrt(r23(1)*r23(1) + r23(2)*r23(2) + r23(3)*r23(3))
            if ((dist13 .lt. 1.0d-5) .or. (dist23 .lt. 1.0d-5)) then   
            else
             r21(:) = r2(:) - r1(:)
             rnabc(:) = rna(:) - (r1(:) + r21(:)/2.0d0)
             x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)

             do inu = 1, num_orb(in2)
              do imu = 1, num_orb(in1)

               sterm = s_mat(imu,inu,ineigh,iatom)
            
               dterm = (dipc(1,imu,inu,ineigh,iatom)*rnabc(1)    &
                &     + dipc(2,imu,inu,ineigh,iatom)*rnabc(2)    &
                &     + dipc(3,imu,inu,ineigh,iatom)*rnabc(3))

               emnpl(imu,inu) = dq3*sterm/x + dq3*dterm/(x*x*x)

               ewaldlr(imu,inu,ineigh,iatom) = ewaldlr(imu,inu,ineigh,iatom)  &
               &                             + emnpl(imu,inu)*eq2

             end do !end do imu = 1, num_orb(in1)
            end do  ! end do inu = 1, num_orb(in2)
           end if
          end do    ! end do ialp = 1, natoms



          end do  ! end do ineigh = 1, neighn(iatom)
         end do   ! end do iatom = 1, natoms
! ========== NEW DIP_XY ============

! Format Statements
! ===========================================================================
200     format (4i,f10.4)
        return
        end
