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

 
! assemble_hxc_2c.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the two-center and degenerate
! two-center interactions for the Horsfield exchange-correlation.
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
! Program Declaration
! ===========================================================================
        subroutine assemble_hxc_2c (nprocs, Kscf, iordern,      &
     &                              itheory, igauss)
        use charges
        use configuration
        use dimensions
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
        integer interaction
        integer inu
        integer isorp
        integer issh
        integer jatom
        integer kforce
        integer matom
        integer mbeta
        integer my_proc
        integer natomsp
        integer itest
 
        real cost
        real dq1
        real dq2
        real y
 
        real, dimension (numorb_max, numorb_max) :: bcxcx
        real, dimension (3, numorb_max, numorb_max) :: bcxcpx
        real, dimension (4) :: dqfact          ! there are four different isorps
        real, dimension (3, 3) :: eps
        real, dimension (3, 3, 3) :: deps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: sighat

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ===========================================================================
! Initialize interactions
        if (Kscf .eq. 1) vxc = 0.0d0
        if (itheory .eq. 1) vxc_ca = 0.0d0
        if (igauss .eq. 1) density_2c = 0.0d0

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

! Loop over the atoms in the central cell.
!$omp parallel do private (in1, in2, in3, interaction, isorp, jatom, kforce) &
!$omp&            private (matom, mbeta, dq1, dq2, y, bcxcx, bcxcpx, dqfact) &
!$omp&            private (eps, deps, r1, r2, r21, sighat)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neigh_self(iatom)
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)

! Find the net charge on iatom
         dq1 = 0.0d0
         do issh = 1, nssh(in1)
          dq1 = dq1 + (Qin(issh,iatom) - Qneutral(issh,in1))
         end do
 
! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)          ! <==== loop over i's neighbors
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)

! Find the net charge on jatom
          dq2 = 0.0d0
          do issh = 1, nssh(in2)
           dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
          end do

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
 
 
! ****************************************************************************
!
! CALL DOSCENTROS AND GET VXC FOR ATM AND ONTOP CASE - NA CASE
! ****************************************************************************
          if (Kscf .eq. 1) then
           if (iatom .eq. jatom .and. mbeta .eq. 0) then
 
! Do nothing here - special case. Interaction already calculated in unocentros.f
 
           else

! This is the atom case for the exchange-correlation energy
            kforce = 0 
            isorp = 0
            interaction = 7
            in3 = in1
            call doscentros (interaction, isorp, kforce, in1, in2, in3, &
     &                       y, eps, deps, bcxcx, bcxcpx)
 
! Note: these loops are both over in1 indeed!  This is because the two
! wavefunctions are located on the same site, but the potential is located
! at a different site.
            do inu = 1, num_orb(in3)
               do imu = 1, num_orb(in1)
                  vxc(imu,inu,matom,iatom) =                            &
     &            vxc(imu,inu,matom,iatom) + bcxcx(imu,inu)
               end do
            end do

! This is the ontop case for the exchange-correlation energy
            isorp = 0
            interaction = 6
            in3 = in2
            call doscentros (interaction, isorp, kforce, in1, in2, in3, &
     &                       y, eps, deps, bcxcx, bcxcpx)
            do inu = 1, num_orb(in3)
             do imu = 1, num_orb(in1)
              vxc(imu,inu,ineigh,iatom) =                               &
     &        vxc(imu,inu,ineigh,iatom) + bcxcx(imu,inu)
             end do
            end do

            end if   ! end if of iatom = jatom
           end if    ! end if of Kscf = 1
! ****************************************************************************
!
! CALL DOSCENTROS AND GET VXC FOR ATM AND ONTOP CASE - CA CASE
! ****************************************************************************
          if (itheory .eq. 1) then

           if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. Interaction already calculated in unocentros.f

           else

! Initialize dqfact; dq1 and dq2 are weighted by dq(in1) and dq(in2) 
            dqfact(1) = -dq1/(2.0d0*dq(in1))
            dqfact(2) =  dq1/(2.0d0*dq(in1)) 
            dqfact(3) = -dq2/(2.0d0*dq(in2)) 
            dqfact(4) =  dq2/(2.0d0*dq(in2)) 
           
! This is the atom case for the exchange-correlation energy 
! Note - the summation over isorp leads to a first order expansion for vxc_ca 
! with respect to the charges dq1 and dq2.  
            interaction = 7
            in3 = in1 
            do isorp = 1, 4 
               call doscentros (interaction, isorp, kforce, in1, in2,   &
     &                          in3, y,  eps, deps, bcxcx, bcxcpx) 
     
! Note: these loops are both over in1 indeed!  This is because the two 
! wavefunctions are located on the same site, but the potential is located at 
! a different site.  
               do inu = 1, num_orb(in3) 
                do imu = 1, num_orb(in1) 
                vxc_ca(imu,inu,matom,iatom) =                           &
     &          vxc_ca(imu,inu,matom,iatom)+dqfact(isorp)*bcxcx(imu,inu)
                end do
               end do
           end do

! This is the ontop case for the exchange-correlation energy 
           interaction = 6 
           in3 = in2 
           do isorp = 1, 4 
            call doscentros (interaction, isorp, kforce, in1, in2, in3, &
     &                       y, eps, deps, bcxcx, bcxcpx) 
 
            do inu = 1, num_orb(in3) 
             do imu = 1, num_orb(in1) 
              vxc_ca(imu,inu,ineigh,iatom) =                            &
     &        vxc_ca(imu,inu,ineigh,iatom)+dqfact(isorp)*bcxcx(imu,inu) 
             end do 
            end do 
            end do 
           end if
          end if 
 
! ****************************************************************************
! End loop over iatom and its neighbors - jatom.
         end do
        end do

! Format Statements
! ===========================================================================
 
        return
        end subroutine assemble_hxc_2c
        
