! copyright info:
!
!                             @Copyright 2002
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega

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

 
! assemble_G_S.f90
! Program Description
! ===========================================================================
!       This routine assembles the vector < G mu | nu > , where G stands for
! Gradient, in the variable gover(ix,imu,inu,ineigh,iatom).
! (G wrt iatom-position)
! The other gradient (with respect to the position of the second atom) 
! is calculated as
! < mu | G nu > = - gover
! (G wrt jatom-position)
!
!
! ===========================================================================
! Code written by:
! Jose Ortega Mateo
! Departmento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
!
! Program Declaration
! ===========================================================================
        subroutine assemble_G_S (nprocs, impi)
        use nonadiabatic
        use configuration
        use constants_fireball
        use density
        use dimensions
        use forces
        use interactions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: impi
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
        integer ix
        integer jatom
        integer iforce
        integer matom
        integer mbeta
        integer my_proc
        integer natomsp
 
        real muxc
        real sumS
        real sumT
        real y
 
        real, dimension (numorb_max, numorb_max) :: sx
        real, dimension (3, numorb_max, numorb_max) :: spx
        real, dimension (3, 3) :: eps
        real, dimension (3, 3, 3) :: deps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: sighat

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ===========================================================================
! Initialize to zero.
        gover = 0.0d0

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
 
! Loop over the atoms in the central cell.
!!$omp parallel do private (matom, r1, in1, ineigh, mbeta, jatom, r2, in2)    &
!!$omp&  private (r21, y, sighat, eps, deps, ix, sumT, sumS, inu, imu, muxc)  &
!!$omp&  private (isorp, kforce, interaction, in3, bcnax, bcnapx)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neigh_self(iatom)
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)
 
! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)     ! <==== loop 2 over iatom's neighbors
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)
 
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
          else
           sighat(:) = r21(:)/y
          end if
 
          call epsilon (r2, sighat, eps)
          call deps2cent (r1, r2, eps, deps)
!  
! ****************************************************************************
!
!  ASSEMBLE GRADIENT-MATRIX ELEMENTS BETWEEN BASIS ORBITALS
! ****************************************************************************
! The derivative of the overlap is spx, where p means derivative and x means crytal
! coordinates. The derivatives are:  d/d(b1) <phi(r-b1) ! phi(r-b2)>
!
! gover is   < d/d(b1) phi(r-b1) ! phi(r-b2)> ,
! which is a vector in crystal
! coordinates and is stored in gover(ix,imu,inu,ineigh,iatom).
! For the calculation of
! < phi(r-b1) ! d/d(b2) phi(r-b2)> 
! we will use 
! < phi(r-b1) ! d/d(b2) phi(r-b2)> = - < d/d(b1) phi(r-b1) ! phi(r-b2)> 
!
! First, assemble one-center terms (calculated before)
!
         if (iatom .eq. jatom .and. mbeta .eq. 0) then
! 1-C : ineigh = matom !
! JOM-info careful with the sign of gover1c
           call build_gover1c(in1)
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             gover(:,imu,inu,matom,iatom) = gover1c(:,imu,inu)
            end do
           end do
         else
! 2-C : ineigh 
!
! For these interactions, there are no subtypes and isorp = 0
! The derivatives are computed only if iforce = 1 !
          iforce = 1
          isorp = 0
          interaction = 1
          in3 = in2
          call doscentros (interaction, isorp, iforce, in1, in2, in3, y, &
     &                     eps, deps, sx, spx)

! write gover to the appropriate array

           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             gover(:,imu,inu,ineigh,iatom) = spx(:,imu,inu)
            end do
           end do
!
! end 1C-2C if
         end if
! ****************************************************************************
! JOM-test
!          write(*,*)'iatom,ineigh=',iatom,ineigh
!          do inu = 1, num_orb(in2)
!           do imu = 1, num_orb(in1)
!          write(*,100)imu,inu,gover(:,imu,inu,ineigh,iatom)
!           end do
!          end do
! ****************************************************************************
! End loop over iatom and its neighbors - jatom.
         end do
        end do

! Format Statements
100     format('gover=',2i4,3f8.4)
! ===========================================================================
        return
      end subroutine assemble_G_S
