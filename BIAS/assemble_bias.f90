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


! assemble_Vbias.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the two-center interactions
! related to an applied bias voltage.
!
! Scheme
!
!       eV/2       zb0
!  -----------------\
!                    \ eV*(1/2-(z-zb0)/(zb1-zb0))
!                     \
!                      \--------------------------
!  0                    zb1          -eV/2      Lz
!  ---> z-axis
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: jelinekp@fzu.cz
! Program Declaration
! ===========================================================================
        subroutine assemble_Vbias (nprocs, iforce, iordern, ioff2c)
        use configuration
        use constants_fireball
        use dimensions
        use forces
        use interactions
        use neighbor_map
        use bias
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iforce
        integer, intent (in) :: iordern
        integer, intent (in) :: nprocs

! 24 different two-center interactions
        integer, intent (in), dimension (1:24) :: ioff2c

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
        integer jatom
        integer kforce
        integer matom
        integer mbeta
        integer mneigh_self
        integer my_proc
        integer natomsp

        real y
        real z1, z2, za
        real Vb
        real, dimension (3, 3, 3) :: deps
        real, dimension (3, 3) :: eps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: sighat
        real, dimension (numorb_max, numorb_max) :: sx
        real, dimension (3, numorb_max, numorb_max) :: spx

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Initialize interactions
        write (*,*)  ' Assembling bias_ham'
! necessary because these arrays will be partially computed and then summed
! over the procs
        Vbias_mat = 0.0d0
        Vbiasp_mat = 0.0d0

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
!!$omp parallel do private (ineigh, in1, in2, in3, jatom, matom, mbeta, kforce) &
!!$omp  private (r1, r2, r21, y, sighat, eps, deps, interaction, isorp) &
!!$omp  private (sx, spx, tx, tpx, bcna, bcnapx, bcnax, imu, inu)

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

! evaluate the bias potential
          z1 = ratom (3,iatom)
          z2 = ratom (3,jatom)
          za = (z1 + z2)/2.0d0
          if (za .gt. zb1) then
           Vb = -0.5d0*Vbias
          else if (za .gt. zb0) then
           Vb = (0.5d0 - (za - zb0)/(zb1-zb0))*Vbias
          else
           Vb = 0.5d0*Vbias
          endif
          write (100,*) 'atms:',iatom, jatom
          write (100,*) 'z* :',z1,z2,za
          write (100,*) 'Vb :',Vb, Vbias
          if (iatom .eq. jatom) write (120,*) ratom(3,iatom),Vb
! ****************************************************************************
!
! CALL DOSCENTROS AND GET S (Overlap matrix)
! ****************************************************************************
! The derivatives are tpx and spx, where p means derivative and x means crytal
! coordinates. The derivatives are:  d/d(b1) <phi(r-b1) ! O ! phi(r-b2)>
! where O = 1 (for overlap), or T.  The derivative is a vector in crystal
! coordinates and is stored in tp(ix,imu,inu,iatom,ineigh). The subroutine
! returns the derivative for just that one value of iatom and ineigh, and the
! result is returned in the arguement list, tpx(3,4,4) or spx(3,4,4).  The
! derivatives are computed only if iforce = 1 !
! For these interactions, there are no subtypes and isorp = 0
          isorp = 0
          interaction = 1
          in3 = in2
          call doscentros (interaction, isorp, iforce, in1, in2, in3, y,&
     &                     eps, deps, sx, spx)

! Write s and t to appropriate arrays
          do inu = 1, num_orb(in2)
           do imu = 1, num_orb(in1)
            Vbias_mat(imu,inu,ineigh,iatom) = sx(imu,inu)*Vb

            if (iforce .eq. 1) then
             Vbiasp_mat(:,imu,inu,ineigh,iatom) = spx(:,imu,inu)*Vb
            end if
           end do
          end do

! ****************************************************************************
! End loop over iatom and its neighbors - jatom.
         end do  ! do ineigh
        end do  ! do i

! Format Statements
! ===========================================================================
100     format(9f8.4)

        return
        end subroutine assemble_Vbias
