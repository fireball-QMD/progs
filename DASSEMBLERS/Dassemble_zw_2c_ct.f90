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

 
! Dassemble_zw_2c_ct.f90
! Program Description
! ===========================================================================
! ===========================================================================!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_zw_2c_ct (nprocs, iordern)
        use charges
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
        integer, intent (in) :: iordern
        integer, intent (in) :: nprocs
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iatomstart
        integer icount
        integer icount_sav
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
        integer ix
        integer jatom
        integer jcount
        integer jcount_sav
        integer jssh
        integer kforce
        integer matom, matom2
        integer mbeta
        integer my_proc
        integer natomsp
        integer :: count_l, count_l_ini, issh1, issh2
 
 
        real dq1
        real dq2
        real dstn_temp
        real dterm_1
        real dterm_2
        real dterm
        real dxn
        real rcutoff_j
        real rend
        real rend1
        real rend2
        real sterm_1
        real sterm_2
        real y
        real rcutoff_i
        real :: A,B
 
        real, dimension (numorb_max, numorb_max) :: bcca
        real, dimension (3, numorb_max, numorb_max) :: bccap
        real, dimension (3, numorb_max, numorb_max) :: bccapx
        real, dimension (numorb_max, numorb_max) :: bccax
        real, dimension (3,numorb_max, numorb_max) :: demnpl
        real, dimension (3, 3, 3) :: deps
        real, dimension (3, numorb_max, numorb_max) :: dewaldsr
        real, dimension (3) :: dA, dB
! JOM-JIMM
!        real, dimension (numorb_max, numorb_max) :: dstn1
!        real, dimension (numorb_max, numorb_max) :: dstn2
        real dstn1
        real dstn2
        real, dimension (numorb_max, numorb_max) :: emnpl
        real, dimension (3, 3) :: eps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: sighat
        real, dimension (3) :: spterm_1
        real, dimension (3) :: spterm_2
! JOM-JIMM
!        real, dimension (numorb_max, numorb_max) :: stn1
!        real, dimension (numorb_max, numorb_max) :: stn2
        real stn1
        real stn2 

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ===========================================================================
! Initialize the forces to zero.
      faxc_ca = 0.0d0
      fotxc_ca = 0.0d0
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



         !The relevant integrals g2nu and its derivatives g2nup have
         !been computed in assemble_zw_2c_ct in the first step of the
         !SCF loop

        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neigh_self(iatom)
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)
         dq1 = 0.0
         do issh = 1, nssh(in1)
          dq1 = dq1 + (Qin(issh,iatom) - Qneutral(issh,in1))
         end do
         do ineigh = 1, neighn(iatom)    ! <==== loop 2 over iatom's neighbors
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          matom2 = neigh_self(jatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)
          dq2 = 0.0d0
          do issh = 1, nssh(in2)
           dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
          end do
          r21(:) = r2(:) - r1(:)
          y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
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
! ATOM CASE
! ****************************************************************************
! Initialize bcca and bccap for the charge atom interactions
          bcca = 0.0d0
          bccap = 0.0d0
          kforce = 1
          interaction = 4
          in3 = in1
          do isorp = 1, nssh(in2)
           dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
           do imu = 1, num_orb(in1)
             !bcca(imu,inu) = bcca(imu,inu) + dxn*bccax(imu,inu)
             issh = orb2shell(imu,in1)
             bccap(:,imu,imu) = bccap(:,imu,imu) + &
             & g2nup(:,issh,isorp,ineigh,iatom)*dxn
           end do !end do imu
          end do !end do isorp
 
 
! Now write bccap to the charged atom force piece.
! Add dewaldsr to the charged force term.  This is added in, since the energy
! is subtracted.
! Notice the explicit - sign, this makes it force like.
          do inu = 1, num_orb(in3)
           do imu = 1, num_orb(in1)
            do ix = 1, 3
!!$omp atomic
              faxc_ca(ix,ineigh,iatom) = faxc_ca(ix,ineigh,iatom)                   &
       &       - rho(imu,inu,matom,iatom)*bccap(ix,imu,inu)           
            end do
           end do
          end do
 
! ****************************************************************************
!
!                               ONTOP CASE
! ****************************************************************************
          if (iatom .eq. jatom .and. mbeta .eq. 0) then
 
! Do nothing here - special case.
 
          else

 
! Now ontop Left.
! Note that we only compute ontop left. That is because we do cross terms.
! If we were to do both ontop left and ontop right, then we would get
! double counting in the forces.
 
! Initialize bccap for the charge atom interactions
           bccap = 0.0d0
 
! Charged atom piece - Left ontop
           interaction = 2
           in3 = in2
           do isorp = 1, nssh(in1)
            !call doscentros (interaction, isorp, kforce, in1, in1, in3, y,   &
            ! &                       eps, deps, bccax, bccapx)
            dxn = (Qin(isorp,iatom) - Qneutral(isorp,in1))
 
! Now add d/dr1  to bccapx.
            do imu = 1, num_orb(in1)
             do inu = 1, num_orb(in3)
              issh1=orb2shell(imu,in1)
              issh2=orb2shell(inu,in2)
              !Here is the split of the INTEGRAL gcab...--->
              A = 0.5*s_mat(imu,inu,ineigh,iatom)-dip(imu,inu,ineigh,iatom)/y   
              B = 0.5*s_mat(imu,inu,ineigh,iatom)+dip(imu,inu,ineigh,iatom)/y  
              dA(:) = 0.5*sp_mat(:,imu,inu,ineigh,iatom)-dipp(:,imu,inu,ineigh,iatom)/y- &
                    & dip(imu,inu,ineigh,iatom)*sighat(:)/(y*y) 
              dB(:) = 0.5*sp_mat(:,imu,inu,ineigh,iatom)+dipp(:,imu,inu,ineigh,iatom)/y+ &
                    &  dip(imu,inu,ineigh,iatom)*sighat(:)/(y*y)  
              !force:
              bccap(:,imu,inu) = bccap(:,imu,inu) + dxn*(dA(:)*g2nu(issh1,isorp,matom,iatom)+ &
              &dB(:)*g2nu(isorp,issh2,ineigh,iatom)+A*g2nup(:,issh1,isorp,matom,iatom)+ &
              & B*g2nup(:,isorp,issh2,ineigh,iatom))
             ! bccap(:,imu,inu) = bccap(:,imu,inu) + dxn*(dA(:)*g2nu(issh1,isorp,ineigh,iatom)+ &
             ! & dB(:)*g2nu(issh2,isorp,matom2,jatom)+A*g2nup(:,issh1,isorp,ineigh,iatom)+      &
             ! & B*g2nup(:,issh2,isorp,matom2,jatom))
             end do !end do inu
            end do !end do imu
           end do !end do isorp
 
! Notice the explicit - sign which makes f force like.
           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
             do ix = 1, 3
!!$omp atomic
             fotxc_ca(ix,ineigh,iatom) = fotxc_ca(ix,ineigh,iatom)                 &
     &        - rho(imu,inu,ineigh,iatom)*bccap(ix,imu,inu)             
             end do
            end do
           end do
          end if
 
! ****************************************************************************
! End loop over iatom and its neighbors - jatom.
         end do
        end do
! Format Statements
      !fotxc_ca = 0.0d0
      !faxc_ca = 0.0d0
! ===========================================================================
        return
        end
