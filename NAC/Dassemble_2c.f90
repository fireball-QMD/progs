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

 
! Dassemble_2c.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the two-center and degenerate
! two-center interactions.
!
! ===========================================================================
! JOM : adapted to also calculate the gradient of the Hamiltonian
! G < mu | H | nu >, 2C-part
! gh_2c (ix,imu,inu,ineigh,iatom)
! gh_2c is the gradient wrt to iatom
! for the gradient wrt to jatom we will use  
! (- gh_2c (ix,imu,inu,ineigh,iatom) )
! gh_atm (ix,imu,inu,ineigh,iatom) is for the ATOM special case
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658 
! Provo, UT 841184602-4658
! FAX 801-422-2265
! Office telephone 801-422-7444
! ===========================================================================
!
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_2c (nprocs, impi, igauss)
        use configuration
        use constants_fireball
        use density
        use dimensions
        use forces
        use interactions
        use neighbor_map
        use nonadiabatic
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: igauss
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
        integer kforce
        integer matom
        integer mbeta
        integer my_proc
        integer natomsp
 
        real muxc
        real sumS
        real sumT
        real y
 
        real, dimension (numorb_max, numorb_max) :: bcnax
        real, dimension (3, numorb_max, numorb_max) :: bcnapx
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
! Initialize forces to zero.
        fana = 0.0d0
        fotna = 0.0d0
        ft = 0.0d0
        fro = 0.0d0
        if (igauss .eq. 1) fxcro = 0.0d0
! JOM : Initialize gh_2c to zero.
        gh_2c = 0.0d0
        gh_atm = 0.0d0

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
 
 
! ****************************************************************************
!
!  ASSEMBLE S AND T FORCES
! ****************************************************************************
! The derivatives are tpx and spx, where p means derivative and x means crytal
! coordinates. The derivatives are:  d/d(b1) <phi(r-b1) ! O ! phi(r-b2)>
! where O = 1 (for overlap), or T.  The derivative is a vector in crystal
! coordinates and is stored in tp(ix,imu,inu,iatom,ineigh). The subroutine
! returns the derivative for just that one value of iatom and ineigh, and the
! result is returned in the arguement list, tpx(3,4,4) or spx(3,4,4).  The
! derivatives are computed only if iforce = 1 !
!
! Overlap and Kinetic.
! No need to call doscentros to get spx or tpx, because it was calculated in
! assemble_2c.f.
! Now do kinetic and repulsive overlap forces.
          do ix = 1, 3
           sumT = 0.0d0
           sumS = 0.0d0
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             sumT = sumT                                                     &
     &             + rho(imu,inu,ineigh,iatom)*tp_mat(ix,imu,inu,ineigh,iatom)
             sumS = sumS                                                     &
     &             + cape(imu,inu,ineigh,iatom)*sp_mat(ix,imu,inu,ineigh,iatom)
            end do
           end do

! Now add sum to appropriate force term. see notes "the total band structure
! force", ofs 9/14/88. Also see notes "total force due to fr offsites"
!                                 and "total force due to fs offsites" 11/10/88
! The (-1.d0) makes it "force-like".
!!$omp critical
! Direct terms.
           ft(ix,iatom) = ft(ix,iatom) + (-1.0d0)*sumT
           fro(ix,iatom) = fro(ix,iatom) +        sumS

! Cross terms.
           ft (ix,jatom) = ft (ix,jatom) - (-1.0d0)*sumT
           fro(ix,jatom) = fro(ix,jatom) -          sumS
!!$omp end critical
          end do ! do ix
 
! JOM : add gradient of kinetic to gh_2c
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             gh_2c(:,imu,inu,ineigh,iatom) = gh_2c(:,imu,inu,ineigh,iatom) &
     &                     +   tp_mat(:,imu,inu,ineigh,iatom)
            end do
           end do
! JOM-end : gradient of kinetic to gh_2c



! Gaussian approximation to three-center exchange-correlation interactions.
          if (igauss .eq. 1) then
           do ix = 1, 3
            do inu = 1, num_orb(in2)
             do imu = 1, num_orb(in1) 
              muxc = - vxc_3c(imu,inu,ineigh,iatom)                          &
     &               + nuxc_3c(imu,inu,ineigh,iatom)                         &
     &                 *bar_density_2c(imu,inu,ineigh,iatom)                 &
     &               + nuxc_total(imu,inu,ineigh,iatom)                      &
     &                 *bar_density_3c(imu,inu,ineigh,iatom) 
     
              fxcro(ix,ineigh,iatom) = fxcro(ix,ineigh,iatom)                &
     &         + muxc*rho(imu,inu,ineigh,iatom)*sp_mat(ix,imu,inu,ineigh,iatom) 
             end do 
            end do 
           end do 
          end if 

 
! ****************************************************************************
!
! ASSEMBLE NEUTRAL ATOM FORCE FOR ATOM CASE
! ****************************************************************************
! The vna 2 centers are: ontop (L), ontop (R), and atm.
! First, do vna_atom case. Here we compute <i | v(j) | i> matrix elements.
 
! If r1 = r2, then this is a case where the two wavefunctions are at the
! same site, but the potential vna is at a different site (atm case).
! The derivative wrt the "atom r1" position (not the NA position) are
! stored in bcnapx.
! Form the "force-like" derivative of the atom terms for NA,
! or -(d/dr1) <phi(mu,r-r1)!h(r-ratm)!phi(nu,r-r1)>.
          isorp = 0
          kforce = 1           ! calculate forces here
          interaction = 4
          in3 = in1
          call doscentros (interaction, isorp, kforce, in1, in2, in3, y,    &
    &                      eps, deps, bcnax, bcnapx)
 
! Note that the loop below involves num_orb(in1) ONLY. Why?
! Because the potential is somewhere else (or even at iatom), but we are
! computing the vna_atom term, i.e. < phi(i) | v | phi(i) > but V=v(j) )
! interactions.
! Notice the explicit - sign, this makes it force like.
          do inu = 1, num_orb(in3)
           do imu = 1, num_orb(in1)
            do ix = 1, 3
             fana(ix,ineigh,iatom) = fana(ix,ineigh,iatom)                   &
     &        - rho(imu,inu,matom,iatom)*bcnapx(ix,imu,inu)*eq2
            end do
           end do
          end do

! JOM : add gradient of vna_atom to gh_atm
! in this case ineigh = matom !
! JOM-info : notice that imu, inu both correspond to iatom, while ineigh
! is the other atom!!!
           do inu = 1, num_orb(in1)
            do imu = 1, num_orb(in1)
             gh_atm(:,imu,inu,ineigh,iatom) = gh_atm(:,imu,inu,ineigh,iatom) &
     &                       + bcnapx(:,imu,inu)*eq2
            end do
           end do
! JOM-end : gradient of vna_atom to gh_2c
 
1001      continue
! ****************************************************************************
!
! ASSEMBLE NEUTRAL ATOM FORCE FOR ONTOP CASE
! ****************************************************************************
! Ontop case.
! We have Ontop Left, and Ontop Right.
! Left is <1|V(1)|2> and Right is <1|V(2)|2>
! Check to make sure we are not doing <1|V(1)|1>. This is done in the atm case.
          if (iatom .eq. jatom .and. mbeta .eq. 0) then
 
! Do nothing here - special case.
 
          else
 
! Form the Left ontop force. Use the derivatives, since the derivative is
! with respect to d/d(ratom) when the atom is ontop atom 1.
! Note that we only compute ontop left. That is because we do cross terms.
! If we were to do both ontop left and ontop right, then we would get
! double counting in the forces.
 
! Neutral atom piece
! Form the Left ontop force. Use the derivatives, since the derivative is
! with respect to d/d(ratom) when the atom is ontop atom 1.
           isorp = 0
           interaction = 2
          kforce = 1           ! calculate forces here
           in3 = in2
           call doscentros (interaction, isorp, kforce, in1, in1, in3, y,   &
     &                      eps, deps, bcnax, bcnapx)
 
! Notice the explicit - sign which makes f force like.
           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
             do ix = 1, 3
              fotna(ix,ineigh,iatom) = fotna(ix,ineigh,iatom)                &
     &         - rho(imu,inu,ineigh,iatom)*bcnapx(ix,imu,inu)*eq2
             end do
            end do
           end do

! JOM : add gradient of vna_ontop_left to gh_2c

           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             gh_2c(:,imu,inu,ineigh,iatom) = gh_2c(:,imu,inu,ineigh,iatom) &
     &                       + bcnapx(:,imu,inu)*eq2
            end do
           end do

! JOM : add gradient of vna_ontop_right to gh_2c

! JOM-info : Now, add gradient of vna_ontop_right to gh_2c
! JOM-info : Here we need the ontop_right case, but for forces (see
! above) the contribution from ontop_right is added later in
! assemble_F.f90 by multiplying by 2.0d0
! Form the Right ontop force. Use the derivatives, since the derivative is
! with respect to d/d(ratom) when the atom is ontop atom 1.
           isorp = 0
           interaction = 3
          kforce = 1           ! calculate forces here
           in3 = in2
           call doscentros (interaction, isorp, kforce, in1, in2, in3, y,   &
     &                      eps, deps, bcnax, bcnapx)
 
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             gh_2c(:,imu,inu,ineigh,iatom) = gh_2c(:,imu,inu,ineigh,iatom) &
     &                       + bcnapx(:,imu,inu)*eq2
            end do
           end do

! JOM-end : gradient of vna_ontop to gh_2c
 
          end if
! ****************************************************************************
! JOM-test
1000      continue
!           write(*,*)'**********GH_2C**************'
            write(*,*)'iatom, ineigh', iatom, ineigh
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
            write(*,100)imu,inu,gh_2c(:,imu,inu,ineigh,iatom)
            end do
           end do
!           write(*,*)'***********GH_ATM************'
!           write(*,*)'iatom, ineigh', iatom, ineigh
!          do inu = 1, num_orb(in2)
!           do imu = 1, num_orb(in1)
!           write(*,*)'gh_atm=', gh_atm(1,imu,inu,ineigh,iatom), imu,inu
!           end do
!          end do
! End loop over iatom and its neighbors - jatom.
         end do
        end do

! Format Statements
100     format('gh_2c='2i4,3f8.4)
! ===========================================================================
        return
      end subroutine Dassemble_2c
