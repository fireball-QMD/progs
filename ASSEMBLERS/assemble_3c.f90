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

 
! assemble_3c.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the true three-center interactions.
! These are true three-center in that iatom .ne. jatom .ne. katom.
! The matrix elements look like <Psi(1)|V(3)|Psi(2)>.
!
! The bulk of the work is done in trescentros.f. This program assembles the
! results.
!
! A third party term is when we consider the NA (or etc.) to be at the origin
! and take the matrix element between a pair of neighbors, neither of which is
! the NA (or etc.), and the pair is truly a pair, and not an atom.
! the results are stored in: f3na(ix,i), f3xc(ix,i), etc.
!
! The generic form of the matrix element is
! v = <mu(r-(r1+ratm)) ! v(r-ratm) ! nu(r-(r2+ratm))>.
!
! See notes
! "general theory of third party terms" for the definition on p. 2 of
! these three terms.
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
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_3c (nprocs, iordern, igauss, itheory_xc)
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
        integer, intent (in) :: itheory_xc
        integer, intent (in) :: igauss
        integer, intent (in) :: iordern
        integer, intent (in) :: nprocs
 
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ialp
        integer iatom
        integer iatomstart
        integer ibeta
        integer ierror
        integer imu
        integer in1
        integer in2
        integer indna
        integer ineigh
        integer interaction
        integer inu
        integer isorp
        integer issh
        integer jatom
        integer jbeta
        integer jssh
        integer mneigh
        integer my_proc
        integer natomsp
        integer jneigh       

        real cost
        real distance_13
        real distance_23
        real dstn_temp
        real stn_temp1
        real stn_temp2
        real x
        real y

        real, dimension (numorb_max, numorb_max) :: bcnax
        real, dimension (3, 3) :: eps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: r31
        real, dimension (3) :: r32
        real, dimension (3) :: r13
        real, dimension (3) :: r23
        real, dimension (3) :: rhat
        real, dimension (3) :: rna
        real, dimension (3) :: rnabc
        real, dimension (3) :: sighat
 
! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ========================================================================
! The array vna will not be initialized to zero here.
! Presumably, the two-center interactions have already been calculated.
! As such, at this point and time these arrays should not be zero.
 
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

! Choose atom ialp in the central cell. This is the atom whose position
! we take the derivative, and is the atom who has the the neutral atom
! potential.
! Loop over the atoms in the central cell.
!$omp parallel do private (iatom, ibeta, in1, in2, indna, interaction, isorp)&
!$omp&            private (jatom, jbeta, mneigh, cost, distance_13)          &
!$omp&            private (distance_23, dstn_temp, stn_temp1, stn_temp2, x)  &
!$omp&            private (y, bcnax, eps, r1, r2, r21, r31, r32, r13, r23)   &
!$omp&            private (rhat, rna, rnabc, sighat)
        do ialp = iatomstart, iatomstart - 1 + natomsp
         rna(:) = ratom(:,ialp)
         indna = imass(ialp)
 
! Loop over the neighbors of each ialp.
! Now look at all common neighbor pairs which were figured out in main.
! The common neighbors were figured out in common_neighbors.f90
         do ineigh = 1, neigh_comn(ialp)
          mneigh = neigh_comm(ineigh,ialp)
 
! The second atom (jatom) is the mneigh'th neighbor of iatom.
          if (mneigh .ne. 0) then
           iatom = neigh_comj(1,ineigh,ialp)
           ibeta = neigh_comb(1,ineigh,ialp)
           r1(:) = ratom(:,iatom) + xl(:,ibeta)
           in1 = imass(iatom)

           jatom = neigh_comj(2,ineigh,ialp)
           jbeta = neigh_comb(2,ineigh,ialp)
           r2(:) = ratom(:,jatom) + xl(:,jbeta)
           in2 = imass(jatom)
           jneigh = neigh_back(iatom,mneigh)

 
! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
           r21(:) = r2(:) - r1(:)
           y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
 
           r13(:) = r1(:) - rna(:)
           r23(:) = r2(:) - rna(:)
           distance_13 = sqrt(r13(1)**2 + r13(2)**2 + r13(3)**2)
           distance_23 = sqrt(r23(1)**2 + r23(2)**2 + r23(3)**2)

! Find the unit vector in sigma direction.
           if (y .lt. 1.0d-05) then
            sighat(1) = 0.0d0
            sighat(2) = 0.0d0
            sighat(3) = 1.0d0
            write (*,*) ' There is an error here in assemble_3c.f '
            write (*,*) ' r1 = r2!!!! BAD!!!! '
           else
            sighat(:) = r21(:)/y
           end if
 
! Find rnabc = vector pointing from center of bondcharge to the neutral atom.
! The center of the bondcharge is at r1 + r21/2.  This gives us the distance
! dnabc (x value in 2D grid).
           rnabc(:) = rna(:) - (r1(:) + r21(:)/2.0d0)
           x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)
 
! Find the unit vector in rnabc direction.
           if (x .lt. 1.0d-05) then
            rhat(1) = 0.0d0
            rhat(2) = 0.0d0
            rhat(3) = 0.0d0
           else
            rhat(:) = rnabc(:)/x
           end if
           cost = sighat(1)*rhat(1) + sighat(2)*rhat(2) + sighat(3)*rhat(3)
           call epsilon (rhat, sighat, eps)
 
! For now we just do the neutral atom interactions.
! Charged atom interactions are assembled in assemble_ca_3c.f
! So set isorp = 0 within this subroutine.
!
!              interaction    subtypes     index
!
!      bcna         1           0..9(max)   1..10
 
! ****************************************************************************
!
! PERFORM ACTUAL CALCULATIONS
! CALL TRESCENTROS FOR NEUTRAL ATOM PIECE
! ****************************************************************************
           isorp = 0
           interaction = 1
           if (igauss .eq. 0) then
            call trescentros (interaction, isorp, isorpmax, in1, in2,   &
     &                        indna, x, y, cost, eps, bcnax, nspecies)
           else 
            call trescentrosG_VNA (in1, in2, indna, x, y, cost, eps,    &
     &                             bcnax, rcutoff)
           end if
 
! Add this piece for iatom, jatom, and katom into the total - bcna and bcca
!$omp critical (vna3)
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             vna(imu,inu,mneigh,iatom) =                                &
     &        vna(imu,inu,mneigh,iatom) + bcnax(imu,inu)*eq2
            !Symmetrize Hamiltonian (April 2018): jneigh is the
            !back_neigh:
              vna(inu,imu,jneigh,jatom) = vna(imu,inu,mneigh,iatom)
            end do
           end do
!$omp end critical (vna3)
! ****************************************************************************
! End loop over ialp and its common neighbors.
          end if
         end do ! do ineigh
        end do ! do ialp

! Format Statements
! ===========================================================================

        return
        end subroutine assemble_3c
