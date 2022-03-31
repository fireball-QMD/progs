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

 
! Dassemble_3c.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the true three-center interactions.
! These are true three-center in that iatom .ne. jatom .ne. katom.
! The matrix elements look like <Psi(1)|V(3)|Psi(2)>.
!
! The bulk of the work is done in Dtrescentros.f. This program assembles the
! results.
!
! A third party term is when we consider the NA (or etc.) to be at the origin
! and take the matrix element between a pair of neighbors, neither of which is
! the NA (or etc.), and the pair is truly a pair, and not an atom.
! the results are stored in: f3na(ix,i) and f3nl(ix,i)
!
! The generic form of the matrix element is
! v = <mu(r-(r1+ratm)) ! v(r-ratm) ! nu(r-(r2+ratm))>.
! We want the derivative with respect to ratm keeping the position of
! orbitals mu, nu fixed (const r1+ratm).
!
! The derivative f3 looks like for the neutral atom case:
! f3na = - sum (all neighbors of atom alpha at (li,bi) but bi.ne.balph)
!    * sum (all neighbors m of (li,bi), and not having b value balph)
!    * rho(mu,nu,i,m)* deriv wrt balpha <i! v(balph) !j>.
! Note the minus sign to make it "force-like".
!
! This program gets fa, fb, and fc pieces from Dtrescentros.
!
! See notes
! "general theory of third party terms" for the definition on p. 2 of
! these three terms. Breifly fa=-d/dratm, fb=-d/d(bi), fc=-d/d(bj)
! where matrix elements are < psi(r-bi) ! v(r-ratm) ! psi(r-bj)>
! ===========================================================================
! JOM : adapted to also calculate the gradient of the Hamiltonian
! G < mu | H | nu >, 3C-part
! gh_3c (ix,katom,imu,inu,ineigh,iatom)
! gh_3c is the gradient wrt to katom of the 3-C contribution to matrix element 
! < mu | H | nu >, mu in iatom, and nu in ineigh of iatom
! ===========================================================================
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
 
! Program Declaration
! ===========================================================================
        subroutine Dassemble_3c (nprocs, iordern, igauss) 
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
        integer ix
        integer jatom
        integer jbeta
        integer jssh
        integer mneigh
        integer my_proc
        integer natomsp
        integer jneigh
 

        real cost
        real x
        real y
 
        real, dimension (numorb_max, numorb_max) :: bcnax
        real, dimension (3, 3, 3) :: depsA
        real, dimension (3, 3, 3) :: depsB
        real, dimension (3, 3) :: eps
        real, dimension (3, numorb_max, numorb_max) :: f3naXa
        real, dimension (3, numorb_max, numorb_max) :: f3naXb
        real, dimension (3, numorb_max, numorb_max) :: f3naXc
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: rhat
        real, dimension (3) :: rna
        real, dimension (3) :: rnabc
        real, dimension (3) :: sighat

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Procedure
! ========================================================================
! Initialize the force contributions to zero.
        f3naa = 0.0d0
        f3nab = 0.0d0
        f3nac = 0.0d0
! JOM : Initialize gh_3c to zero.
        gh_3c = 0.0d0

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
!$omp parallel do private (rna, indna, ineigh, mneigh, iatom, ibeta, r1, in1) &
!$omp&  private (jatom, jbeta, r2, in2, r21, y, sighat, rnabc, x, rhat, cost) &
!$omp&  private (eps, depsA, depsB, isorp, interaction, bcnax, f3naXa, f3naXb)&
!$omp&  private (f3naXc, imu, inu, ix)
        do ialp = iatomstart, iatomstart - 1 + natomsp
         rna(:) = ratom(:,ialp)
         indna = imass(ialp)

 
! Loop over the neighbors of each ialp.
! Now look at all common neighbor pairs which were figured out in main.
! The common neighbors were figured out in common_neighbors.f90
         do ineigh = 1, neigh_comn(ialp)
          mneigh = neigh_comm(ineigh,ialp)
          write(*,*)'mneigh',mneigh
 
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
 
! dera3 = depsA = deps/dratm in the 3-center system
! der13 = depsB = deps/dr1 in the 3-center system
           call deps3center (r1, r2, r21, y, rna, rnabc, eps, depsA, depsB)
 
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
            call Dtrescentros (interaction, isorp, isorpmax, in1,       &
     &                         in2, indna, x, y, cost, eps, depsA,      &
     &                         depsB, rhat, sighat, bcnax, f3naXa,      &
     &                         f3naXb, f3naXc, nspecies)
           else
            call DtrescentrosG_VNA (in1, in2, indna, x, y, cost, eps,   &
     &                              depsA, depsB, rhat, sighat, bcnax,  &
     &                              f3naXa, f3naXb, f3naXc, rcutoff)
           end if

! The terms f3naXa, f3naXb, and f3naXc are already force-like.
!$omp critical (Dna3)
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             do ix = 1, 3
              f3naa(ix,ialp) = f3naa(ix,ialp)                           &
     &         + 2.0d0*rho(imu,inu,mneigh,iatom)*f3naXa(ix,imu,inu)*eq2
              f3nab(ix,iatom) = f3nab(ix,iatom)                         &
     &        + 2.0d0*rho(imu,inu,mneigh,iatom)*f3naXb(ix,imu,inu)*eq2
              f3nac(ix,jatom) = f3nac(ix,jatom)                         &
     &         + 2.0d0*rho(imu,inu,mneigh,iatom)*f3naXc(ix,imu,inu)*eq2
             end do
            end do
           end do
! JOM : gh_3c
! JOM : Notice the minus sign since
! the terms f3naXa, f3naXb, and f3naXc are force-like.
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             do ix = 1, 3
              gh_3c(ix,ialp,imu,inu,mneigh,iatom) =                     &
     &        gh_3c(ix,ialp,imu,inu,mneigh,iatom) -                     &
     &         f3naXa(ix,imu,inu)*eq2

              
              gh_3c(ix,ialp,inu,imu,jneigh,jatom) =                     &
     &        gh_3c(ix,ialp,imu,inu,mneigh,iatom) 

              gh_3c(ix,iatom,imu,inu,mneigh,iatom) =                     &
     &        gh_3c(ix,iatom,imu,inu,mneigh,iatom) -                     &
     &         f3naXb(ix,imu,inu)*eq2


              gh_3c(ix,iatom,inu,imu,jneigh,jatom) =                     &
     &        gh_3c(ix,iatom,imu,inu,mneigh,iatom) 

              gh_3c(ix,jatom,imu,inu,mneigh,iatom) =                     &
     &        gh_3c(ix,jatom,imu,inu,mneigh,iatom) -                     &
     &         f3naXc(ix,imu,inu)*eq2


              gh_3c(ix,jatom,inu,imu,jneigh,jatom) =                     &
     &        gh_3c(ix,jatom,imu,inu,mneigh,iatom) 

             end do
            end do
           end do
! JOM-end : gh_3c
!
!$omp end critical (Dna3) 
! ****************************************************************************
             write(*,*)'**********GH_3C**************'
           write(*,*)'ialp, iatom, jatom', ialp, iatom, jatom
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
      write(*,*)'gh_3c=', gh_3c(1,ialp,imu,inu,mneigh,iatom), imu, inu, ialp            
      write(*,*)'gh_3c=', gh_3c(1,iatom,imu,inu,mneigh,iatom), imu, inu, iatom            
      write(*,*)'gh_3c=', gh_3c(1,jatom,imu,inu,mneigh,iatom), imu, inu, jatom            
            end do
           end do

! End loop over ialp and its common neighbors.
          end if
         end do
        end do

! Format Statements
! ===========================================================================
 
        return
        end
