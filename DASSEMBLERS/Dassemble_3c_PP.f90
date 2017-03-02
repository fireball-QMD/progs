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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
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
! ===========================================================================
 
! Program Declaration
! ===========================================================================
        subroutine Dassemble_3c_PP (nprocs, impi) 
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
        integer inu
        integer issh
        integer ix
        integer jatom
        integer jbeta
        integer jssh
        integer m31
        integer m32
        integer mneigh
        integer my_proc
        integer natomsp
        integer ncc
 
        integer, external :: mpairnay

        real, dimension (numorb_max, numorb_max) :: bcnlx
        real, dimension (numorb_max) :: cl
        real, dimension (3, numorb_max, numorb_max) :: f3nlXa
        real, dimension (3, numorb_max, numorb_max) :: f3nlXb
        real, dimension (3, numorb_max, numorb_max) :: f3nlXc
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r31
        real, dimension (3) :: r32
        real, dimension (3) :: rna

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Procedure
! ========================================================================
! Initialize the force contributions to zero.
        f3nla = 0.0d0
        f3nlb = 0.0d0
        f3nlc = 0.0d0

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

! Choose atom ialp in the central cell. This is the atom whose position
! we take the derivative, and is the atom who has the the neutral atom
! potential.
! Loop over the atoms in the central cell.
!$omp parallel do private (rna, indna, cl, ineigh, mneigh, iatom, ibeta, r1) &
!$omp&  private (in1, jatom, jbeta, r2, in2, r31, r32, m31, m32, inu, imu)   &
!$omp&  private (f3nlXa, f3nlXb, f3nlXc, bcnlx, ncc, ix)
        do ialp = iatomstart, iatomstart - 1 + natomsp
         rna(:) = ratom(:,ialp)
         indna = imass(ialp)

! For the non-local potential pieces call find the coefficients corresponding
! to indna.
         call cl_value (indna, cl)
 
! Loop over the neighbors of each ialp.
! Now look at all common neighbor pairs which were figured out in main.
! The common neighbors were figured out in common_neighbors.f90
         do ineigh = 1, neighPP_comn(ialp)
          mneigh = neighPP_comm(ineigh,ialp)
 
! The second atom (jatom) is the mneigh'th neighbor of iatom.
          if (mneigh .ne. 0) then
           iatom = neighPP_comj(1,ineigh,ialp)
           ibeta = neighPP_comb(1,ineigh,ialp)
           r1(:) = ratom(:,iatom) + xl(:,ibeta)
           in1 = imass(iatom)

           jatom = neighPP_comj(2,ineigh,ialp)
           jbeta = neighPP_comb(2,ineigh,ialp)
           r2(:) = ratom(:,jatom) + xl(:,jbeta)
           in2 = imass(jatom)
 
! ****************************************************************************
!
! PERFORM ACTUAL CALCULATIONS
! ASSEMBLE NON-LOCAL PSEUDOPOTENTIAL PIECE.
! ****************************************************************************
! Here 1=iatom, 2=jatom, 3=ialp:  <1 | V(3) | 2>.
           r31(:) = rna(:) - r1(:)
           r32(:) = rna(:) - r2(:)
 
! Find m value for iatom, ialp pair, and for jatom, ialp pair.
           m31 = mpairnay (iatom, ialp, r31)
           m32 = mpairnay (jatom, ialp, r32)
 
! A "general" matrix element is:
! <iatom|VNL(3)|jatom> = <iatom|V(3)><V(3)|jatom> =
! SUM_ncc cl(ncc)*sVNL(mu,ncc,iatom,m31)*sVNL(nu,ncc,jatom,m32)
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             f3nlXb(:,imu,inu) = 0.0d0
             f3nlXc(:,imu,inu) = 0.0d0
             bcnlx(imu,inu) = 0.0d0
             do ncc = 1, num_orbPP(indna)
              do ix = 1, 3
               f3nlXb(ix,imu,inu) = f3nlXb(ix,imu,inu)                       &
     &          - cl(ncc)*spVNL(ix,imu,ncc,m31,iatom)*sVNL(inu,ncc,m32,jatom)
               f3nlXc(ix,imu,inu) = f3nlXc(ix,imu,inu)                       &
     &         - cl(ncc)*sVNL(imu,ncc,m31,iatom)*spVNL(ix,inu,ncc,m32,jatom)
              end do
              bcnlx(imu,inu) = bcnlx(imu,inu)                                &
     &         + cl(ncc)*sVNL(imu,ncc,m31,iatom)*sVNL(inu,ncc,m32,jatom)
             end do
             f3nlXa(:,imu,inu) = - f3nlXb(:,imu,inu) - f3nlXc(:,imu,inu)
            end do
           end do

! Now write the forces to f3nla, f3nlb, and f3nlc. Multiply by rho.
!$omp critical (Dnl3)
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             do ix = 1, 3
              f3nla(ix,ialp) = f3nla(ix,ialp)                                &
     &         + rhoPP(imu,inu,mneigh,iatom)*f3nlXa(ix,imu,inu)
              f3nlb(ix,iatom) = f3nlb(ix,iatom)                              &
     &         + rhoPP(imu,inu,mneigh,iatom)*f3nlXb(ix,imu,inu)
              f3nlc(ix,jatom) = f3nlc(ix,jatom)                              &
     &         + rhoPP(imu,inu,mneigh,iatom)*f3nlXc(ix,imu,inu)
             end do
            end do
           end do
!$omp end critical (Dnl3)
! ****************************************************************************
! End loop over ialp and its common neighbors.
          end if
         end do
        end do

! Format Statements
! ===========================================================================
 
        return
        end subroutine Dassemble_3c_PP
