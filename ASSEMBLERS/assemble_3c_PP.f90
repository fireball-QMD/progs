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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
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
!
! Program Declaration
! ===========================================================================
        subroutine assemble_3c_PP (nprocs, iordern)
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
        integer m31
        integer m32
        integer mneigh
        integer my_proc
        integer natomsp
        integer ncc
 
        integer, external :: mpairnay

        real, dimension (numorb_max, numorb_max) :: bcnlx
        real, dimension (numorb_max) :: cl
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
!$omp parallel do private (rna, indna, cl, ineigh, mneigh, iatom, ibeta)     &
!$omp&  private (r1, in1, jatom, jbeta, r2, in2, r31, r32, m31, m32, inu)    &
!$omp&  private (imu, ncc, bcnlx)  
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
             bcnlx(imu,inu) = 0.0d0
             do ncc = 1, num_orbPP(indna)
              bcnlx(imu,inu) = bcnlx(imu,inu)                                &
     &         + cl(ncc)*sVNL(imu,ncc,m31,iatom)*sVNL(inu,ncc,m32,jatom)
             end do ! do ncc
            end do ! do imu
           end do ! do inu
 
! Add this piece for iatom, jatom, and katom into the total (bcnlx ===> vnl)
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
!$omp atomic
             vnl(imu,inu,mneigh,iatom) =                                   &
     &        vnl(imu,inu,mneigh,iatom) + bcnlx(imu,inu)

            end do ! do imu
           end do ! do inu


! ****************************************************************************
! End loop over ialp and its common neighbors.
          end if
         end do ! do ineigh
        end do ! do ialp

! Format Statements
! ===========================================================================

        return
        end subroutine assemble_3c_PP
