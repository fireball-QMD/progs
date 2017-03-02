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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
! assemble_lr.f90
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
        subroutine assemble_lr (nprocs, iordern) 
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
 
        real distance12
        real dq1
        real dterm
        real sterm
 
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (natoms) :: sub_ewald

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ===========================================================================
! Initialize interactions to zero.
        ewaldlr = 0.0d0

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

! Initialize here because it will be summed over procs
        sub_ewald = 0.0d0
 
! First a preliminary quantity. ewald(i,j) = sum(L) 1/(| bi-bj +/- L|
! We need to calculate SUMS of ewald sums in which the charge is included.
! These will also be used later in Pnanl3pP.f.  The details are on page 10 of
! the "Derivative of the long range matrix element". In the nutshell we are
! calculating:
!              Sum_(i,L) q(i)/|b(i)-b(alpha)+L|
 
!$omp parallel do private (jatom, in2, dq1, issh)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         do jatom = 1, natoms
          in2 = imass(jatom)
 
! Calculate the charge on jatom
          dq1 = 0.0d0
          do issh = 1, nssh(in2)
           dq1 = dq1 + (Qin(issh,jatom) - Qneutral(issh,in2))
          end do
          sub_ewald(iatom) = sub_ewald(iatom) + dq1*ewald(iatom,jatom)
         end do
        end do

! sum sub_ewald over procs
        if (iordern .eq. 1) call assemble_ordern_sub_ewald (natoms, sub_ewald)
 
! Now the meat of the calculation.  Construct ewaldlr(mu,nu,i,m) ===>
! the matrix elements of the long-range parts of the Hamiltonian.
! We make matrix elements for the Long Range Ewald according to our theory:
! ewaldlr(mu,nu,ineigh,iatom) =
! {s(mu,nu,ineigh,iatom)/2}*SUM(j_basis)(Qin(jatom) - Qneutral(jatom))
!                                        *(ewald(iatom,jatom)
!                                          + ewald(ineigh,jatom))*eq2
! eq2 makes it into the units of eV.
! Loop over the atoms in the central cell.

!$omp parallel do private (r1, in1, ineigh, jatom, mbeta, r2, in2, distance12) &
!$omp&  private (inu, imu, sterm, dterm)
        do iatom = iatomstart, iatomstart - 1 + natomsp
 
! Position of iatom.
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)
 
! Loop over the neighbors of the atom i.
         do ineigh = 1, neighn(iatom)
          jatom = neigh_j(ineigh,iatom)
          mbeta = neigh_b(ineigh,iatom)
 
! Position of jatom.
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)
          distance12 = sqrt((r2(1) - r1(1))**2 + (r2(2) - r1(2))**2          &
     &                                         + (r2(3) - r1(3))**2)
 
! "Charge" on each atom of the bondcharge. We split the charge S and
! dipole p to be S/2-p/d on atom 1 and S/2+p/d on atom 2. If atom 1 is
! equal to atom 2, then drop the p term.
          do inu = 1, num_orb(in2)
           do imu = 1, num_orb(in1)
            sterm = s_mat(imu,inu,ineigh,iatom)/2.0d0
 
            if (distance12 .gt. 1.0d-4) then
             dterm = dip(imu,inu,ineigh,iatom)/distance12
            else
             dterm = 0.0d0
            endif

            ewaldlr(imu,inu,ineigh,iatom) = ewaldlr(imu,inu,ineigh,iatom)    &
     &       + (sterm - dterm)*sub_ewald(iatom)*eq2                          &
     &       + (sterm + dterm)*sub_ewald(jatom)*eq2
           end do
          end do
 
! End loop over neighbors
         end do
 
! End loop over atoms
        end do
 
! Format Statements
! ===========================================================================
        return
        end
