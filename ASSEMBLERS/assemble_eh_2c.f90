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
 
! assemble_eh_2c.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the two-center extended hubbard
! interactions.
!
! ===========================================================================
! Original code written by:
! John Tomfohr and Otto F. Sankey
! Department of Physics & Astronomy
! Campus Box 871504
! Arizona State University
! Tempe, AZ 85287-1504
! (480) 965-0667 (office)      email: John.Tomfohr@asu.edu
! (480) 965-4334 (office)      email: Otto.Sankey@asu.edu
!
! Code rewritten by:
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
        subroutine assemble_eh_2c (nprocs, iordern, ioff2c)
        use charges
        use configuration
        use constants_fireball
        use integrals
        use interactions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iordern
        integer, intent (in) :: nprocs

        integer, intent (in), dimension (1:24) :: ioff2c
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iatomstart
        integer ideriv
        integer ierror
        integer imu
        integer in1
        integer in2
        integer index
        integer index_coulomb
        integer ineigh
        integer interaction
        integer inu
        integer issh
        integer jatom
        integer jssh
        integer kforce
        integer matom
        integer mbeta
        integer my_proc
        integer n1
        integer n2
        integer natomsp
 
        real distance
        real dq1
 
        real, dimension (nsh_max, nsh_max) :: coulomb
        real, dimension (nsh_max, nsh_max) :: coulombD
        real, dimension (ME2c_max) :: dslist
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (ME2c_max) :: slist
        real, dimension (natoms) :: sub_ewald
        real, dimension (nsh_max, nsh_max) :: xcnu
        real, dimension (nsh_max, nsh_max) :: xcnuD

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Procedure
! ===========================================================================
! A note about ioff2c: Using this for interaction = 22 sets all 1/r 
! interations equal to zero. However, it does not set the correct short 
! range coulomb interations (e.g. coulomb.14.14.dat, etc.) equal to zero.
! These can be set to zero by setting ioff2c = 0 for interaction = 24.
! This sets the U hubbard integrals (integral phi**2 1/r-r' phi**2)
! equal to zero in any hubbard calculation. You might say why not just set
! the coulomb integrals equal to zero when we read them in. I found out
! the hard way, that this does not work because assemble_usr uses the
! coulomb integrals even for a neutral atom.

! Initialize interactions to zero.
        vca = 0.0d0
        Vcoulomb = 0.0d0
        Vxcnu = 0.0d0
        vxc_ca = 0.0d0
        ewaldlr = 0.0d0
        ewaldsr = 0.0d0
        Vewaldsr = 0.0d0

! Determine which atoms are assigned to this processor.
        if (iordern .eq. 1) then
         call MPI_COMM_RANK (MPI_BTN_WORLD, my_proc, ierror)
         natomsp = natoms/nprocs
         if (my_proc .lt. mod(natoms,nprocs)) then
          natomsp = natomsp + 1
          iatomstart = natomsp*my_proc + 1
         else
          iatomstart = (natomsp + 1)*mod(natoms,nprocs)                     &  
                      + natomsp*(my_proc - mod(natoms,nprocs)) + 1
         end if
        else
         iatomstart = 1
         natomsp = natoms
        end if
        sub_ewald = 0.0d0

! First a preliminary quantity. ewald(i,j) = sum(L) 1/(| bi-bj +/- L|
! We need to calculate SUMS of ewald sums in which the charge is included.
! These will also be used later in Pnanl3pP.f.  The details are on page 10 of
! the "Derivative of the long range matrix element". In the nutshell we are
! calculating:
!       Sum_(i,L) q(i)/|b(i)-b(alpha)+L|

!$omp parallel do private (in2, dq1)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         do jatom = 1, natoms
          in2 = imass(jatom)

! Calculate the charge on jatom
          dq1 = 0.0d0
          do issh = 1, nssh(in2)
           dq1 = dq1 + (Qin(issh,jatom) - Qneutral(issh,in2))
          end do
!$omp atomic
          sub_ewald(iatom) = sub_ewald(iatom) + dq1*ewald(iatom,jatom)
         end do
        end do

! sum sub_ewald over procs
        if (iordern .eq. 1) call assemble_ordern_sub_ewald (natoms, sub_ewald)
 
! ****************************************************************************
! Loop over the atoms in the central cell.
!$omp parallel do private (in1, in2, ideriv, index_coulomb, interaction)     &
!$omp&            private (jatom, kforce, matom, mbeta, n1, n2, distance)    &
!$omp&            private (coulomb, coulombD, dslist, r1, r2, r21, slist)    &
!$omp&            private (xcnu, xcnuD) 
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neigh_self(iatom)
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)
 
! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)         ! <==== loop 2 over i's neighbors
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)
 
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
          r21(:) = r2(:) - r1(:)
          distance = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
 
! ****************************************************************************
!
! GET EXTENDED-HUBBARD COULOMB INTERACTIONS - STORE IN VCOULOMB.
! ****************************************************************************
! Now find the coulomb integrals need to evaluate the extended-hubbard
! interaction.
          index_coulomb = nssh(in1)*nssh(in2)
          
! Skip this if ioff2c = 0 for the coulomb (extended Hubbard) interactions. 
          if (ioff2c(24) .eq. 1) then

! Loop over all the non-zero integrals for this interaction:
           interaction = 12
           ideriv = 0
           kforce = 0
           
           do index = 1, index_coulomb
            call interpolate_1d (interaction, ideriv, in1, in2, index,       &
     &                           kforce, distance, slist(index), dslist(index))
           end do
 
! We have the data, it is stored in the following way: v(1,1), v(1,2),
! v(1,3)...v(1,n2max), v(2,1), v(2,2), ..., but it is a 1D array,
! ordered according to the index_coulomb. Restore this to true 2x2 format:
           n1 = nssh(in1)
           n2 = nssh(in2)
           call recoverC (n1, n2, slist, dslist, coulomb, coulombD)
 
           do issh = 1, nssh(in1)
            do jssh = 1, nssh(in2)
!$omp atomic
             Vcoulomb(issh,iatom) = Vcoulomb(issh,iatom)                     &
     &        + coulomb(issh,jssh)*(Qin(jssh,jatom) - Qneutral(jssh,in2))
            end do
           end do
          end if
 
 
! ****************************************************************************
!
! GET EXTENDED-HUBBARD EXCHANGE-CORRELATION INTERACTIONS - STORE IN VXCNU
! ****************************************************************************
! Now find the exchange-correlation integrals need to evaluate the
! extended-hubbard interaction.
! Loop over all the non-zero integrals for this interaction:
          interaction = 14
          ideriv = 0
          kforce = 0

          do index = 1, index_coulomb
           call interpolate_1d (interaction, ideriv, in1, in2, index,        &
     &                          kforce, distance, slist(index), dslist(index))
          end do
 
! We have the data, it is stored in the following way: v(1,1), v(1,2),
! v(1,3)...v(1,n2max), v(2,1), v(2,2), ..., but it is a 1D array,
! ordered according to the index_coulomb. Restore this to true 2x2 format:
          call recoverC (n1, n2, slist, dslist, xcnu, xcnuD)
 
          if (matom .ne. ineigh) then
!There are two different atoms.
           do issh = 1, nssh(in1)
            do jssh = 1, nssh(in2)
!$omp atomic
             Vxcnu(issh,iatom) = Vxcnu(issh,iatom)                           &
     &        + xcnu(issh,jssh)*(Qin(jssh,jatom) - Qneutral(jssh,in2))
            end do
           end do
          else
           do issh = 1, nssh(in1)
            do jssh = 1, nssh(in2)
!$omp atomic

! ============================================================================
! NOTE for SNXC and OLSXC in Extended Hubbard theory. (OCT. 15. 04. MHL)
! xcnu1c is the matrix elements calculated for Horsfield, but we use the same
! matrix elements for SNXC and OLSXC for the time being.
! It needs to be fixed.
! ============================================================================

             Vxcnu(issh,iatom) = Vxcnu(issh,iatom)                           &
     &        + xcnu1c(issh,jssh,in1)*(Qin(jssh,jatom) - Qneutral(jssh,in2))
            end do
           end do
          end if
 
! ****************************************************************************
!
! GET EXTENDED-HUBBARD SHORT-RANGE EWALD INTERACTIONS - STORE IN EWALDSR
! ****************************************************************************
          if (matom .ne. ineigh .and. ioff2c(22) .eq. 1) then
           do issh = 1, nssh(in1)
            do jssh = 1, nssh(in2)
!$omp atomic
             Vewaldsr(issh,iatom) = Vewaldsr(issh,iatom)                     &
     &        + (Qin(jssh,jatom) - Qneutral(jssh,in2))/distance
            end do
           end do
          end if
 
! ****************************************************************************
! End loop over neighbors of iatom - jatom.
         end do

! End loop iatom.
        end do
 
! ****************************************************************************
!
! STORE VCOULOMB, VXCNU, AND VEWALDSR IN VCA AND VXC_CA
! ****************************************************************************
! Loop over the atoms in the central cell.
!$omp parallel do private (in1, jatom, in2, jssh, issh)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         in1 = imass(iatom)
 
! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)         ! <==== loop 2 over i's neighbors
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom)
 
          do inu = 1, num_orb(in2)
           jssh = mu2shell(inu,in2)
           do imu = 1, num_orb(in1)
            issh = mu2shell(imu,in1)
!$omp atomic
            vca(imu,inu,ineigh,iatom) = vca(imu,inu,ineigh,iatom)            &
     &       + (eq2/2.0d0)*s_mat(imu,inu,ineigh,iatom)                       &
     &         *(Vcoulomb(issh,iatom) + Vcoulomb(jssh,jatom))
!$omp atomic
            vxc_ca(imu,inu,ineigh,iatom) = vxc_ca(imu,inu,ineigh,iatom)      &
     &       + s_mat(imu,inu,ineigh,iatom)                                   &
     &         *(Vxcnu(issh,iatom) + Vxcnu(jssh,jatom))/2.0d0
!$omp atomic
!This will be zero if ioff2c(22) = 0
            ewaldsr(imu,inu,ineigh,iatom) = ewaldsr(imu,inu,ineigh,iatom)    &
     &       + (eq2/2.0d0)*s_mat(imu,inu,ineigh,iatom)                       &
     &         *(Vewaldsr(issh,iatom) + Vewaldsr(jssh,jatom))
!$omp atomic
            ewaldlr(imu,inu,ineigh,iatom) = ewaldlr(imu,inu,ineigh,iatom)    &
     &       + (eq2/2.0d0)*s_mat(imu,inu,ineigh,iatom)                       &
     &         *(sub_ewald(iatom) + sub_ewald(jatom))
           end do
          end do
 
! End loop over iatom and its neighbors - jatom.
         end do
        end do

! ****************************************************************************
!
! NOW SUBTRACT OFF EWALD TERMS FROM VCOULOMB
! ****************************************************************************
! The new Vcoulomb will be passed back and used in assemble_eh_usr.f90
! Later, take Vcoulomb*Qout and subtract that from the uiiuee energy.
!$omp parallel do private(in1)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         in1 = imass(iatom)
         do issh = 1, nssh(in1)
!$omp atomic
          Vcoulomb(issh,iatom) = Vcoulomb(issh,iatom) - Vewaldsr(issh,iatom) &
     &                           + sub_ewald(iatom) 
         end do
        end do

        if (iordern .eq. 1) call assemble_eh_2c_ordern_final (natoms)
 
! ****************************************************************************
! Format Statements
! ===========================================================================
 
        return
        end
