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
 
! Dassemble_eh_2c.f90
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
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_eh_2c (nprocs, my_proc, iordern)
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
        integer, intent (in) :: my_proc
        integer, intent (in) :: nprocs 
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iatomstart
        integer ideriv
        integer iforce
        integer imu
        integer in1
        integer in2
        integer index
        integer index_coulomb
        integer ineigh
        integer interaction
        integer inu
        integer issh
        integer ix
        integer jatom
        integer jssh
        integer matom
        integer mbeta
        integer n1
        integer n2
        integer natomsp
 

        real distance
        real force
        real dq2
 
        real, dimension (nsh_max, nsh_max) :: coulomb
        real, dimension (nsh_max, nsh_max) :: coulombD
        real, dimension (ME2c_max) :: dslist
        real, dimension (3) :: eta 
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (ME2c_max) :: slist
        real, dimension (natoms) :: sub_ewald
        real, dimension (3, natoms) :: sub_dewald
        real, dimension (nsh_max, nsh_max) :: xcnu
        real, dimension (nsh_max, nsh_max) :: xcnuD
        real, dimension (nsh_max) :: dqi
        real, dimension (nsh_max) :: dqj
 
! Procedure
! ===========================================================================
! Initialize forces to zero.
        fcoulomb = 0.0d0
        flrew = 0.0d0
        fxcnu = 0.0d0
 
! Determine which atoms are assigned to this processor.
        if (iordern .eq. 1) then
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
 
! Initialize some arrays to zero.
        sub_ewald = 0.0d0
        sub_dewald = 0.0d0
 
! First a preliminary quantity. ewald(i,j) = sum(L) 1/(| bi-bj +/- L|
! We need to calculate SUMS of ewald sums in which the charge is included.
! These will also be used later in Pnanl3pP.f.  The details are on page 10 of
! the "Derivative of the long range matrix element". In the nutshell we are
! calculating:
!              Sum_(i,L) q(i)/|b(i)-b(alpha)+L|

!$omp parallel do private (in2, dq2)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         do jatom = 1, natoms
          in2 = imass(jatom)

! Calculate the charge on jatom
          dq2 = 0.0d0
          do issh = 1, nssh(in2)
           dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
          end do
          sub_ewald(iatom) = sub_ewald(iatom) + dq2*ewald(iatom,jatom)
          sub_dewald(:,iatom) = sub_dewald(:,iatom) + dq2*dewald(:,iatom,jatom)
         end do
        end do

! Sum sub_ewald, sub_dewald over procs
         if (iordern .eq. 1)                                                 &
     &    call assemble_ordern_sub_dewald (natoms, sub_ewald, sub_dewald)

! ****************************************************************************
! Loop over the atoms in the central cell.

!$omp parallel do private (matom, r1, in1, dq1, mbeta, jatom, r2, in2, dq2)  &
!$omp&            private (r21, distance, eta, index_coulomb, interaction)   &
!$omp&            private (ideriv, iforce, slist, dslist, n1, n2, coulomb)   &
!$omp&            private (coulombD, issh, jssh, xcnu, xcnuD)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neigh_self(iatom)
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)
 
! Calculate the charge on iatom in shell issh
         dqi = 0.0d0
         do issh = 1, nssh(in1)
          dqi(issh) = Qin(issh,iatom) - Qneutral(issh,in1)
         end do

! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)         ! <==== loop 2 over i's neighbors
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)
 
! Find charge on jatom in shell issh
          dqj = 0.0d0
          do issh = 1, nssh(in2)
           dqj(issh) = Qin(issh,jatom) - Qneutral(issh,in2)
          end do

! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
          r21(:) = r2(:) - r1(:)
          distance = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))
          if (matom .ne. ineigh) then
           eta(:) = (r2(:) - r1(:))/distance
          else
           eta(:) = 0.0d0
          end if


! ****************************************************************************
!
! GET EXTENDED-HUBBARD COULOMB FORCES - STORE IN FCOULOMB.
! ****************************************************************************
! Now find the coulomb integrals need to evaluate the extended-hubbard  
! interaction.
! Loop over all the non-zero integrals for this interaction:
          index_coulomb = nssh(in1)*nssh(in2)
          interaction = 12
          ideriv = 0
          iforce = 1
          do index = 1, index_coulomb
           call interpolate_1d (interaction, ideriv, in1, in2, index,  &
     &                      iforce, distance, slist(index), dslist(index))
          end do

! We have the data, it is stored in the following way: v(1,1), v(1,2),
! v(1,3)...v(1,n2max), v(2,1), v(2,2), ..., but it is a 1D array,
! ordered according to the index_coulomb. Restore this to true 2x2 format:
          n1 = nssh(in1)
          n2 = nssh(in2)
          call recoverC (n1, n2, slist, dslist, coulomb, coulombD)

          do issh = 1, nssh(in1)
           do jssh = 1, nssh(in2)
            do ix = 1, 3
!$omp atomic
! The exact integral for neighbors. This is the 1/2 J * Q *Q term.
! d=|bj+L - bi|, U(d). 
! f(k)=-dU/dbk =-dU/dd*(delta(k,j)*d|bj+L-bi|/dbj+delta(k,i)*d|bj+L-bi|/dbk)
!       = -dU/dd *(delta(k,j) eta * 1    + delta(k,i) eta * -1      )
! Note: If i=j, then f gives exactly zero for f(k), AS IT SHOULD.(U(L) only.)

	         force=(eq2/2.0d0)*eta(ix)*(-coulombD(issh,jssh))*dqi(issh)*dqj(jssh)
! k=jatom
	         fcoulomb(ix,jatom)=fcoulomb(ix,jatom)+force*(+1.0d0)
! k=iatom
	         fcoulomb(ix,iatom)=fcoulomb(ix,iatom)+force*(-1.0d0)
            end do
           end do
          end do

          do inu = 1, num_orb(in2)
           jssh = mu2shell(inu,in2)
           do imu = 1, num_orb(in1)
            issh = mu2shell(imu,in1)
            do ix = 1, 3
!$omp atomic
! The exact integral for neighbors. 
! This is the rho (Valpha+Vbeta)/2 * S' term.
             force= -rho(imu,inu,ineigh,iatom)                     &
     &         *sp_mat(ix,imu,inu,ineigh,iatom)*(eq2/2.0d0)   &
     &         *(Vcoulomb(issh,iatom) + Vcoulomb(jssh,jatom))
 
	         fcoulomb(ix,iatom) = fcoulomb(ix,iatom) + force
 	         fcoulomb(ix,jatom) = fcoulomb(ix,jatom) - force
            end do
           end do
          end do


! ****************************************************************************
!
! GET EXTENDED-HUBBARD EXCHANGE-CORRELATION FORCES - STORE IN FXCNU 
! ****************************************************************************
! Now find the exchange-correlation integrals need to evaluate the 
! extended-hubbard interaction.
! Loop over all the non-zero integrals for this interaction:
          interaction = 14
          ideriv = 0
          do index = 1, index_coulomb
           call interpolate_1d (interaction, ideriv, in1, in2, index,  &
     &                      iforce, distance, slist(index), dslist(index))
          end do

! We have the data, it is stored in the following way: v(1,1), v(1,2),
! v(1,3)...v(1,n2max), v(2,1), v(2,2), ..., but it is a 1D array,
! ordered according to the index_coulomb. Restore this to true 2x2 format:
          call recoverC (n1, n2, slist, dslist, xcnu, xcnuD)
          if (matom .eq. ineigh) then
           
          else
           do issh = 1, nssh(in1)
            do jssh = 1, nssh(in2)
             do ix = 1, 3
!$omp atomic
	          force = 0.5d0*eta(ix)*(-1.0d0*xcnuD(issh,jssh))         &
     &               *(Qin(issh,iatom) - Qneutral(issh,in1))  &
     &               *(Qin(jssh,jatom) - Qneutral(jssh,in2))
! k=jatom
              fxcnu(ix,jatom)=fxcnu(ix,jatom) + force
! k=iatom
              fxcnu(ix,iatom)=fxcnu(ix,iatom) - force
             end do
            end do
           end do

           do inu = 1, num_orb(in2)
            jssh = mu2shell(inu,in2)
            do imu = 1, num_orb(in1)
             issh = mu2shell(imu,in1)
             do ix = 1, 3
!$omp atomic
             force =  - rho(imu,inu,ineigh,iatom)              &
     &         *sp_mat(ix,imu,inu,ineigh,iatom)                &
     &         *(Vxcnu(issh,iatom) + Vxcnu(jssh,jatom))/2.0d0

              fxcnu(ix,iatom)=fxcnu(ix,iatom) + force
              fxcnu(ix,jatom)=fxcnu(ix,jatom) - force
             end do ! do ix
            end do ! do imu
           end do ! do inu
          endif ! if (matom .eq. ineigh)
          
! ****************************************************************************
!
! GET EXTENDED-HUBBARD SHORT_RANGE EWALD FORCES - STORE IN FCOULOMB.
! ****************************************************************************
          if (matom .ne. ineigh) then
           do issh = 1, nssh(in1)
            do jssh = 1, nssh(in2)
             do ix = 1, 3
!$omp atomic
! The ewald-like term for near neighbors. This is the 1/2 J Q Q term.
! Note: If i=j, then f gives exactly zero for f(k), AS IT SHOULD.(U(L) only.)
	          force = -(eq2/2.0d0)*eta(ix)*(-1.0d0)*((-1.0d0)/(distance**2))        &
     &         *(Qin(issh,iatom) - Qneutral(issh,in1))  &
     &         *(Qin(jssh,jatom) - Qneutral(jssh,in2))

! k=jatom
              fcoulomb(ix,jatom) = fcoulomb(ix,jatom) + force

! k=iatom
              fcoulomb(ix,iatom) = fcoulomb(ix,iatom) - force
             end do
            end do
           end do

           do inu = 1, num_orb(in2)
            jssh = mu2shell(inu,in2)
            do imu = 1, num_orb(in1)
             issh = mu2shell(imu,in1)
             do ix = 1, 3
!$omp atomic
! The ewald-like term for near neighbors. 
! This is the rho (Valpha+Vbeta)/2 S' term
! Note. Vewaldsr=0 if ioff_ewald.ne.1.

! Special Note. ewaldsr is SUBTRACTED from H. This adds another minus sign.
! Thiis minus is at the end of force =.....*(-1.)
! From buildh:  h_mat = h_mat + vca + vxc_ca + ewaldlr - ewaldsr

              force= -rho(imu,inu,ineigh,iatom)                     &
     &         *sp_mat(ix,imu,inu,ineigh,iatom)*(eq2/2.0d0)   &
     &         *(Vewaldsr(issh,iatom) + Vewaldsr(jssh,jatom))*(-1.0d0)

              fcoulomb(ix,iatom) = fcoulomb(ix,iatom) + force
              fcoulomb(ix,jatom) = fcoulomb(ix,jatom) - force
             end do
            end do
           end do
          end if


! ****************************************************************************
!
! GET EXTENDED-HUBBARD LONG_RANGE EWALD FORCES - STORE IN FCOULOMB.
! ****************************************************************************
          do inu = 1, num_orb(in2)
           do imu = 1, num_orb(in1)
            do ix = 1, 3
!$omp atomic
             flrew(ix,iatom) = flrew(ix,iatom)                               &
     &        - rho(imu,inu,ineigh,iatom)*sp_mat(ix,imu,inu,ineigh,iatom)    &
     &         *(eq2/2.0d0)*(sub_ewald(iatom) + sub_ewald(jatom))
            end do
           end do
          end do

! ****************************************************************************
! End loop over iatom and its neighbors - jatom.
         end do
         do ix = 1, 3
!$omp atomic
! Ewald-like term for all neighbors. The 1/2dr**(-1)/dd Q Q term.
          flrew(ix,iatom) = flrew(ix,iatom) + fewald(ix,iatom)*(eq2/2.0d0)
         end do
        end do
 
        if (iordern .eq. 1) call Dassemble_eh_2c_ordern_final (natoms)

! Format Statements
! ===========================================================================
 
        return
        end
