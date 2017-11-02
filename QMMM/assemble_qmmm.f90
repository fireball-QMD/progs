! Program Declaration
! ===========================================================================
        subroutine assemble_qmmm (nprocs, iordern) 
        use charges
        use configuration
        use constants_fireball
        use dimensions
        use interactions
        use neighbor_map
        use energy
        use qmmm_module, only : qmmm_struct, qmmm_nml
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
        integer in3
        integer ineigh
        integer issh
        integer jatom
        integer katom
        integer gatom
        integer mbeta
        integer my_proc
        integer natomsp
 
        real distance12
        real dij
        real dterm
        real sterm
        real out_charge
        real dq3
        real dq4
 
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (natoms) :: sub_ewaldqmmm
!        real, dimension (qmmm_struct%qm_mm_pairs) :: mm_charges
        real, dimension (qmmm_struct%qm_mm_pairs) :: charge_new

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE


! Procedure
! ===========================================================================
! Initialize interactions to zero.

        ewaldqmmm = 0.0d0

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

        sub_ewaldqmmm = 0.0d0

!        out_charge = 0
!        do katom = 1, qmmm_struct%qm_mm_pairs
!          out_charge = out_charge + qmmm_struct%qm_xcrd(4,katom)
!        end do
!        write(*,*) 'Total mm charge = '
!        out_charge = out_charge / qmmm_struct%qm_mm_pairs
!        write(*,*) 'Delta charge = ', out_charge
!        do katom = 1, qmmm_struct%qm_mm_pairs
!          mm_charges(katom) = qmmm_struct%qm_xcrd(4,katom) - out_charge
!        end do

!        do gatom = 1, qmmm_struct%qm_mm_pairs
!          if ( qmmm_struct%qm_xcrd_dist(gatom) > 0.8*qmmm_nml%qmcut ) then
!            qmmm_struct%qm_xcrd(4,gatom) = qmmm_struct%qm_xcrd(4,gatom)*(((qmmm_nml%qmcut - qmmm_struct%qm_xcrd_dist(gatom))/(0.2*qmmm_nml%qmcut)))
!            write(*,*) 'prueba',gatom,sqrt(qmmm_struct%qm_xcrd_dist(gatom)), qmmm_struct%qm_xcrd(4,gatom)!
!          end if
!        end do

! Calculate the charge on jatom

        do iatom = iatomstart, iatomstart - 1 + natomsp
          do katom = 1, qmmm_struct%qm_mm_pairs

            dij = sqrt ( (qmmm_struct%qm_xcrd(1,katom)-ratom(1,iatom))*(qmmm_struct%qm_xcrd(1,katom)-ratom(1,iatom)) + &
                         (qmmm_struct%qm_xcrd(2,katom)-ratom(2,iatom))*(qmmm_struct%qm_xcrd(2,katom)-ratom(2,iatom)) + & 
                         (qmmm_struct%qm_xcrd(3,katom)-ratom(3,iatom))*(qmmm_struct%qm_xcrd(3,katom)-ratom(3,iatom)) )

            sub_ewaldqmmm(iatom) = sub_ewaldqmmm(iatom) - (qmmm_struct%qm_xcrd(4,katom) / dij)
!            sub_ewaldqmmm(iatom) = sub_ewaldqmmm(iatom) - (mm_charges(katom) / dij)
          end do
        end do


!$imp parallel do private (r1, in1, ineigh, jatom, mbeta, r2, in2, distance12) &
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
!            write(*,*) 'iatom',iatom,'ineigh',ineigh,'inu',inu,'imu',imu
!            write(*,*) 'sterm',sterm,'dterm',dterm
            ewaldqmmm(imu,inu,ineigh,iatom) =  ewaldqmmm(imu,inu,ineigh,iatom)    &
     &       + (sterm - dterm)*sub_ewaldqmmm(iatom)*eq2                          &
     &       + (sterm + dterm)*sub_ewaldqmmm(jatom)*eq2
!           write(*,*) 'imu inu ineigh iatom ewaldqmmm', imu, inu, ineigh, iatom, ewaldqmmm(imu,inu,ineigh,iatom)
!           write(*,*) 'ewaldqmmm=',  ewaldqmmm(imu,inu,ineigh,iatom)
           end do
          end do
 
! End loop over neighbors
         end do
 
! End loop over atoms
        end do

        eqmmm = 0.0d0
        do iatom = iatomstart, iatomstart - 1 + natomsp
          in3 = imass(iatom)
          dq3 = 0.0d0
          do issh = 1, nssh(in3)
               dq3 = dq3  + Qneutral(issh,in3)
          end do
          eqmmm = eqmmm - dq3*sub_ewaldqmmm(iatom)*eq2
        end do

! Format Statements
! ===========================================================================
        return
        end
