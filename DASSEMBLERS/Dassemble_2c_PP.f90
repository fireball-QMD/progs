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
 
! Dassemble_2c.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the two-center and degenerate
! two-center interactions of PP matrices.
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
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_2c_PP (nprocs, impi)
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
        integer iatom
        integer iatomstart
        integer ierror
        integer imu
        integer in1
        integer in2
        integer ineigh
        integer inu
        integer ix
        integer jatom
        integer jneigh
        integer kneigh
        integer matom
        integer mbeta
        integer mneigh_self
        integer my_proc
        integer natomsp
        integer ncc
 
        real, dimension (numorb_max) :: cl
        real, dimension (3, numorb_max, numorb_max) :: fnlb

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ===========================================================================
! Initialize forces to zero.
        fanl = 0.0d0
        fotnl = 0.0d0

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


! ****************************************************************************
!
! ASSEMBLE VNL ATM CASE  <phi_i|Psi_j><Psi_j|phi_i>
! **************************************************************************** 
! Loop over the atoms in the central cell.
!!$omp parallel do private (matom, in1, ineigh, mbeta, jatom, in2, cl)        &
!!$omp&  private (mneigh_self, inu, imu, fnlb, ncc, ix, jneigh, kneigh)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neighPP_self(iatom)
         in1 = imass(iatom)
 
! Loop over the neighbors of each iatom.
         do ineigh = 1, nPPn(iatom)   ! <==== loop 2 over iatom's neighbors
          mbeta = nPP_b(ineigh,iatom)
          jatom = nPP_j(ineigh,iatom)
          in2 = imass(jatom)
  
! ****************************************************************************
!
! CALL DPSEUDOP AND GET VNL DERIVATIVE FOR ATM and ONTOP CASES.
! ****************************************************************************
! DO THE ATM CASE:
!
!           * R2 (NL)                   * R3 (NL)
!          /                           / \
!         /                           /   \
!        /dd1,dd2     compare to:    /     \
!       /                           /dd1    \ dd2
!      /                           /         \
!   <mu|nu>                      <mu|        |nu>
!    R1,R1 (WF)                  R1 (WF)      R2 (WF)
!
! ****************************************************************************
          call cl_value(in2,cl)
          mneigh_self = nPP_self(iatom)
          do inu = 1, num_orb(in1)
           do imu = 1, num_orb(in1)
            fnlb(:,imu,inu) = 0.0d0
 
! The result is zero if we have <1|V(1)|1>, which means 
! ineigh = neigh_self(iatom)
!
! fnlb   = n1Xn2 array in which the derivatives d/d(r-atom1) of
!            pseudopotential matrix elements (crystal coordinates)
!            are stored.
 
            if (ineigh .ne. mneigh_self) then
             do ncc = 1, num_orbPP(in2)
              do ix = 1,3
               fnlb(ix,imu,inu) = fnlb(ix,imu,inu)                           &
     &          + cl(ncc)*(spVNL(ix,imu,ncc,ineigh,iatom)                    &
     &                    *sVNL(inu,ncc,ineigh,iatom)                        &
     &                    + sVNL(imu,ncc,ineigh,iatom)                       &
     &                      *spVNL(ix,inu,ncc,ineigh,iatom))
              end do ! do ix
             end do ! do ncc
            end if ! if mneigh_self
           end do ! do imu
          end do ! do inu
 
! Notice the minus sign (-). This makes it force like.
          do inu = 1, num_orb(in1)
           do imu = 1, num_orb(in1)
            do ix = 1,3
             fanl(ix,ineigh,iatom) =  fanl(ix,ineigh,iatom)                  &
      &                - rhoPP(imu,inu,matom,iatom)*fnlb(ix,imu,inu)
            end do ! do ix
           end do ! do imu
          end do ! do  inu

         enddo ! do ineigh

!
! ASSEMBLE VNL ONTOP LEFT CASE   <phi_i|Psi_i><Psi_i|phi_j>
! ****************************************************************************
! Loop over the neighbors of each iatom.
         do ineigh = 1, nPPxn(iatom)   ! <==== loop 2 over iatom's neighbors
          mbeta = nPPx_b(ineigh,iatom)
          jatom = nPPx_j(ineigh,iatom)
          in2 = imass(jatom)

! ****************************************************************************
! DO THE ONTOP CASE:
!
!           R2 (WF)
!           |nu>
!          /
!         /
!        /
!       /
!      /
!  <mu|*
!    R1,R1 (WF;NL)
!
! ****************************************************************************
          if (iatom .eq. jatom .and. mbeta .eq. 0) then
 
! Do nothing here - special case.
 
          else
           call cl_value(in1,cl)
           mneigh_self = nPP_self(iatom)
           jneigh = nPPx_point(ineigh,iatom)
!
! fnlb   = n1Xn2 array in which the derivatives d/d(r-atom1) of
!          : pseudopotential matrix elements (crystal coordinates)
!            are stored.
 
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             fnlb(:,imu,inu) = 0.0d0
             do ncc = 1, num_orbPP(in1)
 
! Note - spVNL are always the derivative of the orbital and NOT the potential.
! Repeat. sVNL(i,m) = <i|V(m)> and spVNL = d/dri <i|VNL(m)>.
! What we need in the next line is
! d/drm <i|VNL(m)> = -d/dr1<i|VNL(m)> = -spVNL(i,m).
              do ix = 1,3
               fnlb(ix,imu,inu) = fnlb(ix,imu,inu)                           &
     &          - cl(ncc)*sVNL(imu,ncc,mneigh_self,iatom)                    &
     &                  *spVNL(ix,inu,ncc,jneigh,jatom)
              end do ! do ix
             end do ! do ncc
            end do ! do imu
           end do ! do inu

           kneigh = nPPx_map(ineigh,iatom)
! Note the (-1) sign, to make it force like.
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             do ix = 1,3
              fotnl(ix,ineigh,iatom) = fotnl(ix,ineigh,iatom)                &
     &         - rhoPP(imu,inu,kneigh,iatom)*fnlb(ix,imu,inu)
             end do ! do ix
            end do ! do imu
           end do ! do inu
          end if ! if (iatom)

! This was ontop right.
! We do not do ontop left, because we would double count it.
! We effectively do it by doing direct and cross terms in assemble_F.f
 
! ****************************************************************************
! End loop over iatom and its neighbors - jatom.
         end do ! do ineigh
        end do ! do iatom


! Format Statements
! ===========================================================================
        return
        end subroutine Dassemble_2c_PP
