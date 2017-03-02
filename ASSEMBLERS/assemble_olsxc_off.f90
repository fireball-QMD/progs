! copyright info:
!
!                             @Copyright 2005
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio State University - Dave Drabold
! University of Regensburg - Juergen Fritsch

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
! assemble_olsxc_2c.f90
! Program Description
! ===========================================================================
!       This routine assembles the two & three-center exchange-correlation
! for the average density approximation. 
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
        subroutine assemble_olsxc_off (nprocs, my_proc, iordern, itheory)
        use charges
        use density
        use dimensions
        use interactions
        use neighbor_map
	use configuration
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: iordern
        integer, intent (in) :: itheory
        integer, intent (in) :: my_proc
        integer, intent (in) :: nprocs


! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iatomstart
        integer imu
        integer in1, in2, in3
        integer ineigh
        integer interaction
        integer inu
        integer isorp
        integer jatom
        integer kforce
        integer matom
        integer mbeta
        integer natomsp

        real  y 
        real dxn
        real, dimension (numorb_max, numorb_max) :: bcxcx
        real, dimension (numorb_max, numorb_max) :: denmx
        real, dimension (numorb_max, numorb_max) :: den1x
        real, dimension (numorb_max, numorb_max) :: rhomx
        real, dimension (3, numorb_max, numorb_max) :: rhompx
        real, dimension (3, 3) :: eps
        real, dimension (3, 3, 3) :: deps
        real, dimension (3) :: r1, r2, r21
        real, dimension (3) :: sighat
        real, dimension (numorb_max, numorb_max) :: sx

! Procedure
! ===========================================================================
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

! Loop over the atoms in the central cell.
!$omp parallel do private ()
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neigh_self(iatom)
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)

! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)       ! <==== loop over i's neighbors
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)

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
! CALL DOSCENTROS AND GET VXC FOR ATM CASE - AVERAGE DENSITY APPROXIMATION 
! ****************************************************************************
! If r1 .ne. r2, then this is a case where the potential is located
! at one of the sites of a wavefunction (ontop case).
          if (iatom .eq. jatom .and. mbeta .eq. 0) then
 
! Do nothing here - special case. Interaction already calculated in atm case.
 
          else

! This is the ontop case for the exchange-correlation energy
! Horsfield like term <i mu|Vxc(rho_i+j)| j nu>
           isorp = 0
           interaction = 6
           in3 = in2
           call doscentros (interaction, isorp, kforce, in1, in2, in3, y,    &
     &                      eps, deps, rhomx, rhompx)
           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
!$omp atomic
             vxc(imu,inu,ineigh,iatom) =                                     &
     &        vxc(imu,inu,ineigh,iatom) + rhomx(imu,inu)
            end do
           end do

!dani.JOM[
          if ( itheory.eq.1 ) then
! the vxc_ontopl case, transfer of charge (McWeda)
           interaction = 18
           in3 = in2
           do isorp = 1, nssh(in1)
            call doscentros (interaction, isorp, kforce, in1, in1, in3, y,   &
     &                       eps, deps, rhomx, rhompx)

            dxn = (Qin(isorp,iatom) - Qneutral(isorp,in1))
            do inu = 1, num_orb(in3)
             do imu = 1, num_orb(in1)
             vxc_ca(imu,inu,ineigh,iatom) =                                 &
     &         vxc_ca(imu,inu,ineigh,iatom) + rhomx(imu,inu)*dxn
             end do
            end do
           end do


! the vxc_ontopr case, transfer of charge (McWeda)
           interaction = 19
           in3 = in2
           do isorp = 1, nssh(in2)
            call doscentros (interaction, isorp, kforce, in1, in2, in3, y,   &
     &                       eps, deps, rhomx, rhompx)

            dxn = (Qin(isorp,jatom) - Qneutral(isorp,in2))
            do inu = 1, num_orb(in3)
             do imu = 1, num_orb(in1)
              vxc_ca(imu,inu,ineigh,iatom) =                                 &
     &         vxc_ca(imu,inu,ineigh,iatom) + rhomx(imu,inu)*dxn
             end do
            end do
           end do
          endif
!dani.JOM]
!

! Restore density and overlap matrices
!  den1x    .... <i|n_i|j> + <i|n_j|j> (on-top interactions)
!  denmx    .... <i|n|j> = <i|n_i|j> + <i|n_j|j> + S_{i.ne.j.ne.k}<i|n_k|j>
!  sx       .... <i|j>
           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
!$omp atomic
             denmx(imu,inu) = rho_off(imu,inu,ineigh,iatom)
             den1x(imu,inu) = rhoij_off(imu,inu,ineigh,iatom)
             sx(imu,inu) = s_mat(imu,inu,ineigh,iatom)
            end do
           end do

! Calculate <i| V_xc(n) |j> and <i|V_xc(n_i+n_j)|j>              
           call build_olsxc_off (in1, in2, den1x, denmx, sx, ineigh, iatom,  &
     &                           bcxcx)

! now complete 'non-diagonal' terms <i|V(n)|j>
           if (itheory .eq. 0) then
            do inu = 1, num_orb(in2)
             do imu = 1, num_orb(in1)
              vxc(imu,inu,ineigh,iatom) =                                    &
     &         vxc(imu,inu,ineigh,iatom) + bcxcx(imu,inu)
             end do
            end do
           else
            do inu = 1, num_orb(in2)
             do imu = 1, num_orb(in1)
              vxc_ca(imu,inu,ineigh,iatom) =                                 &
     &         vxc_ca(imu,inu,ineigh,iatom) + bcxcx(imu,inu)
             end do
            end do
           end if

! End if for r1 .ne. r2 case
          end if
! ****************************************************************************
! End loop over iatom and its neighbors - jatom.
         end do
        end do

! Format Statements
! ===========================================================================

        return
        end subroutine assemble_olsxc_off
