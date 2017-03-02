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
 
! Dassemble_hxc_2c.f90
! Program Description
! ===========================================================================
!       This routine assembles all of the two-center and degenerate
! two-center interactions for the Horsfield exchange-correlation.
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
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_hxc_2c (nprocs, iordern, itheory, igauss)
        use charges
        use configuration
        use density
        use dimensions
        use forces
        use interactions
        use neighbor_map
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: igauss
        integer, intent (in) :: iordern
        integer, intent (in) :: itheory
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
        integer in3
        integer ineigh
        integer interaction
        integer inu
        integer isorp
        integer issh
        integer ix
        integer jatom
        integer kforce
        integer matom
        integer mbeta
        integer my_proc
        integer natomsp
 
        real cost
        real dq1
        real dq2
        real y
 
        real, dimension (numorb_max, numorb_max) :: bcxcx
        real, dimension (3, numorb_max, numorb_max) :: bcxcpx
        real, dimension (3, 3) :: eps
        real, dimension (3, 3, 3) :: deps
        real, dimension (4) :: dqfact          ! there are four different isorps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: sighat

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Procedure
! ===========================================================================
! Initialize forces to zero.
        faxc = 0.0d0
        fotxc = 0.0d0
        if (itheory .eq. 1) faxc_ca = 0.0d0
        if (itheory .eq. 1) fotxc_ca = 0.0d0
        if (igauss .eq. 1) then
        end if

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
 
! Loop over the atoms in the central cell.
!$omp parallel do private (in1, in2, in3, interaction, isorp, jatom, kforce) &
!$omp &           private (matom, mbeta, dq1, dq2, y, bcxcx, bcxcpx, eps)    &
!$omp &           private (deps, dqfact, r1, r2, r21, sighat) 
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neigh_self(iatom)
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)

! Find charge on iatom.
         dq1 = 0.0d0
         do issh = 1, nssh(in1)
          dq1 = dq1 + (Qin(issh,iatom) - Qneutral(issh,in1))
         end do

! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)     ! <==== loop 2 over iatom's neighbors
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)
 
! Find charge on jatom.
          dq2 = 0.0d0
          do issh = 1, nssh(in2)
           dq2 = dq2 + (Qin(issh,jatom) - Qneutral(issh,in2))
          end do

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
          else
           sighat(:) = r21(:)/y
          end if
 
          call epsilon (r2, sighat, eps)
          call deps2cent (r1, r2, eps, deps)
 
 
! ****************************************************************************
!
! ASSEMBLE EXCHANGE-CORRELATION FOR ATM AND ONTOP CASE - NA CASE
! ****************************************************************************
! Check to make sure we are not doing <1|V(1)|1>. This is done in the atm case.
          if (iatom .eq. jatom .and. mbeta .eq. 0) then
 
! Do nothing here - special case. This interaction comes from unicentro.f
 
          else
 
! This is the atm case for the exchange-correlation energy
           isorp = 0
           kforce = 1              ! calculate forces here
           interaction = 7
           in3 = in1
           call doscentros (interaction, isorp, kforce, in1, in2, in3,  &
     &                      y, eps, deps, bcxcx, bcxcpx)
 
! Note: these loops are both over in1 indeed!  This is because the two
! wavefunctions are located on the same site, but the potential is located
! at a different site.
! Notice the explicit - sign which makes F force like.
           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
             do ix = 1, 3
!$omp atomic
              faxc(ix,ineigh,iatom) = faxc(ix,ineigh,iatom)             &
     &         - rho(imu,inu,matom,iatom)*bcxcpx(ix,imu,inu)
             end do
            end do
           end do

! This is the ontop case for the exchange-correlation energy
           isorp = 0
           interaction = 6
           in3 = in2
           call doscentros (interaction, isorp, kforce, in1, in2, in3,  &
     &                      y, eps, deps, bcxcx, bcxcpx)
 
! Notice the explicit - sign which makes f force like.
           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
             do ix = 1, 3
!$omp atomic
              fotxc(ix,ineigh,iatom) = fotxc(ix,ineigh,iatom)           &
     &         - rho(imu,inu,ineigh,iatom)*bcxcpx(ix,imu,inu)
             end do
            end do
           end do

! Use the gaussian approximation for the three-center interactions.
           if (igauss .eq. 1) then

! Form the right ontop force. Use the derivatives, since the derivative is
! with respect to d/d(ratom) when the atom is ontop iatom.
            cost = 1.0d0
            in3 = in2
            call dosgaussians (in1, in2, in3, y, cost, eps, deps, bcxcx,&
     &                         bcxcpx, rcutoff)
            do inu = 1, num_orb(in3)
             do imu = 1, num_orb(in1)
              do ix = 1, 3
!$omp atomic
               fotxc(ix,ineigh,iatom) = fotxc(ix,ineigh,iatom)          &
     &          - rho(imu,inu,ineigh,iatom)                             &
     &            *nuxc_3c(imu,inu,ineigh,iatom)*bcxcpx(ix,imu,inu)
              end do
             end do
            end do
           end if
          end if


! ****************************************************************************
!
! ASSEMBLE EXCHANGE-CORRELATION FOR ATM AND ONTOP CASE - CA CASE
! ****************************************************************************
! Check to make sure we are not doing <1|V(1)|1>. This is done in the atm case.
          if (itheory .eq. 1) then
           if (iatom .eq. jatom .and. mbeta .eq. 0) then

! Do nothing here - special case. This interaction comes from unicentro.f

           else

! Initialize dqfact; dq1 and dq2 are weighted by dq(in1) and q(in2)
            dqfact(1) = -dq1/(2.0d0*dq(in1))
            dqfact(2) =  dq1/(2.0d0*dq(in1))
            dqfact(3) = -dq2/(2.0d0*dq(in2))
            dqfact(4) =  dq2/(2.0d0*dq(in2))

! This is the atm case for the exchange-correlation energy
! Note - the summation over isorp leads to a first order expansion for vxc_ca
! with respect to the charges dq1 and dq2.
            interaction = 7
            in3 = in1
            do isorp = 1, 4
             call doscentros (interaction, isorp, kforce, in1, in2, &
     &                        in3, y, eps, deps, bcxcx, bcxcpx)

! Note: these loops are both over in1 indeed!  This is because the two
! wavefunctions are located on the same site, but the potential is located
! at a different site.
! Notice the explicit - sign which makes f force like.
             do inu = 1, num_orb(in3)
              do imu = 1, num_orb(in1)
               do ix = 1, 3
!$omp atomic
               faxc_ca(ix,ineigh,iatom) = faxc_ca(ix,ineigh,iatom) -    &
     &         rho(imu,inu,matom,iatom)*bcxcpx(ix,imu,inu)*dqfact(isorp)
               end do
              end do
             end do
            end do

! This is the ontop case for the exchange-correlation energy
            interaction = 6
            in3 = in2
            do isorp = 1, 4
             call doscentros (interaction, isorp, kforce, in1, in2, in3,&
     &                        y, eps, deps, bcxcx, bcxcpx)

! Notice the explicit - sign which makes f force like.
             do inu = 1, num_orb(in3)
              do imu = 1, num_orb(in1)
               do ix = 1, 3
!$omp atomic
               fotxc_ca(ix,ineigh,iatom) = fotxc_ca(ix,ineigh,iatom) -  &
     &         rho(imu,inu,ineigh,iatom)*bcxcpx(ix,imu,inu)*dqfact(isorp)
               end do
              end do
             end do
            end do
           end if
          end if
 
! ****************************************************************************
! End loop over iatom and its neighbors - jatom.
         end do
        end do

! Format Statements
! ===========================================================================
        return
        end
