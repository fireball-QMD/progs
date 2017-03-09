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

 
! Dassemble_olsxc_2c.f90
! Program Description
! ===========================================================================
!       This routine assembles the OFF-SITE OLSXC interactions and
! two-center forces. The matrix elements look like <Psi(1)|V(1+2)|Psi(2)>.
!
!  OFF-SITE means that (1) and (2) are at different atoms.
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
! J.Ortega & J.P.Lewis
! ===========================================================================
 
! Program Declaration
! ===========================================================================
        subroutine Dassemble_olsxc_2c (nprocs, iordern )
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
        integer, intent (in) :: nprocs
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer iatomstart
        integer ierror
        integer inu
        integer imu
        integer in1
        integer in2
        integer in3
        integer index1
        integer index2
        integer ineigh
        integer interaction
        integer isorp
        integer issh
        integer ix
        integer jatom
        integer jssh
        integer kforce
        integer l1
        integer l2
        integer mbeta
        integer matom
        integer my_proc
        integer n1
        integer n2
        integer natomsp
 

        real y
        real muxc
        real dmuxc
        real d2muxc
        real exc
        real dexc
        real d2exc
        real sx
        real rho_av
        real rhoin

        real, dimension (numorb_max, numorb_max) :: bcxcx
        real, dimension (3, numorb_max, numorb_max) :: bcxcpx
        real, dimension (3, numorb_max, numorb_max) :: mxcb
        real, dimension (3) :: rhoinp
        real, dimension (3) :: rhop_av
        real, dimension (3) :: spx
        real, dimension (3, 3) :: eps
        real, dimension (3, 3, 3) :: deps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: sighat

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Procedure
! ========================================================================
! Initialize the force contributions to zero.
! JOM info : the notation in Dassemble_olsxc_2c and Dassemble_olsxc_on
! is misleading. we should use "faxc" in Dassemble_olsxc_on (atm case)
! and "fotxc" in Dassemble_olsxc_2c (ontop cases)

! JOM I change the notation
!       faxc = 0.0d0
        fotxc = 0.0d0

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

!***************************************************************************
!                        T W O - C E N T E R   P A R T
! 
!***************************************************************************  
! We need to calculate a "Horsfield"-like contribution, its two-center 
! correction in the average density approx., and the snxc-2c
! contributions
!--------------------------------------------------------------------------
! Loop over the atoms in the central cell.
!$omp parallel do private (matom, r1, in1, ineigh, mbeta, jatom, r2, in2, r21) &
!$omp&  private (y, eps, deps, isorp, kforce, interaction, in3, bcxcx, bcxcpx) &
!$omp&  private (imu,inu,ix, n1, issh, l1, n2, jssh, l2, rho_av, rhop_av, exc) &
!$omp&  private (muxc, dexc, d2exc, dmuxc, d2muxc, index1, index2, sx, spx)    &
!$omp&  private (rhoin, rhoinp, mxcb, sighat)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         matom = neigh_self(iatom)
         r1(:) = ratom(:,iatom)
         in1 = imass(iatom)

! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)          ! <==== loop over i's neighbors
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          r2(:) = ratom(:,jatom) + xl(:,mbeta)
          in2 = imass(jatom)
 
          if (iatom .eq. jatom .and. mbeta .eq. 0) then
 
! Do nothing here - special case. Interaction already calculated in on-site
! stuff, i.e. assemble_olsxc_on.f90
! We calculate only off diagonal elements here.
 
          else

! SET-UP STUFF
!****************************************************************************
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
 
!***************************************************************************  
! "HORSFIELD"
!***************************************************************************  
! This is the ontop case for the exchange-correlation energy
           kforce = 1
           isorp = 0
           interaction = 6
           in3 = in2
           call doscentros (interaction, isorp, kforce, in1, in2, in3, y,   &
     &                      eps, deps, bcxcx, bcxcpx)
 
! Notice the explicit - sign which makes f force like.
           do inu = 1, num_orb(in3)
            do imu = 1, num_orb(in1)
             do ix = 1, 3
!             faxc(ix,ineigh,iatom) = faxc(ix,ineigh,iatom)            &
!    &         - rho(imu,inu,ineigh,iatom)*bcxcpx(ix,imu,inu)
              fotxc(ix,ineigh,iatom) = fotxc(ix,ineigh,iatom)            &
     &         - rho(imu,inu,ineigh,iatom)*bcxcpx(ix,imu,inu)
             end do
            end do
           end do 

!***************************************************************************  
! OLSXC-CORRECTION TO "HORSFIELD"
!***************************************************************************  
! Loop over shells.
           n1 = 0
           do issh = 1, nssh(in1)
            l1 = lssh(issh,in1)
            n1 = n1 + l1 + 1
            n2 = 0
            do jssh = 1, nssh(in2)
             l2 = lssh(jssh,in2)
             n2 = n2 + l2 + 1
             rho_av =  arhoij_off(issh,jssh,ineigh,iatom)
!JOM: I assume the following: check with Jimmy the farmer
             rhop_av(:) =  arhopij_off(:,issh,jssh,ineigh,iatom)
             call cepal (rho_av, exc, muxc, dexc, d2exc, dmuxc, d2muxc)

             do index1 = -l1, l1
              do index2 = -l2, l2
               imu = n1 + index1
               inu = n2 + index2
               sx = s_mat(imu,inu,ineigh,iatom)
               spx(:) = sp_mat(:,imu,inu,ineigh,iatom)
               rhoin =  rhoij_off(imu,inu,ineigh,iatom)
               rhoinp(:) =  rhopij_off(:,imu,inu,ineigh,iatom)

! mxcb is the derivative of XC matrix element w.r.t. atom 1
               mxcb(:,imu,inu) = spx(:)*(rho_av*dmuxc - muxc)                &
     &          + rhop_av(:)*d2muxc*(rho_av*sx - rhoin) - rhoinp(:)*dmuxc
              end do ! do index
             end do ! do index
! End loop over shells.
             n2 = n2 + l2
            end do
            n1 = n1 + l1
           end do
 
!***************************************************************************  
! SNXC-2C
!***************************************************************************  
! Loop over shells.
           n1 = 0
           do issh = 1, nssh(in1)
            l1 = lssh(issh,in1)
            n1 = n1 + l1 + 1
            n2 = 0
            do jssh = 1, nssh(in2)
             l2 = lssh(jssh,in2)
             n2 = n2 + l2 + 1
             rho_av =  arho_off(issh,jssh,ineigh,iatom)
!JOM: I assume the following: check with Jimmy the farmer
             rhop_av =  arhop_off(:,issh,jssh,ineigh,iatom)
             call cepal (rho_av, exc, muxc, dexc, d2exc, dmuxc, d2muxc)
             do index1 = -l1, l1
              do index2 = -l2, l2
               imu = n1 + index1
               inu = n2 + index2 
 
               sx = s_mat(imu,inu,ineigh,iatom)
               spx(:) = sp_mat(:,imu,inu,ineigh,iatom)
               rhoin =  rho_off(imu,inu,ineigh,iatom)
               rhoinp(:) =  rhop_off(:,imu,inu,ineigh,iatom)

! mxcb is the derivative of XC matrix element w.r.t. atom 1
               mxcb(:,imu,inu) = mxcb(:,imu,inu)                             &
     &          + spx(:)*(muxc - rho_av*dmuxc)                               &
     &          + rhop_av(:)*d2muxc*( rhoin - rho_av*sx ) + rhoinp(:)*dmuxc

!              faxc(:,ineigh,iatom) = faxc(:,ineigh,iatom)                   &
!    &          - rho(imu,inu,ineigh,iatom)*mxcb(:,imu,inu)
               fotxc(:,ineigh,iatom) = fotxc(:,ineigh,iatom)                   &
     &          - rho(imu,inu,ineigh,iatom)*mxcb(:,imu,inu)
              end do ! do index
             end do ! do index
! End loop over shells.
             n2 = n2 + l2
            end do
            n1 = n1 + l1
           end do
 
!***************************************************************************  
!
!****************************************************************************
! End loop over iatom and its neighbors - jatom.
          end if
         end do
        end do
 
!
! Format Statements
! ===========================================================================
 
        return
      end subroutine Dassemble_olsxc_2c
