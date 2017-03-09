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

 
! Dassemble_olsxc_on_mdet.f90
! Program Description
! ===========================================================================
!       This routine assembles the forces for the McWEDA on-site two-center 
! exchange correlation interactions. 
!
! JOM : adapted to calculate also gh_atm ( < Grad mu | H | nu > ) for
! the non-adiabatic calculation
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
!
! Program Declaration
! ===========================================================================
        subroutine Dassemble_olsxc_on_mdet (nprocs, iordern)
        use charges
        use configuration
        use constants_fireball
        use density
        use interactions
        use forces
        use neighbor_map
        use nonadiabatic
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
        integer imu
        integer in1
        integer in2
        integer in3
        integer index
        integer ineigh
        integer interaction
        integer inu
        integer isorp
        integer issh
        integer ix
        integer jatom
        integer kforce
        integer jndex
        integer jssh
        integer l1 
        integer l2
        integer matom
        integer mbeta
        integer my_proc
        integer n1
        integer n2
        integer natomsp

        real exc
        real dexc
        real d2exc
        real d2muxc
        real dmuxc
        real muxc
        real rhoxc
        real rhoxc_av
        real y
        real q_mu

        real, dimension (3, numorb_max, numorb_max) :: bcxcpx
        real, dimension (numorb_max, numorb_max) :: delta
        real, dimension (3) :: drhoxc
        real, dimension (3) :: drhoxc_av

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE

! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Initialize forces to zero. 
! JOM info : the notation in Dassemble_olsxc_2c and Dassemble_olsxc_on
! is misleading. we should use "faxc" in Dassemble_olsxc_on (atm
! case)
! and "fotxc" in Dassemble_olsxc_2c (ontop cases)

! JOM I change the notation
!       fotxc = 0.0d0
        faxc = 0.0d0
        dxcdcc  = 0.0d0

! Comment about 'delta':
! Don't use 'delk' variable, dimension is different, 'delta' goes over orbitals
        delta = 0.0d0
        do imu = 1,numorb_max
         delta(imu,imu) = 1.0d0
        enddo

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


! ****************************************************************************
!
! ASSEMBLE SN EXCHANGE CORRELATION FORCE FOR ON-SITE PART  <i|rho|i> 
! ****************************************************************************
! We calculate onsite (two-center part) coming from assemble_oslxc_on ()
! terms:
!  f(r1) = - <mu i|nu i>*arho*Vxc''(arho)*d(arho)/dr1   
!          + <mu i|rho| nu i>*Vxc''(arho)*d(arho)/dr1  
!          + Vxc'(arho)*d(<mu i|rho| nu i>)/dr1
!
! Variables:
! <mu i|rho| nu i>   ... rho_on()
! arho               ... arho_on()
! d(arho)/dr1        ... arhop_on()

! Loop over the atoms in the central cell.
!$omp parallel do private (in1, matom, mbeta, jatom, in2, bcxcpx, n1, l1)     &
!$omp&    private (issh, rhoxc_av, drhoxc_av, exc, muxc, dexc, d2exc, dmuxc)  &
!$omp&    private (d2muxc, q_mu, index, imu, rhoxc, drhoxc, n2, l2)   &
!$omp&    private (jssh, jndex, inu, ix)
        do iatom = iatomstart, iatomstart - 1 + natomsp
         in1 = imass(iatom)
         matom = neigh_self(iatom)
 
! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)
          mbeta = neigh_b(ineigh,iatom)
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom)

! Clean derivative matrix
          bcxcpx = 0.0d0
 
! We have ontop left, and ontop right.
! Left is <1|V(1)|2> and Right is <1|V(2)|2>
! Check to make sure we are not doing <1|V(1)|1>. This is done in the atm case.
          if (iatom .eq. jatom .and. mbeta .eq. 0) then
 
! Do nothing here - special case.
 
          else

! Evaluate derivative of average density arho_2c with respect to r2,r1
! First we need the diagonal terms, then we can calculate off-diag terms
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! DIAGONAL TERMS 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! d(arho)/dR = 1/(2l+1) Sum_{-l,+l} d(<mu i|rho_j|mu i>)/dR where:
!      d(<mu i|rho_j|nu i>)/dR  .... denmpx (in CRYSTAL coord!!)
!      d(arho)/dR               .... darho                    
!      darho = 0.0d0
           n1 = 0

! Loop over shell
           do issh = 1, nssh(in1)
            l1 = lssh(issh,in1)
            n1 = n1 + l1 + 1

! ------------------------------------------------------------------------
! Double counting correction forces (dxcv)
! ------------------------------------------------------------------------
            rhoxc_av = arho_on(issh,issh,iatom)
            drhoxc_av(:) = arhop_on(:,issh,issh,ineigh,iatom)
            call cepal(rhoxc_av, exc, muxc, dexc, d2exc, dmuxc, d2muxc)

            q_mu = Qneutral(issh,in1) / (2.0d0*l1 + 1)

! Loop over orbitals in the x-shell
            do index = -l1, l1
             imu = n1 + index
             rhoxc  = rho_on(imu,imu,iatom)
             drhoxc(:) = rhop_on(:,imu,imu,ineigh,iatom)

             bcxcpx(:,imu,imu) =                                             &
     &        ((dexc - dmuxc)*drhoxc_av(:)                                   &
     &          + (d2exc - d2muxc)*(rhoxc - rhoxc_av)*drhoxc_av(:)           &
     &          +  (dexc - dmuxc)*(drhoxc(:) - drhoxc_av(:)) )*q_mu

! Convert to forces (we use '-' in sum)
             dxcdcc(:,ineigh,iatom) =                                        &
     &        dxcdcc(:,ineigh,iatom) -  bcxcpx(:,imu,imu)
            end do !do index = -l1, l1
            n1 = n1 + l1
!-------------------------------------------------------------------------
! End Double counting correction forces (dxcv)
!-------------------------------------------------------------------------
           end do ! end do issh


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! OFF-DIAGONAL TERMS d(<mu i|rho|nu i>)/dR
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! d(arho(shi,shj))/dR = ( d(arho(shi,shi))/dR + d(arho(shj,shj))/dR ) / 2.0
           n1 = 0
           do issh = 1, nssh(in1)
            l1 = lssh(issh,in1)
            n1 = n1 + l1 + 1
            n2 = 0
            do jssh = 1, nssh(in1)
             l2 = lssh(jssh,in1)
             n2 = n2 + l2 + 1
             ! calculate XC potentials 
             rhoxc_av     = arho_on(issh,jssh,iatom)
             drhoxc_av(:) = arhop_on(:,issh,jssh,ineigh,iatom) 
             call cepal(rhoxc_av, exc, muxc, dexc, d2exc, dmuxc, d2muxc)

! Loop over orbitals in the x-shell
             do index = -l1, l1
              imu = n1 + index
! Loop over orbitals in the y-shell
              do jndex = -l2, l2
               inu = n2 + jndex
                       
               rhoxc  = rho_on(imu,inu,iatom)
               drhoxc(:) = rhop_on(:,imu,inu,ineigh,iatom)
! Assemble matrices                       
               do ix = 1,3
                bcxcpx(ix,imu,inu) =                                         &
     &           drhoxc(ix)*dmuxc + rhoxc*drhoxc_av(ix)*d2muxc               &
     &           - delta(imu,inu)*drhoxc_av(ix)*d2muxc*rhoxc_av  

! Convert to forces (we use '-' in sum)
                faxc(ix,ineigh,iatom) = faxc(ix,ineigh,iatom)              &
     &            -  bcxcpx(ix,imu,inu)*rho(imu,inu,matom,iatom)
! JOM 
! add it to the atom-case
!
                gh_atm(ix,imu,inu,ineigh,iatom) =                            &
     &           gh_atm(ix,imu,inu,ineigh,iatom) + bcxcpx(ix,imu,inu)
! JOM-end
               end do ! do ix =1,3
              end do !do jndex = -l2, l2
             end do !do index = -l1, l1
             n2 = n2 + l2
            end do !do jssh = 1, nssh(in1)
            n1 = n1 + l1
           end do ! do issh = 1, nssh(in1)
          end if

! ****************************************************************************
! End loop over iatom and its neighbors - iniegh.
         end do ! end ineigh
        end do ! end iatom

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
        return
        end subroutine Dassemble_olsxc_on_mdet
