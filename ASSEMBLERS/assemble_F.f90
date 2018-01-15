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


! assemble_F.f90
! Program Description
! ===========================================================================
!       This program assembles all of the forces and writes them out
! if desired.
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
! Program Declaration
! ===========================================================================
        subroutine assemble_F (natoms, itheory, itheory_xc, igauss, ivdw,    &
       &                       iharmonic, ibias, iwrt_fpieces)
        use dimensions
        use forces
        use neighbor_map
        use options, only : idftd3
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: igauss
        integer, intent(in) :: itheory
        integer, intent(in) :: itheory_xc
        integer, intent(in) :: ivdw
        integer, intent(in) :: iharmonic
        integer, intent(in) :: ibias
        integer, intent(in) :: iwrt_fpieces
        integer, intent(in) :: natoms

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ineigh
        integer jatom

        integer ix

        real rms
        real maximum
        real minimum
        real, dimension (3, natoms) :: f3ca       ! three-center charged atom
        real, dimension (3, natoms) :: f3na       ! three-center neutral atom
        real, dimension (3, natoms) :: f3nl       ! three-center non-local
        real, dimension (3, natoms) :: f3xc       ! three-center xc
        real, dimension (3, natoms) :: f3xc_ca    ! three-center xc - charges
        real, dimension (3, natoms) :: fbs        ! total band-structure force
        real, dimension (3, natoms) :: fca        ! total charged atom force
        real, dimension (3, natoms) :: fcaatm     ! charged atom/atom force
        real, dimension (3, natoms) :: fcaot      ! charged atom/ontop force
        real, dimension (3, natoms) :: fna        ! total neutral atom force
        real, dimension (3, natoms) :: fnaatm     ! neutral atom/atom force
        real, dimension (3, natoms) :: fnaot      ! neutral atom/ontop force
        real, dimension (3, natoms) :: fnl        ! total non-local force
        real, dimension (3, natoms) :: fnlatm     ! non-local/atom force
        real, dimension (3, natoms) :: fnlot      ! non-local/ontop force
        real, dimension (3, natoms) :: fxc        ! total xc force
        real, dimension (3, natoms) :: fxc_ca     ! total xc force - charges
        real, dimension (3, natoms) :: fxcatm     ! xc atm/atm force
        real, dimension (3, natoms) :: fxcatm_ca  ! xc atm/atm force - charges
        real, dimension (3, natoms) :: fxcot      ! xc ontop force
        real, dimension (3, natoms) :: fxcot_ca   ! xc ontop force - charges

! Procedure
! ===========================================================================
        write (*,*) ' Welcome to assemble_F - ftot assembled here. '

! ****************************************************************************
! Assemble three-center forces from Dassemble_3c:
! f3naa, f3nab, f3nac, f3nla, f3nlb, f3nlc, f3xca, f3xcb, f3xcc
! ****************************************************************************
        do iatom = 1, natoms

! non-local pseudopotential interactions
         f3nl(:,iatom) = f3nla(:,iatom) + f3nlb(:,iatom) + f3nlc(:,iatom)

! neutral atom interactions
         f3na(:,iatom) = f3naa(:,iatom) + f3nab(:,iatom) + f3nac(:,iatom)

! charged atom interactions
         if (itheory .eq. 1) then
          f3ca(:,iatom) = f3caa(:,iatom) + f3cab(:,iatom) + f3cac(:,iatom)
         else
          f3ca(:,iatom) = 0.0d0
         end if

! exchange-correlation interactions
         f3xc(:,iatom) = f3xca(:,iatom) + f3xcb(:,iatom) + f3xcc(:,iatom)

! charged exchange-correlation interactions
         if (itheory .eq. 1) then
          f3xc_ca(:,iatom) = f3xca_ca(:,iatom) + f3xcb_ca(:,iatom)           &
     &                      + f3xcc_ca(:,iatom)
         else
          f3xc_ca(:,iatom) = 0.0d0
         end if
        end do

! ****************************************************************************
! Compute the atomic terms and the ontop terms.
! Initialize fnaatm, fnaot, fnlatm, fnlot, fxcatm, and fxcot.
! ****************************************************************************
        fcaatm = 0.0d0
        fnaatm = 0.0d0
        fcaot = 0.0d0
        fnaot = 0.0d0
        fxcatm = 0.0d0
        fxcatm_ca = 0.0d0
        fxcot = 0.0d0
        fxcot_ca = 0.0d0
        fnlatm = 0.0d0
        fnlot = 0.0d0
        if((itheory_xc.eq.1) .or. (itheory_xc.eq.2)) dxcv = 0.0d0
! ****************************************************************************
! Loop over all atoms iatom in the central cell.
! Loop over all neighbors of iatom
        do iatom = 1, natoms
         do ineigh = 1, neighn(iatom)
          jatom = neigh_j(ineigh,iatom)

! There may be a confusion here about the signs.
! The variable "sum" is already force-like because it came from fana, etc.
! another +/- comes from direct (+) and cross (-) terms.
! See p. 5, "the total band structure force".
! The atom terms.
! neutral atom forces - atm case
          fnaatm(:,iatom) = fnaatm(:,iatom) + fana(:,ineigh,iatom)
          fnaatm(:,jatom) = fnaatm(:,jatom) - fana(:,ineigh,iatom)

          if (itheory .eq. 1) then
           fcaatm(:,iatom) = fcaatm(:,iatom) + faca(:,ineigh,iatom)
           fcaatm(:,jatom) = fcaatm(:,jatom) - faca(:,ineigh,iatom)
          end if

! exchange-correlation forces - atm case
! JOM info : the notation in Dassemble_olsxc_2c and Dassemble_olsxc_on
! was misleading. I have now corrected it.               
! we now use "faxc" in Dassemble_olsxc_on (atm case)
! and "fotxc" in Dassemble_olsxc_2c (ontop cases)
          fxcatm(:,iatom) = fxcatm(:,iatom) + faxc(:,ineigh,iatom)
          fxcatm(:,jatom) = fxcatm(:,jatom) - faxc(:,ineigh,iatom)

! charged exchange-correlation forces - atm case
          if (itheory .eq. 1) then
           fxcatm_ca(:,iatom) = fxcatm_ca(:,iatom) + faxc_ca(:,ineigh,iatom)
           fxcatm_ca(:,jatom) = fxcatm_ca(:,jatom) - faxc_ca(:,ineigh,iatom)
          end if

! exchange-correlation three-center interactions - gaussian approximation
          if (igauss .eq. 1) then
           f3xc(:,iatom) = f3xc(:,iatom) + fxcro(:,ineigh,iatom)
           f3xc(:,jatom) = f3xc(:,jatom) - fxcro(:,ineigh,iatom)
          end if

! The ontop terms.
! The factor 2.0d0 comes from ontop the bra or ontop the ket.
! neutral atom forces - ontop case
          fnaot(:,iatom) = fnaot(:,iatom) + 2.0d0*fotna(:,ineigh,iatom)
          fnaot(:,jatom) = fnaot(:,jatom) - 2.0d0*fotna(:,ineigh,iatom)

! charged atom forces - ontop case
          if (itheory .eq. 1) then
           fcaot(:,iatom) = fcaot(:,iatom) + 2.0d0*fotca(:,ineigh,iatom)
           fcaot(:,jatom) = fcaot(:,jatom) - 2.0d0*fotca(:,ineigh,iatom)
          end if

! exchange-correlation - ontop case
! JOM info : the notation in Dassemble_olsxc_2c and Dassemble_olsxc_on
! was misleading. I have now corrected it.               
! we now use "faxc" in Dassemble_olsxc_on (atm case)
! and "fotxc" in Dassemble_olsxc_2c (ontop cases)
! Does this term have a factor 2.0d0 ? No :
! The ontop-left and on-top-right contributions to the derivative
! of the matrix elements have already been added, so no 2.0d0 factor
          fxcot(:,iatom) = fxcot(:,iatom) + fotxc(:,ineigh,iatom)
          fxcot(:,jatom) = fxcot(:,jatom) - fotxc(:,ineigh,iatom)

! charged exchange-correlation - ontop case
          if (itheory .eq. 1) then
           fxcot_ca(:,iatom) = fxcot_ca(:,iatom) + fotxc_ca(:,ineigh,iatom)
           fxcot_ca(:,jatom) = fxcot_ca(:,jatom) - fotxc_ca(:,ineigh,iatom)
          end if

          if((itheory_xc.eq.1) .or. (itheory_xc.eq.2)) then
! double counting correction (only OLSXC && SNXC)
           dxcv(:,iatom) = dxcv(:,iatom) + dxcdcc(:,ineigh,iatom)
           dxcv(:,jatom) = dxcv(:,jatom) - dxcdcc(:,ineigh,iatom)

          endif

         end do     ! end loop over neighbors
        end do      ! end loop over atoms


! ****************************************************************************
!    PSEUDOPOTENTIAL  2-CENTER   PART
! ****************************************************************************
! Loop over all atoms iatom in the central cell.
! Loop over all PP-neighbors of iatom
        do iatom = 1, natoms
         do ineigh = 1, nPPn(iatom)
          jatom = nPP_j(ineigh,iatom)

! There may be a confusion here about the signs.
! The variable "sum" is already force-like because it came from fana, etc.
! another +/- comes from direct (+) and cross (-) terms.
! See p. 5, "the total band structure force".
! The atom terms.
! non-local forces - atm case
          fnlatm(:,iatom) = fnlatm(:,iatom) + fanl(:,ineigh,iatom)
          fnlatm(:,jatom) = fnlatm(:,jatom) - fanl(:,ineigh,iatom)

         end do     ! end loop over PP-neighbors
        end do      ! end loop over atoms

! Loop over all PP-neighbors of iatom
        do iatom = 1, natoms
         do ineigh = 1, nPPxn(iatom)
          jatom = nPPx_j(ineigh,iatom)

! The ontop terms.
! The factor 2.0d0 comes from ontop the bra or ontop the ket.
! non-local forces - ontop case
          fnlot(:,iatom) = fnlot(:,iatom) + 2.0d0*fotnl(:,ineigh,iatom)
          fnlot(:,jatom) = fnlot(:,jatom) - 2.0d0*fotnl(:,ineigh,iatom)
         end do     ! end loop over PP-neighbors
        end do      ! end loop over atoms

! Now put together all the pieces of the NA, XC, and NL forces.
! Form fna = f3na + fnaatm + fnaot, and similarly for the others.
        do iatom = 1, natoms

! charged atom total force
         if (itheory .eq. 1) then
          fca(:,iatom) = f3ca(:,iatom) + fcaatm(:,iatom) + fcaot(:,iatom)
         else if (itheory .eq. 2) then
          fca(:,iatom) = fcoulomb(:,iatom)
         else
          fca(:,iatom) = 0.0d0
         end if

! neutral atom total force
         fna(:,iatom) = f3na(:,iatom) + fnaatm(:,iatom) + fnaot(:,iatom)
! non-local total force
         fnl(:,iatom) = f3nl(:,iatom) + fnlatm(:,iatom) + fnlot(:,iatom)

! exchange-correlation total force
         fxc(:,iatom) = f3xc(:,iatom) + fxcatm(:,iatom) + fxcot(:,iatom)

! charged atom exchange-correlation total force
         if (itheory .eq. 1) then
          fxc_ca(:,iatom) = fxcatm_ca(:,iatom) + fxcot_ca(:,iatom)           &
     &                     + f3xc_ca(:,iatom)
         else if (itheory .eq. 2) then
          fxc_ca(:,iatom) = fxcnu(:,iatom)
         else
          fxc_ca(:,iatom) = 0.0d0
         end if
        end do


! ****************************************************************************
! If iwrt_fpieces = 1, then write out all of the pieces of the
! band-structure force.
! ****************************************************************************
        if (iwrt_fpieces .eq. 1) then
         write (*,*) '  '
         write (*,*) ' ******************************************************* '
         write (*,*) ' Write out the contributions to the band-structure force.'
         write (*,*) ' ******************************************************* '
         write (*,*) '  '
         write (*,*) ' The kinetic force: '
         do iatom = 1, natoms
          write (*,101) iatom, ft(:,iatom)
         end do

         write (*,*) '  '
         write (*,*) ' The neutral atom force: '
         do iatom = 1, natoms
          write (*,102) iatom, fna(:,iatom)
         end do

         write (*,*) '  '
         write (*,*) ' The neutral atom force ATM : '
         do iatom = 1, natoms
          write (*,108) iatom, fnaatm(:,iatom)
         end do

         write (*,*) '  '
         write (*,*) ' The neutral atom force ONTop : '
         do iatom = 1, natoms
          write (*,109) iatom, fnaot(:,iatom)
         end do

         write (*,*) '  '
         write (*,*) ' The neutral atom force, 3-center: '
         do iatom = 1, natoms
          write (*,301) iatom, f3na(:,iatom)
         end do

         if (itheory .eq. 1 .or. itheory .eq. 2) then
          write (*,*) '  '
          write (*,*) ' The charged atom force: '
          do iatom = 1, natoms
           write (*,103) iatom, fca(:,iatom)
          end do

          write (*,*) '  '
          write (*,*) ' The 3c charged atom force: '
          do iatom = 1, natoms
           write (*,302) iatom, f3ca(:,iatom)
          end do

          write (*,*) '  '
          write (*,*) ' The atomic charged atom force: '
          do iatom = 1, natoms
           write (*,303) iatom, fcaatm(:,iatom)
          end do

          write (*,*) '  '
          write (*,*) ' The ontop charged atom force: '
          do iatom = 1, natoms
           write (*,304) iatom, fcaot(:,iatom)
          end do
         end if

         write (*,*) '  '
         write (*,*) ' The non-local force: '
         do iatom = 1, natoms
          write (*,104) iatom, fnl(:,iatom)
         end do

         if (itheory .eq. 0) then
          write (*,*) '  '
          write (*,*) ' The exchange-correlation force: '
          do iatom = 1, natoms
           write (*,105) iatom, fxc(:,iatom)
          end do

          write (*,*) '  '
          write (*,*) ' The exchange-correlation force, 3c-center: '
          do iatom = 1, natoms
           write (*,305) iatom, f3xc(:,iatom)
          end do

          write (*,*) '  '
          write (*,*) ' The exchange-correlation force atm: '
          do iatom = 1, natoms
           write (*,306) iatom, fxcatm(:,iatom)
          end do

          write (*,*) '  '
          write (*,*) ' The exchange-correlation force ontop: '
          do iatom = 1, natoms
           write (*,307) iatom, fxcot(:,iatom)
          end do
         end if

         if (itheory .eq. 1 .or. itheory .eq. 2) then
            write (*,*) '  '
            write (*,*) ' The charged atom exchange-correlation force: '
            do iatom = 1, natoms
               write (*,106) iatom, fxc_ca(:,iatom)
            end do

            write (*,*) '  '
            write (*,*) ' The exchange-correlation force 3c: '
            do iatom = 1, natoms
               write (*,308) iatom, f3xc_ca(:,iatom)
            end do

            write (*,*) '  '
            write (*,*) ' The exchange-correlation force atm: '
            do iatom = 1, natoms
               write (*,309) iatom, fxcatm_ca(:,iatom)
            end do

            write (*,*) '  '
            write (*,*) ' The exchange-correlation force ontop: '
            do iatom = 1, natoms
               write (*,310) iatom, fxcot_ca(:,iatom)
            end do

         end if

         if (itheory .eq. 1 .or. itheory .eq. 2) then
          write (*,*) '  '
          write (*,*) ' The long-range electrostatic force: '
          do iatom = 1, natoms
           write (*,107) iatom, flrew(:,iatom)
          end do
         end if
         write (*,*) ' *********************************************** '
         write (*,*) '  '
        end if

! ****************************************************************************
! Put together the entire bandstructure force. Add in the LR BS force.
! ****************************************************************************
        do iatom = 1, natoms
         fbs(:,iatom) = ft(:,iatom) + fna(:,iatom) + fnl(:,iatom) + fxc(:,iatom)
         if (itheory .eq. 1 .or. itheory .eq. 2) then
          fbs(:,iatom) = fbs(:,iatom) + fca(:,iatom) + fxc_ca(:,iatom)       &
     &                                + flrew(:,iatom) + flrew_qmmm(:,iatom)
         end if
        end do
        if (itheory .eq. 3) fbs = 0.0d0

! ****************************************************************************
! If iwrt_fpieces = 1, then write out all of the pieces of the total force.
! ****************************************************************************
        if (iwrt_fpieces .eq. 1) then
         write (*,*) '  '
         write (*,*) ' *********************************************** '
         write (*,*) ' Write out the contributions to the total force. '
         write (*,*) ' *********************************************** '
         write (*,*) '  '
         write (*,*) ' The band-structure force: '
         do iatom = 1, natoms
          write (*,200) iatom, fbs(:,iatom)
         end do

         write (*,*) '  '
         write (*,*) ' The short-range force: '
         do iatom = 1, natoms
          write (*,201) iatom, dusr(:,iatom)
         end do

         write (*,*) '  '
         write (*,*) ' The exchange-correlation double-counting force: '
         do iatom = 1, natoms
          write (*,202) iatom, dxcv(:,iatom)
         end do

         write (*,*) '  '
         write (*,*) ' The overlap-repulsive force: '
         do iatom = 1, natoms
          write (*,203) iatom, fro(:,iatom)
         end do

         if (ivdw .eq. 1) then
          write (*,*) '  '
          write (*,*) ' The van der Waals force: '
          do iatom = 1, natoms
           write (*,204) iatom, fvdw(:,iatom)
          end do
         end if

         if (iharmonic .eq. 1) then
          write (*,*) '  '
          write (*,*) ' The external field force: '
          do iatom = 1, natoms
           write (*,205) iatom, fharmonic(:,iatom)
          end do
         end if

         if (ibias .eq. 1) then
          write (*,*) '  '
          write (*,*) ' The bias field force: '
          do iatom = 1, natoms
           write (*,206) iatom, fbias(:,iatom)
          end do
         end if


        end if

! ****************************************************************************
! Now the total force: ftot
! ****************************************************************************
        write (*,*) '  '
        write (*,*) ' The grand total force (eV/A): '

        do iatom = 1, natoms
         ftot(:,iatom) = fbs(:,iatom) + dusr(:,iatom) + dxcv(:,iatom)        &
     &                                + fro(:,iatom)
         if (ivdw .eq. 1) ftot(:,iatom) = ftot(:,iatom) + fvdw(:,iatom)
         if (idftd3 .ge. 1) ftot(:,iatom) = ftot(:,iatom) + ftot_dftd3(:,iatom)
         if (iharmonic .eq. 1) ftot(:,iatom) = ftot(:,iatom)                 &
     &                                        + fharmonic(:,iatom)
         if (ibias .eq. 1) ftot(:,iatom) = ftot(:,iatom)                     &
     &                                        + fbias(:,iatom)
         write (*,100) iatom, ftot(:,iatom)
        end do

        rms = 0.0d0
        do iatom = 1, natoms
         do ix = 1,3
          rms = rms + ftot(ix,iatom)**2
         end do
        end do
        rms = sqrt(rms/(3*natoms))

        maximum = maxval(ftot)
        minimum = abs(minval(ftot))
        if (minimum .gt. maximum) maximum = minimum

        write (*,*) '  '
        write (*,400) maximum, rms
        write (*,*) '  '

! Format Statements
! ===========================================================================
100     format (2x, ' iatom = ', i4, ' ftot      = ', 3e14.6)
101     format (2x, ' iatom = ', i4, ' ft        = ', 3e14.6)
102     format (2x, ' iatom = ', i4, ' fna       = ', 3e14.6)
108     format (2x, ' iatom = ', i4, ' fnaatm    = ', 3e14.6)
109     format (2x, ' iatom = ', i4, ' fnaot     = ', 3e14.6)
103     format (2x, ' iatom = ', i4, ' fca       = ', 3e14.6)
104     format (2x, ' iatom = ', i4, ' fnl       = ', 3e14.6)
105     format (2x, ' iatom = ', i4, ' fxc       = ', 3e14.6)
106     format (2x, ' iatom = ', i4, ' fxc_ca    = ', 3e14.6)
107     format (2x, ' iatom = ', i4, ' flrew     = ', 3e14.6)
200     format (2x, ' iatom = ', i4, ' fbs       = ', 3e14.6)
201     format (2x, ' iatom = ', i4, ' dusr      = ', 3e14.6)
202     format (2x, ' iatom = ', i4, ' dxcv      = ', 3e14.6)
203     format (2x, ' iatom = ', i4, ' fro       = ', 3e14.6)
204     format (2x, ' iatom = ', i4, ' fvdw      = ', 3e14.6)
205     format (2x, ' iatom = ', i4, ' fharmonic = ', 3e14.6)
206     format (2x, ' iatom = ', i4, ' fbias     = ', 3e14.6)
301     format (2x, ' iatom = ', i4, ' f3na      = ', 3e14.6)
302     format (2x, ' iatom = ', i4, ' f3ca      = ', 3e14.6)
303     format (2x, ' iatom = ', i4, ' fcaatm    = ', 3e14.6)
304     format (2x, ' iatom = ', i4, ' fcaot     = ', 3e14.6)
305     format (2x, ' iatom = ', i4, ' f3xc      = ', 3e14.6)
306     format (2x, ' iatom = ', i4, ' fxcatm    = ', 3e14.6)
307     format (2x, ' iatom = ', i4, ' fxcot     = ', 3e14.6)
308     format (2x, ' iatom = ', i4, ' f3xc_ca   = ', 3e14.6)
309     format (2x, ' iatom = ', i4, ' fxcatm_ca = ', 3e14.6)
310     format (2x, ' iatom = ', i4, ' fxcot_ca  = ', 3e14.6)

400     format (2x, ' Cartesian Forces:  Max = ', f16.8, '    RMS = ', f16.8)

        return
        end
