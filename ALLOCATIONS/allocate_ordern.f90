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

 
! allocate_ordern.f90
! Program Description
! ===========================================================================
!       These allocations are for the slave process routine.
 
! ===========================================================================
! Code written by:
! Spencer Shellman 
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine allocate_ordern (natoms, nspecies)
        use charges
        use constants_fireball
        use integrals
        use interactions
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms
        integer, intent (in) :: nspecies

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
        allocate (Qneutral (nsh_max, nspecies))
        allocate (Qin (nsh_max, natoms))
        allocate (Qout (nsh_max, natoms))
        allocate (QLowdin_TOT (natoms))
        allocate (QMulliken_TOT (natoms))
        allocate (dq (nspecies))
        allocate (mu2shell (numorb_max, nspecies))
        allocate (xcnu1c (nsh_max, nsh_max, nspecies))
        allocate (num_orbPP (nspecies))
        allocate (nsshPP (nspecies))
        allocate (icon3c (nspecies, nspecies, nspecies))
        allocate (numx3c_bcna (0:isorpmax, nspecies**3))
        allocate (numy3c_bcna (0:isorpmax, nspecies**3))
        allocate (hx_bcna (0:isorpmax, nspecies**3))
        allocate (hy_bcna (0:isorpmax, nspecies**3))
        allocate (x3cmax_bcna (0:isorpmax, nspecies**3))
        allocate (y3cmax_bcna (0:isorpmax, nspecies**3))
        allocate (z2cmax (interactions2c_max, nspecies, nspecies))
        allocate (numz2c (interactions2c_max, nspecies, nspecies))
        allocate (numx3c_xc3c (0:ideriv_max, nspecies**3))
        allocate (numy3c_xc3c (0:ideriv_max, nspecies**3))
        allocate (hx_xc3c (0:ideriv_max, nspecies**3))
        allocate (hy_xc3c (0:ideriv_max, nspecies**3))
        allocate (x3cmax_xc3c (0:ideriv_max, nspecies**3))
        allocate (y3cmax_xc3c (0:ideriv_max, nspecies**3))
        allocate (index_max2c (nspecies, nspecies))
        allocate (index_max3c (nspecies, nspecies))
        allocate (index_maxPP (nspecies, nspecies))
        allocate (lsshPP (nsh_max, nspecies))
        allocate (cl_PP (0:nsh_max - 1, nspecies))

        allocate (bcna_01 (numXmax, numYmax, ME3c_max, 0:isorpmax, nspecies**3))
        allocate (bcna_02 (numXmax, numYmax, ME3c_max, 0:isorpmax, nspecies**3))
        allocate (bcna_03 (numXmax, numYmax, ME3c_max, 0:isorpmax, nspecies**3))
        allocate (bcna_04 (numXmax, numYmax, ME3c_max, 0:isorpmax, nspecies**3))
        allocate (bcna_05 (numXmax, numYmax, ME3c_max, 0:isorpmax, nspecies**3))
        allocate (xc3c_01 (numXmax, numYmax, ME3c_max, 0:ideriv_max,         &
     &                     nspecies**3))
        allocate (xc3c_02 (numXmax, numYmax, ME3c_max, 0:ideriv_max,         &
     &                     nspecies**3))
        allocate (xc3c_03 (numXmax, numYmax, ME3c_max, 0:ideriv_max,         &
     &                    nspecies**3))
        allocate (xc3c_04 (numXmax, numYmax, ME3c_max, 0:ideriv_max,         &
     &                    nspecies**3))
        allocate (xc3c_05 (numXmax, numYmax, ME3c_max, 0:ideriv_max,         &
     &                    nspecies**3))

        if (.not. superspline)                                               &
     &   allocate (xintegral_2c (ME2c_max, nfofx, interactions2c_max,        &
     &                              nspecies, nspecies))
        if (superspline)                                                     & 
         allocate (splineint_2c(1:4, ME2c_max, nfofx, interactions2c_max,    &
     &                          nspecies, nspecies))

        allocate (mu (ME3c_max, nspecies, nspecies))
        allocate (mvalue (ME3c_max, nspecies, nspecies))
        allocate (nu (ME3c_max, nspecies, nspecies))
        allocate (muPP (ME3c_max, nspecies, nspecies))
        allocate (nuPP (ME3c_max, nspecies, nspecies))

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
        return
        end
