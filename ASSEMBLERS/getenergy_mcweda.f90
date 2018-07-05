! copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Institute of Physics, Czech Republic - Pavel Jelinek

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


! readdata.f90
! Program Description
! ===========================================================================
!       This routine calculates the total energy for McWeda approx
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine getenergy_mcweda (itime_step)

        use options 
        use energy
        use mpi_main
        use scf
        use configuration
        use interactions
        use charges, only : efermi
        
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
       integer, intent (in) :: itime_step
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
    
 
! Procedure
! ===========================================================================

! Now call the stuff for the short-range interactions and the ewald terms
! for the final time to get the forces.
         if (itheory .eq. 1 .or. itheory .eq. 2) then
          call get_ewald (nprocs, my_proc, iforce, icluster, itheory,&
     &                    iordern)
         end if
 
! Now call the short-ranged potential to get: u0(iatom,ineigh) and uee00(iatom).
! Here iatom is a basis atom in the unit cell and ineigh is the ineigh'th
! neighbor to iatom.  The total energy per unit cell is
! sum(iatom,ineigh) u0(iatom,ineigh) - sum(iatom) uee00(iatom).
! The energy/cell uii_uee is returned in the calling statment.
         call assemble_usr (itheory, itheory_xc, iforce, uxcdcc_hf,  &
     &                      uiiuee)

! If using the average density approximation or the orbital occupancy
! formalism, then the double-counting correction term must be different.
        if (itheory_xc .eq. 1) then 
           uxcdcc = uxcdcc_sn
           etotxc_1c = 0.0d0
        endif
        if (itheory_xc .eq. 2) then 
! to avoid confusion here we add etotxc_1c to double counting term
! and set etotxc_1c to zero to do not double it in final print
           uxcdcc = uxcdcc_ols + etotxc_1c 
           etotxc_1c = 0.0d0
        endif

! ===========================================================================
         etot = ebs + uiiuee + uxcdcc + etotxc_1c
         if (ivdw .eq. 1) then
          call get_vdw () 
          etot = etot + vdw
         end if

         if (idftd3 .ge. 1) then
          call dftd3_corrections
          etot = etot + etot_dftd3
         end if

         if (iharmonic .eq. 1) then
          call getHarmonic() 
          etot = etot + enHarmonic
         end if
         etot = etot + eqmmm
!         write(*,*) 'eqmmm =', eqmmm
         etotper = etot/natoms
!         please leave this until I finish--Brandon
!         etotPerfectLattice = -322.813965
!         sys1En = (1.0 - internalLambda(tau(itime_step)))*(etotper - etotPerfectLattice)
!         sys2En = internalLambda(tau(itime_step))*enHarmonic
!         tiEn = etotPerfectLattice + sys1En + sys2En
!         print*,'tiEn',tiEn
!         print*,'tau(itime_step)',tau(itime_step)
!         print*,'enHarmonic',enHarmonic
!         print*,'etotPerfectLattice',etotPerfectLattice
 
         if (wrtout) then
          write (*,*) ' ---------- T H E  T O T A L  E N E R G Y ----------- '
          if (itheory .ne. 0) write (*,500) itime_step, Kscf 
          if (itheory .eq. 0) write (*,501) itime_step

          write (*,*) '  '
          write (*,502) ebs
          write (*,503) uiiuee
          write (*,504) etotxc_1c
          write (*,505) uxcdcc
          if (ivdw .eq. 1) write (*,506) vdw
          write (*,507) etot
          write (*,508) etotper
          write (*,509) atomic_energy
          write (*,510) etot - atomic_energy
          write (*,512) efermi
          write (*,*) '  '
          write (*,511) (etot - atomic_energy)/natoms
          write (*,*) ' ----------------------------------------------------- '
         end if
         etotold = etotnew
         etotnew = etotper
!-----------------------------------------------------------------------------
!export energies for thermodynamic integration of <dE/dlambda>
!----------------------------------------------------------------------------
         if (ithermoint .eq. 1) then
          open(unit=400,file='energyValues.txt',status='unknown')
          write(400,*)etot
          close(400)
         endif

! Deallocate Arrays
! ===========================================================================

 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
500     format (2x, ' Time step = ', i6, ' SCF step = ', i3)
501     format (2x, ' Time step = ', i6)
502     format (2x, '           ebs = ', f15.6)
503     format (2x, '     uii - uee = ', f15.6)
504     format (2x, '     etotxc_1c = ', f15.6)
505     format (2x, '        uxcdcc = ', f15.6)
506     format (2x, '           vdw = ', f15.6)
507     format (2x, '          ETOT = ', f15.6)
508     format (2x, '     Etot/atom = ', f15.6)
509     format (2x, ' Atomic Energy = ', f15.6)
510     format (2x, '     CohesiveE = ', f15.6)
511     format (2x, ' Cohesive Energy per atom  = ', f15.6)
512     format (2x, '   Fermi Level = ', f15.6)

        return
        end subroutine getenergy_mcweda
 
