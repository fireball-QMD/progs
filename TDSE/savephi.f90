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


! psi2phi.f90
! Program Description
! ===========================================================================
!       This routine projects coefficients of the actual time-dependent wave-function
! onto eigenfunction for given Rn (ionic position)
! ===========================================================================
! Code rewritten by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine getpsi2phi ()

        use tdse
        use density
        use dimensions
        use interactions
        use neighbor_map
        use configuration
        use kpoints

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
!       integer, intent (in) :: itime_step

! Local Parameters and Data Declaration
! ===========================================================================

       integer  ielec
       integer  iatom
       integer  jatom
       integer  ineigh
       integer  iband
       integer  imu
       integer  jmu
       integer  inu
       integer  jnu
       integer  ikpoint
       integer  in1
       integer  in2
       complex  a0
       complex  phi
       complex, dimension (norbitals) :: spsi


! Local Variable Declaration and Description
! ===========================================================================


! Procedure
! ===========================================================================

        a0 = cmplx(0.0d0,0.0d0)

        do ikpoint = 1, nkpoints
         do ielec = 1, nelec
! get S*psi
          spsi = a0
          do iatom = 1, natoms
           in1 = imass(iatom)
           do ineigh = 1, neighn (iatom)
            jatom = neigh_j(iatom, ineigh)
            in2 = imass(jatom)
            do imu = 1, num_orb(in1)
             jmu = imu + degelec(iatom)
             do inu = 1, num_orb(in2)
              jnu = inu + degelec(jatom)
              spsi(jmu) = spsi(jmu) + s_mat(imu,inu,iatom,ineigh)*psi(jnu,ielec,ikpoint)
             enddo ! do inu
            enddo ! do imu
           enddo ! ineigh
          enddo ! iatoms
! get phi+ * (S*psi)
          do iband = 1,norbitals_new
           do imu = 1, norbitals
            phi = cmplx(bbnkre(imu,iband,ikpoint),-1.0d0*bbnkim(imu,iband,ikpoint))
            psi2phi(iband,ielec,ikpoint) = psi2phi(iband,ielec,ikpoint)*phi*spsi(imu)
           enddo ! imu
! writeout control output
           write (*,*) 'k-point no.',ikpoint
           write (*,*) 'projection of ',ielec,'-e onto',iband,'-band =',psi2phi(iband,ielec,ikpoint)
           write (*,100)
          enddo ! iband
         enddo ! ielec
        enddo ! ikpoints

! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine getpsi2phi

