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


! wrtout_psiT.f90
! Program Description
! ===========================================================================
!       This subroutine writes out an actual time snapshot of wavefunction
! into files for each electron
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
        subroutine wrtout_psiT (itime_step, ietime)

        use tdse
        use dimensions
        use interactions
        use kpoints
        use MD
!rhoT
        use density
        use neighbor_map
        use configuration

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
       integer, intent (in) :: itime_step
       integer, intent (in) :: ietime

! Local Parameters and Data Declaration
! ===========================================================================


! Local Variable Declaration and Description
! ===========================================================================
       integer  ielec
       integer  imu
       integer  jmu
       integer  inu
       integer  ikpoint
!rho T
       integer iatom
       integer jatom
       integer ineigh
       integer in1
       integer in2

       character(4)    :: ext
       character(80)   :: fname

       real aux1

       complex, dimension (norbitals) :: psix
       complex a0

! Procedure
! ===========================================================================
! NOTE: s12psi has been calculated in the subroutine tddenmat for
! former ionic configuration
       a0 = cmplx(0.0d0, 0.0d0)

! Loop over electrons
       do ielec = 1, nelec

! project wf from reciprocal to real space
        psix = a0
! Loop over the special k points.
        do ikpoint = 1, nkpoints
         aux1 = weight_k(ikpoint)
! Finally the imu loop
         imu = 0
         do imu = 1, norbitals
          psix(imu) = psix(imu) + aux1*psi(imu,ielec,ikpoint)
         end do ! imu
! End loop over orbitals and kpoints
        end do ! ikpoint

! create an unique filename for each ielec
        write (ext,'(i4.4)') ielec
        fname = 'T_psi_'//ext//'.dat'

! open output file
        if ((itime_step .eq. 1) .and. (ietime .eq. 1)) then
         open (unit = 17, file = fname, status = 'unknown')
        else
         open (unit = 17, file = fname, status = 'unknown',        &
    &            position = 'append')
        end if
        write (17,201,advance="no") itime_step    
         do imu=1,norbitals
           write (17,200,advance="no") real(psi(imu,ielec,ikpoint))**2
         enddo
          write (17,*) ' '
        close (17)

       enddo ! do ielec

! write out rho
       if ((itime_step .eq. 1) .and. (ietime .eq. 1)) then
         open (unit = 20, file = 'rhoT.dat', status = 'unknown')
       else
         open (unit = 20, file = 'rhoT.dat', status = 'unknown',        &
    &            position = 'append')
       end if
       write (20,801) real((itime_step-1)*dt)+real(dte*ietime)
       do iatom = 1,natoms
         in1 = imass(iatom)
         do ineigh = 1, neighn(iatom)
          jatom = neigh_j(ineigh,iatom)
          in2 = imass(jatom)
          do imu = 1,num_orb(in1)
           do inu = 1, num_orb(in2)
            write (20,800) iatom, ineigh, imu, inu, rho(imu,inu,ineigh,iatom)
           enddo
          enddo
         enddo
       enddo
       close (20)
! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))
200     format (f12.4)
201     format (i6)
800     format (2x,4i3,f12.6)
801     format (' dt = ',f16.8)

        return
        end subroutine wrtout_psiT

