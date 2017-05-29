! copyright info:
!
!                             @Copyright 2005
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Universidad de Madrid - Jose Ortega
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


! diag_k.f90
! Program Description
! ===========================================================================
!       This is  a subroutine of k-loop
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
! Program Declaration
! ===========================================================================
        subroutine diag_k ( )

        use configuration
        use density
        use dimensions
        use interactions
        use neighbor_map
        use kpoints
        use options
        use outputs
        use mpi_main
        use transport
        use scf
        use hartree_fock
        use module_dos
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer ikpoint
        integer imu
        real, dimension (3) :: k_temp

! Procedure
! ===========================================================================

! Now we have the real space hamiltonian and overlap. Compute the k-space
! Hamiltonian and overlap and diagonalize. First put the k-points into the
! KPOINTS file.
        if (iwrteigen .eq. 1) then
          open (unit = 19, file = 'eigen.dat', status = 'unknown')
          open (unit = 20, file = 'ek.dat', status = 'unknown')
        endif

! The eigenvectors of the overlap are always computed.
! Allocate the arrays for wavefunction coefficients
!       if (itdse .eq. 1) then
! JOM
! jel-nac
!        if (itdse .eq. 1 .or. imdet .eq. 1) then
        if (itdse .eq. 1 .or. imdet .ne. 0 .or. icDFT .eq. 1) then
! nothing to do, we already allocated these arrays
! JOM-info : for mdet they are allocated in init_mdet.f90
          write (*,*) ' save eigenstuff '
        else
         if (iqout .ne. 2) allocate (blowre (norbitals, norbitals, nkpoints))
         if (iqout .ne. 2 .and. icluster .ne. 1)                           &
     &      allocate (blowim (norbitals, norbitals, nkpoints))
         allocate (bbnkre (norbitals, norbitals, nkpoints))
         if (icluster .ne. 1)                                              &
     &      allocate (bbnkim (norbitals, norbitals, nkpoints))
         allocate (eigen_k (norbitals, nkpoints))
        endif

! reset Green func. for DOS
        if (iwrtdos .eq. 1) green = (0.0d0, 0.0d0)

!$omp parallel do private (k_temp)
        do ikpoint = 1, nkpoints

! The subroutine kspace wants the k-vector in inverse angstrom units.
! NOT pi/alat units.
         k_temp(:) = special_k(:,ikpoint)

         call kspace (nprocs, my_proc, Kscf, iqout, icluster,  &
     &                iwrteigen, ikpoint, k_temp,       &
     &                nkpoints, iwrtdos, iwrthop, iwrtatom, itrans, igap)

! CGP
! if iwrtdos greatter than 1 then we calculate the DOS CGP
         if (iwrtdos .ge. 1) then
           ! GAP ENRIQUE-FF
           if ( .not.((igap .eq. 1).and.(Kscf .eq. 1)) ) then
           ! end GAP ENRIQUE-FF
              call dos (weight_k(ikpoint), k_temp)
           ! GAP ENRIQUE-FF
           end if
           ! end GAP ENRIQUE-FF
         endif
         if (iwrtatom .ge. 1) then
          call hamilt_atom (weight_k(ikpoint), k_temp)
         endif
! end CGP

         if (itrans .eq. 1) then
          call assemble_t12_bare (k_temp, ikpoint, weight_k(ikpoint))
         endif

! Write out the eigenvalues for the band structure.  The first time write out
! nkpoints, and norbitals.
         if (ikpoint .eq. 1) then
!         if (iPRESSURE .eq. 1) write (18,*) nkpoints
          if (iwrteigen .eq. 1) write (19,*) nkpoints, norbitals_new
         end if
!        if (iPRESSURE .eq. 1)                                            &
!    &       write (18,100) special_k(:,ikpoint), weight_k(ikpoint)

         if (iwrteigen .eq. 1) then
           write (19,*) ' ------ the energy eigenvalues ----'
           write (19,100) (eigen_k(imu,ikpoint), imu = 1, norbitals_new)
           write (20,101, advance="no") ikpoint
           do imu = 1, norbitals_new
            write (20,102, advance="no") eigen_k(imu,ikpoint)
           enddo
           write (20,*)
!          write (20,101) ikpoint,(eigen_k(imu,ikpoint), imu=1,norbitals_new)
         end if

        end do ! do ikpoint

! Close output file
!        close (unit = 18)
        if (iwrteigen .eq. 1) then
         close (unit = 19)
         close (unit = 20)
        endif


! Format Statements
! ===========================================================================
100     format (2x, 4(2x,f11.5))
101     format (i4)
102     format (f11.5)



        return
      end subroutine diag_k
