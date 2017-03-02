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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! tddiag_k.f90
! Program Description
! ===========================================================================
!       This is  a time-dependent version of subroutine of k-loop
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
        subroutine tddiag_k ( isH )

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
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        logical, intent (in) :: isH

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

!$omp parallel do private (k_temp)
        do ikpoint = 1, nkpoints

! The subroutine kspace wants the k-vector in inverse angstrom units.
! NOT pi/alat units.
         k_temp(:) = special_k(:,ikpoint)
!  diagonalize S
         call diag_Sk (iqout, icluster, iwrteigen, ikpoint, k_temp, nkpoints)

! diagonalize H (optional)
         if (isH) then
          call diag_Hk (iqout, icluster, iwrteigen, ikpoint, k_temp, nkpoints)
          if (iwrteigen .eq. 1) then
           write (19,*) ' ------ the energy eigenvalues ----'
           write (19,100) (eigen_k(imu,ikpoint), imu = 1, norbitals_new)
           write (20,*) ikpoint,(eigen_k(imu,ikpoint), imu=1,norbitals_new)
!          write (20,101) ikpoint,(eigen_k(imu,ikpoint), imu=1,norbitals_new)
          else
! orthogonalize H
          end if
         end if ! isH

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
101     format (i4, 4f11.5)
!101     format (i4, <norbitals_new>f11.5)


        return
      end subroutine tddiag_k
