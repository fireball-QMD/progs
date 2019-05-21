! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
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

 
! fermie.f90
! Program Description
! ===========================================================================
!       This routine calculates the fermi energy.
!
! ===========================================================================
! Code written by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine fermie (norbitals, ztot, eigen_k, efermi, &
      &                    ioccupy_k, foccupy)
        use dimensions
        use constants_fireball
        use kpoints
        use scf 

        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: norbitals
 
        real, intent(in) :: ztot
        real, intent(inout), dimension (norbitals, nkpoints) :: eigen_k
 
! Output
        integer, intent(out), dimension (norbitals, nkpoints) :: ioccupy_k
        real, intent(out) :: efermi
        real, intent(out), dimension (norbitals, nkpoints) :: foccupy
 
! Local Parameters
! ===========================================================================
        integer, parameter :: imax = 1000 ! maximum sc iterations
        integer, parameter :: nmax = 5000 ! cutoff for degeneracy check
 
        real, parameter :: tol = 1.0d-10
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ikpoint
        integer imu
        integer inu
        integer iter
        integer jkpoint
 
        real delta
        real emin
        real emax
        real qcharge
        real qztot
        real temp
 
! Procedure
! ===========================================================================
! Add in the qstate to the total charge
! HAO--possible fix

         qztot = ztot 

! Initialize the occupation numbers and foccupy to zero.
        ioccupy_k = 0
        foccupy = 0.0d0
   
! The subroutine fermie needs a temperature to calculate the occupations of
! the states so set temperature to some low value (1 eV = 11604 K).
        temp = tempfe/kconvert
 
! Find emin and emax. Also make sure degenerate eigenvalues are truly
! degenerate. However, if norbitals*nkpts*norbitals*nkpts is larger than nmax,
! then skip the degeneracy checking. Otherwise, the checking can take a while.
        if (norbitals**2*nkpoints**2 .lt. nmax) then
         emin = eigen_k(1,1)
         emax = eigen_k(norbitals,1)
         do ikpoint = 1, nkpoints
          do imu = 1, norbitals
           if (eigen_k(imu,ikpoint) .lt. emin) emin = eigen_k(imu,ikpoint)
           if (eigen_k(imu,ikpoint) .gt. emax) emax = eigen_k(imu,ikpoint)
           do jkpoint = ikpoint, nkpoints
            do inu = imu, norbitals
             if (abs(eigen_k(imu,ikpoint) - eigen_k(inu,jkpoint)) .lt. tol) then
              eigen_k(inu,jkpoint) =                                         &
     &         (eigen_k(imu,ikpoint) + eigen_k(inu,jkpoint))/2.0d0
              eigen_k(imu,ikpoint) = eigen_k(inu,jkpoint)
             end if
            end do
           end do
          end do
         end do
        else
!         write (*,*) '  '
!         write (*,*) ' ************ WARNING ******** WARNING ************* '
!         write (*,*) '          skipping the degeneracy checking  '
!         write (*,*) '               in subroutine fermie'
!         write (*,*) ' *************************************************** '
         emin = eigen_k(1,1)
         emax = eigen_k(norbitals,1)
         do ikpoint = 1, nkpoints
          do imu = 1, norbitals
           if (eigen_k(imu,ikpoint) .lt. emin) emin = eigen_k(imu,ikpoint)
           if (eigen_k(imu,ikpoint) .gt. emax) emax = eigen_k(imu,ikpoint)
          end do
         end do
        end if

        emax=emax+0.20d0 
! The value of efermi must be between emin and emax
        iter = 0
        qcharge = 0.0d0
        do while (abs(qcharge - qztot) .gt. tol .and. iter .le. imax)
         iter = iter + 1
! Make a guess at efermi
         efermi = (emax + emin)/2.0d0
         qcharge = 0.0d0
         do ikpoint = 1, nkpoints
          do imu = 1, norbitals
           delta = (eigen_k(imu,ikpoint) - efermi)/temp
! Skip exponential for big -/+ delta
           if (delta .gt. 10.0d0) then
            foccupy(imu,ikpoint) = 0.0d0
            ioccupy_k(imu,ikpoint) = 0
           else if (delta .lt. -10.0d0) then
            foccupy(imu,ikpoint) = 1.0d0
            ioccupy_k(imu,ikpoint) = 1
           else
            foccupy(imu,ikpoint) = 1.0d0/(1.0d0 + exp(delta))
            if (foccupy(imu,ikpoint) .gt. 1.0d-5) then
             ioccupy_k(imu,ikpoint) = 1
            else
             ioccupy_k(imu,ikpoint) = 0
            end if
           end if
           qcharge = qcharge + spin*foccupy(imu,ikpoint)*weight_k(ikpoint)
          end do ! do imu
         end do ! do ikpoint
! Narrow the range that efermi can fall into
         if (qcharge .gt. qztot) then
          emax = efermi
         else
          emin = efermi
         end if
        end do

! cDFT vlada
! Create electron and hole
   !     if (cDFT_active) then
   !      do ikpoint = 1, nkpoints
   !       call project_eh (ioccupy_k, foccupy, ikpoint)
   !      end do
   !     end if 
! cDFT vlada
 
! Print warning for going over maximum iterations.
        if (iter .gt. imax) then
         write (*,*) '  '
         write (*,*) ' ************ WARNING ******** WARNING ************* '
         write (*,*) '        not under tolerance (toll) after ',imax
         write (*,*) '          iterations in subroutine fermie'
         write (*,*) '  '
         write (*,*) '          qcharge = ', qcharge
         write (*,*) '          qztot = ', qztot
         write (*,*) '          emax = ', emax
         write (*,*) '          emin = ', emin
         write (*,*) ' *************************************************** '
         write (*,*) '  '
        end if
 
! Format Statements
! ===========================================================================
 
        return
        end
