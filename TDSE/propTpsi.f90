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


! propT_psi.f90
! Program Description
! ===========================================================================
!       This routine subroutine propagates wave-function in time base on
! the Crank-Nicholson scheme.
! Not: both psi and Hamiltonian are in MO(Lowdin) representation
!
! we propagate the wf solving a linear equation  L*psi(d+dt)=b
! where:
! L = I + img*(dt/2)*H
! b = (I - img*(dt/2)*H)*psi(t)
!
! to solve the linear equation we use LU-decomposition
!
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
        subroutine propTpsi ()

        use tdse
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
       real, parameter :: hbar = 0.65822d0


! Local Variable Declaration and Description
! ===========================================================================
       integer  ielec
       integer  imu
       integer  jmu
       integer  inu
       integer  jnu
       integer  ikpoint
       integer  iatom
       integer  in1
       integer  in2
       integer  ineigh
       integer  jatom
       integer  mbeta
       integer  info
       integer, dimension (norbitals) :: ipivot


       real fact
       real dot
       real norm
       real nfact
       real rval
       real ival
       real, dimension (3) :: vec

       complex, dimension (norbitals) :: psix
       complex, dimension (norbitals) :: b
       complex, dimension (norbitals,norbitals) :: mat
       complex, dimension (norbitals,norbitals) :: work
       complex, dimension (norbitals,norbitals) :: Lmat_save
       complex, dimension (norbitals,norbitals) :: Lmat
       complex phase
       complex mphase
       complex ai
       complex a0
       complex a1
       complex cnorm

! Procedure
! ===========================================================================
! NOTE: s12psi has been calculated in the subroutine tddenmat for
! former ionic configuration
        ai = cmplx(0.0d0, -1.0d0)
        a0 = cmplx(0.0d0, 0.0d0)
        a1 = cmplx(1.0d0, 0.0d0)

! ide matrix
        Lmat = a0
        work = a0
        do inu = 1, norbitals
         Lmat(inu,inu) = a1
         work(inu,inu) = a1
        enddo ! inu

        fact = 0.5d0*dte/hbar
        phase = fact*ai
        mphase = -1.0d0*fact*ai

! Loop over kpoints
        do ikpoint = 1, nkpoints

! restore orthogonal (MO) Hamiltonian
         do imu = 1, norbitals
          do inu = 1, norbitals
           mat(inu,imu) = HLow(inu,imu,ikpoint)
          enddo
         enddo

! calculate L = I + const*H
! input:
! work = I
! Lmat = I
! mat = H

         call zhemm ( 'L', 'U', norbitals, norbitals, phase, work, &
     &               norbitals, mat, norbitals, a1, Lmat, norbitals )
! output:
! Lmat = L
         Lmat_save = Lmat

! Note: for large system we probably need iterative solvers??
! factorize general matrix Lmat (LU decomposition)
         call zgetrf (norbitals, norbitals, Lmat, norbitals, ipivot, info)
         if (info .ne. 0) then
          write (*,*) '  '
          write (*,*) ' Factorization not successful, info = ', info
          if (info .lt. 0) then
           write (*,*) ' The ', info, '-th argument had an illegal '
           write (*,*) ' value. '
          endif
          stop
         endif

! calculate b = (I - const*HLow)*psi
! step 1.
!  input:
!   work = I
!   mat = H

         call zhemm ( 'L', 'U', norbitals, norbitals, mphase, work, &
     &               norbitals, mat, norbitals, a1, work, norbitals )
! output:
! work = (I - const*H)


! loop over electrons
         do ielec = 1, nelec
! copy psi of an electron to temp array
          do imu = 1, norbitals
           psix(imu) = psi(imu,ielec,ikpoint)
          enddo ! imu

! INPUT
! work = (I - const*H)
! psix = psi(to)
! work no longer HERMITIAN
! get b = (I - const*H)*psi for given electron
!          call zhemm ('L', 'U', norbitals, 1, a1, work, &
!     &               norbitals, psix, norbitals, a0, b, norbitals)
          call zgemm ('N', 'N', norbitals, 1, norbitals, a1, work,&
     &               norbitals, psix, norbitals, a0, b, norbitals)

! solve linera equation L*x=b (using LU-decomposition) for general matrix
          call zgetrs ('N', norbitals, 1, Lmat, norbitals, ipivot, b, &
     &                norbitals, info)
          if (info .ne. 0) then
           write (*,*) '  '
           write (*,*) ' LU solution not successful, info = ', info
           if (info .lt. 0) then
            write (*,*) ' The ', info, '-th argument had an illegal '
            write (*,*) ' value. '
           endif
           stop
          endif

! save new wavefunction
          do imu = 1, norbitals
           psi(imu,ielec,ikpoint) = b(imu)
          enddo ! imu

! check normalization of wf
!          norm = 0.0d0
!          do imu = 1, norbitals
!           norm = norm + real(psi(imu,ielec,ikpoint))**2
!          norm = norm + imag(psi(imu,ielec,ikpoint))**2
!         end do ! imu
!          write (*,500) ielec,norm
!          nfact = 1.0d0/sqrt(norm)
!          norm = 0.0d0
!          do imu = 1, norbitals
!           rval = real(psi(imu,ielec,ikpoint))*nfact
!           ival = imag(psi(imu,ielec,ikpoint))*nfact
!           psi(imu,ielec,ikpoint) = cmplx(rval,ival)
!           norm = norm + real(psi(imu,ielec,ikpoint))**2
!          norm = norm + imag(psi(imu,ielec,ikpoint))**2
!          enddo
!          write (*,501) ielec,norm
!           write (*,502) psi(:,ielec,ikpoint)
         end do ! ielec
        end do ! ikpoint


! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))
200     format (8f14.6)
500     format (3x, ' before: electron no.',i4,' norm  =',f14.8)
501     format (3x, ' after: electron no.',i4,' norm  =',f14.8)
502     format ('wfe: ',4f14.8)


        return
        end subroutine propTpsi

