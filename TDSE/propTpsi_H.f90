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
! we propagate the wf solving a linear equation  L*psi(d+dt)=b
! where:
! L = I + img*(dt/2)*(S^-1)*H
! b = (I + img*(dt/2)*(S^-1)*H)*psi(t)
!
! to solve the linear equation we use Cholesky decomposition
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
       real, dimension (3) :: vec

       complex, dimension (norbitals) :: psix
       complex, dimension (norbitals) :: b
       complex, dimension (norbitals,norbitals) :: mat
       complex, dimension (norbitals,norbitals) :: work
       complex, dimension (norbitals,norbitals) :: swork
       complex, dimension (norbitals,norbitals) :: Lmat
       complex phase
       complex mphase
       complex ai
       complex a0
       complex a1


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

        fact = 0.5d0*hbar*dte
        phase = cmplx(cos(fact),sin(fact))
        fact = 0.5d0*hbar*dte*ai
        mphase = cmplx(cos(fact),-sin(fact))

        write (*,*) 'phase =',phase
        write (*,*) 'mphase =',mphase

        do ikpoint = 1, nkpoints

! We now loop over all neighbors jatom of iatom.
         do iatom = 1, natoms
          in1 = imass(iatom)

! Now loop over all neighbors jatom of iatom.
          do ineigh = 1, neighn(iatom)
           mbeta = neigh_b(ineigh,iatom)
           jatom = neigh_j(ineigh,iatom)
           in2 = imass(jatom)

! So this matrix element goes in the i,j slot.
           vec(:) = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
           dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2) &
    &            + special_k(3,ikpoint)*vec(3)
           phase = cmplx(cos(dot),sin(dot))
           do inu = 1, num_orb(in2)
            jnu = inu + degelec(jatom)
            do imu = 1, num_orb(in1)
             jmu = imu + degelec(iatom)
             mat(jmu,jnu) = mat(jmu,jnu) + phase*h_mat(imu,inu,ineigh,iatom)
            end do ! do inu
           end do ! do imu
          end do ! do ineigh

! Now loop over all neighbors jatom of iatom VNL
          do ineigh = 1, neighPPn(iatom)
           mbeta = neighPP_b(ineigh,iatom)
           jatom = neighPP_j(ineigh,iatom)
           in2 = imass(jatom)

! So this matrix element goes in the i,j slot.
           vec(:) = xl(:,mbeta) + ratom(:,jatom) - ratom(:,iatom)
           dot = special_k(1,ikpoint)*vec(1) + special_k(2,ikpoint)*vec(2) &
    &            + special_k(3,ikpoint)*vec(3)
           phase = cmplx(cos(dot),sin(dot))

           do inu = 1, num_orb(in2)
            jnu = inu + degelec(jatom)
            do imu = 1, num_orb(in1)
             jmu = imu + degelec(iatom)
             mat(jmu,jnu) = mat(jmu,jnu) + phase*vnl(imu,inu,ineigh,iatom)
            end do ! do imu
           end do ! do inu
          end do ! do inegh
         end do ! do iatom

! copy S-1 to local array
         swork(:,:) = sm1(:,:,ikpoint)

         write (*,*) 'H'
         do imu = 1, norbitals
          write (*,200) (real(mat(imu,inu)), inu=1,norbitals)
         enddo
         write (*,*)
         write (*,*) 'S-1'
         do imu = 1, norbitals
          write (*,200) (real(swork(imu,inu)), inu=1,norbitals)
         enddo
         write (*,*)
         write (*,*) 'I'
         do imu = 1, norbitals
          write (*,200) (real(Lmat(imu,inu)), inu=1,norbitals)
         enddo
         write (*,*)
! calculate L = I + const*S^(-1)*H
! input:
! Lmat = I
! swork = S^-1
! mat = H
         call zhemm ( 'L', 'U', norbitals, norbitals, phase, swork, &
     &               norbitals, mat, norbitals, a1, Lmat, norbitals )
! output:
! Lmat = L

         write (*,*) 'Lmat'
         do imu = 1, norbitals
          write (*,200) (real(Lmat(imu,inu)), inu=1,norbitals)
         enddo
         write (*,*)
! factorize matrix Lmat for later (Hermitian)
!         call zpotrf ('U', norbitals, Lmat, norbitals, info)
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

! calculate b = (I - const*S(^-1)*H)*psi
! step 1.
!  input:
!   work = I
!   sm1 = S^-1
!   mat = H
         call zhemm ( 'L', 'U', norbitals, norbitals, mphase, swork, &
     &               norbitals, mat, norbitals, a1, work, norbitals )
! output:
! work = (I - const*S(^-1)*H)

! loop over electrons
         do ielec = 1, nelec
! copy psi of an electron to temp array
          do imu = 1, norbitals
           psix(imu) = psi(imu,ielec,ikpoint)
          enddo ! imu
! get b for given electron
          call zhemm ( 'L', 'U', norbitals, 1, a1, work, &
     &               norbitals, psix, norbitals, a0, b, norbitals )

! solve linear equation (Cholesky)
!          call zpotrs ('U', norbitals, 1, Lmat, norbitals, b, &
!     &                norbitals, info)
! solve linera equation (LU) general matrix
          call zgetrs ('N', norbitals, 1, Lmat, norbitals, ipivot, b, &
     &                norbitals, info)
          if (info .ne. 0) then
           write (*,*) '  '
           write (*,*) ' Chlesky solution not successful, info = ', info
           if (info .lt. 0) then
            write (*,*) ' The ', info, '-th argument had an illegal '
            write (*,*) ' value. '
           endif
           stop
          endif

          do imu = 1, norbitals
           psi(imu,ielec,ikpoint) = b(imu)
          end do ! imu
         end do ! ielec
        end do ! ikpoint



! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))
200     format (8f12.6)


        return
        end subroutine propTpsi

