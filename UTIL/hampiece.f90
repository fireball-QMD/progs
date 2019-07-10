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


! hampiece.f90
! Program Description
! ===========================================================================
!       This routine writes out pieces of the hamiltonian.
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
        subroutine hampiece (itheory)
        use configuration
        use dimensions
        use interactions
        use neighbor_map
        use density
        use options, only: V_intra_dip
        implicit none

! Argument Declaration and Description
! ===========================================================================
        integer, intent (in) :: itheory

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer jatom
        integer in1
        integer in2
        integer ineigh
        integer imu
        integer inu
        integer issh
        integer mbeta

        real distance

        real, dimension (3) :: dvec

! Procedure
! ===========================================================================
        write (*,*) ' '
        write (*,*) ' In hampiece.f90 '
        write (*,*) ' Writing out pieces of the Hamiltonian. '
        write (*,*)

! ===========================================================================
        do iatom = 1, natoms
         in1 = imass(iatom)
         write (*,*) '  '
         write (*,100)
         write (*,100)
         write (*,*) ' There are ', neighn(iatom), ' neighbors to atom ', iatom

         do ineigh = 1, neighn(iatom)
          write (*,*) '  '
          write (*,*) ' Matrices connecting to neighbor ', ineigh
          jatom = neigh_j(ineigh,iatom)
          mbeta = neigh_b(ineigh,iatom)
          in2 = imass(jatom)

          dvec = ratom(:,jatom) + xl(:,mbeta) - ratom(:,iatom)
          distance = sqrt(dvec(1)**2 + dvec(2)**2 + dvec(3)**2)
          write (*,200) jatom, mbeta, dvec, distance

          write (*,*) ' Interactions: s, t, vna, vxc, vxc_1c, vnl '
          write (*,*) '  '
          write (*,*) ' s: overlap '
          write (*,300)
          do imu = 1, num_orb(in1)
           write(*,301) (s_mat(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
          end do
          write (*,*) '  '
          write (*,*) ' t: kinetic '
          write (*,300)
          do imu = 1, num_orb(in1)
           write(*,301) (t_mat(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
          end do
          write (*,*) '  '
          write (*,*) ' vna: neutral atom '
          write (*,300)
          do imu = 1, num_orb(in1)
           write(*,301) (vna(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
          end do
          write (*,*) '  '
          write (*,*) ' vxc: exchange/correlation'
          write (*,300)
          do imu = 1, num_orb(in1)
           write(*,301) (vxc(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
          end do
          if (iatom .eq. jatom .and. mbeta .eq. 0 .and. itheory .ne. 3) then
           write (*,*) '  '
           write (*,*) ' vxc_1c: one-center exchange/correlation'
           write (*,300)
           do imu = 1, num_orb(in1)
            write(*,301) (vxc_1c(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
           end do
          end if

          if (itheory .eq. 1) then
           write (*,*) '  '
           write (*,*) '  '
           write (*,*) ' DOGS Interactions: vca, vxc_ca, ewaldlr, and ewaldsr '
           write (*,*) '  '
           write (*,*) ' vca: charged atom '
           write (*,300)
           do imu = 1, num_orb(in1)
            write(*,301) (vca(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
           end do
           write (*,*) '  '
           write (*,*) ' vxc_ca: charged exchange/correlation'
           write (*,300)
           do imu = 1, num_orb(in1)
            write(*,301) (vxc_ca(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
           end do
           write (*,*) '  '
           write (*,*) ' ewaldlr: long-range ewald '
           write (*,300)
           do imu = 1, num_orb(in1)
            write(*,301) (ewaldlr(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
           end do
           write (*,*) '  '
           write (*,*) ' ewaldsr: short-range ewald '
           write (*,300)
           do imu = 1, num_orb(in1)
            write(*,301) (ewaldsr(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
           end do
          end if

          if (itheory .eq. 2) then
           write (*,*) '  '
           write (*,*) '  '
           write (*,*)                                                      &
     &      ' Extended-Hubbard Interactions: vca, vxc_ca, ewaldlr, and ewaldsr '
           write (*,*) ' vca: charged atom '
           write (*,300)
           do imu = 1, num_orb(in1)
            write(*,301) (vca(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
           end do
           write (*,*) '  '
           write (*,*) ' vxc_ca: charged exchange/correlation'
           write (*,300)
           do imu = 1, num_orb(in1)
            write(*,301) (vxc_ca(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
           end do
           write (*,*) '  '
           write (*,*) ' ewaldlr: long-range ewald '
           write (*,300)
           do imu = 1, num_orb(in1)
            write(*,301) (ewaldlr(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
           end do
           write (*,*) '  '
           write (*,*) ' ewaldsr: short-range ewald '
           write (*,300)
           do imu = 1, num_orb(in1)
            write(*,301) (ewaldsr(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
           end do
          end if

! KS theory
          if (itheory .eq. 3) then
           write (*,*) '  '
           write (*,*) '  '
           write (*,*) ' KS interaction: dVh '
           write (*,*) '  '
           write (*,*) ' vca: charged atom '
           write (*,300)
           do imu = 1, num_orb(in1)
            write(*,301) (vca(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
           end do
           write (*,*) ' rho: from grid '
           write (*,300)
           do imu = 1, num_orb(in1)
            write(*,301) (vxc_ca(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
           end do
           write (*,*) '  '
           write (*,*) ' vna: from grid'
           write (*,300)
           do imu = 1, num_orb(in1)
            write(*,301) (vxc_1c(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
           end do
          endif

          write (*,*) '  '
          write (*,*) ' Total Hamiltonian: h '
          write (*,300)
          do imu = 1, num_orb(in1)
           write(*,301) (h_mat(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
          end do
         end do
        end do
        write (*,100)
        write (*,100)


! PP matrix
        do iatom = 1, natoms
         in1 = imass(iatom)
         write (*,*) '  '
         write (*,*) '       VNL matrices    '
         write (*,*) '  '
         write (*,100)
         write (*,100)
         write (*,*) ' There are ',neighPPn(iatom),' neighbors to atom ',iatom

         do ineigh = 1, neighPPn(iatom)
          write (*,*) '  '
          write (*,*) ' Matrices connecting to neighbor ', ineigh
          jatom = neighPP_j(ineigh,iatom)
          mbeta = neighPP_b(ineigh,iatom)
          in2 = imass(jatom)

          dvec = ratom(:,jatom) + xl(:,mbeta) - ratom(:,iatom)
          distance = sqrt(dvec(1)**2 + dvec(2)**2 + dvec(3)**2)
          write (*,200) jatom, mbeta, dvec, distance


          write (*,*) '  '
          write (*,*) ' vnl: non-local'
          write (*,300)
          do imu = 1, num_orb(in1)
           write(*,301) (vnl(imu,inu,ineigh,iatom), inu = 1, num_orb(in2))
          end do

         end do
        end do
        write (*,100)
        write (*,100)
! V_intra_dip
          if (V_intra_dip .eq. 1) then
          write (*,*) '  '
          write (*,*) ' V_intra_dip_1c'
          write (*,300)
          do iatom = 1, natoms
          do imu = 1, num_orb(in1)
           write(*,301) (Vdip_1c(imu,inu,iatom), inu = 1, num_orb(in2))
          end do !end do imu
          end do !end do iatom
          end if !end if V_intra_dip .eq. 1
          
! Format Statements
! ===========================================================================
100     format (75('*'))
200     format (2x, ' jatom, mbeta = ', 2i3, ' r = ', 3f9.4,                &
     &              ' distance = ', f10.6)
300     format (75('='))
301     format (9f9.4)
302     format (3i4, 3f7.3, f10.4)

        return
        end
