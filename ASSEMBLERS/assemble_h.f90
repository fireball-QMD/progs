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
!       This routine reads the different data file.
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine assemble_h ()

        use options
        use interactions
        use scf
        use hartree_fock

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Output

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================

       integer iatom, jatom, iorb, jorb
       logical readnij



! Procedure
! ===========================================================================

!        write (*,*) '  '
!        write (*,100)
!        write (*,*) ' Now we are assembling piecies of Hamiltonian. '
!        write (*,*) '  '

! read extended Hubbard
        if (itheory .eq. 2) then
         call assemble_eh ()
!         return
!        endif

! read Horsfield data
        elseif (itheory_xc .eq. 0) then
          call assemble_hxc ()
!          return
!        endif

! read McWeda data
        elseif (itheory_xc .ne. 0) then
          call assemble_mcweda ()
!          return
        endif

! GAP ENRIQUE-FF
        if ((igap.eq.1) .or. (igap.eq.2)) then
! read the previous nij if it exists, if not, consider nij=0
	     inquire (file = './nijmatrix', exist = readnij)
	     if (readnij) then
          open(271, file = './nijmatrix', status='old')
	      write (*,*) 'nijmatrix exists'
          do iatom = natomhf_beg, natomhf_end
            do jatom = natomhf_beg, natomhf_end
              do iorb = 1, numorb_max
                do jorb = 1, numorb_max
                  read(271,*) nij(iorb,iatom,jorb,jatom)
                end do
              end do
            end do
          end do
	      close(271)
         else
	       nij = 0.0
	     end if
        end if

        if (igap.eq.1) then
	     call assemble_hartree()
        end if

        if ((igap.eq.3).and.(Kscf.eq.1)) then
         call assemble_scissor()
        end if
! end GAP ENRIQUE-FF



! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))


        return
        end subroutine assemble_h

