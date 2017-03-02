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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

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

        write (*,*) '  '
        write (*,100)
        write (*,*) ' Now we are assembling piecies of Hamiltonian. '
        write (*,*) '  '

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

