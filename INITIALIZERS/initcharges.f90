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

 
! initcharges.f90
! Program Description
! ===========================================================================
!       This routine initializes the charges for the simulation. If the
! old charge file exists, then this will be used to determine the atom
! charges for each shell. If not, then the neutral atom charges will be
! used - obtained from the info.dat file.
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
        subroutine initcharges (natoms, nspecies, itheory, ifixcharge, symbol)
        use dimensions
        use interactions
        use charges
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ifixcharge
        integer, intent (in) :: itheory
        integer, intent (in) :: natoms
        integer, intent (in) :: nspecies
 
        character (len = 2), intent (in), dimension(natoms) :: symbol
 
! Output
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer in1
        integer issh
        integer num_atoms
 
        character (len=40) fromwhichfile
 
        logical chargefile

        real, parameter ::  qparm = 100.0d0

! Allocate arrays
! ===========================================================================
        allocate (nelectron(natoms))
        allocate (Qin(nsh_max, natoms))
        allocate (Qinmixer(nsh_max*natoms))
        allocate (QLowdin_TOT (natoms))
        allocate (QMulliken_TOT (natoms))
        allocate (Qout(nsh_max, natoms))
        allocate (Qoutmixer(nsh_max*natoms))
        allocate (dq(nspecies))
        allocate (Q0_TOT(natoms))
        allocate (Qin_es(nsh_max, natoms))
        allocate (QLowdin_TOT_es (natoms))
        allocate (Qout_es(nsh_max, natoms))

! NPA
        allocate (qaux(nsh_max,nspec_max))
 
! Procedure
! ===========================================================================
! Qneutral_total
      do iatom = 1, natoms
         Q0_TOT(iatom) = 0
         in1 = imass(iatom)
         do issh = 1, nssh(in1)
         Q0_TOT(iatom) = Q0_TOT(iatom) + Qneutral(issh,in1)
         end do
        end do


! By default the input charges are initialized to the neutral atom charges
        Qin = 0.0d0
 
        do iatom = 1, natoms
         in1 = imass(iatom)
         do issh = 1, nssh(in1)
          Qin(issh,iatom) = Qneutral(issh,in1)
         end do
        end do
!NPA 
! setup aux array for Natual Population Analysis
        qaux = qparm
         do in1 = 1,nspecies
          do issh = 1, nssh(in1)
!            if (Qneutral(issh,in1) .gt. 0.0d0) qaux(issh,in1) = 1.0d0
          enddo
         enddo

! If there is a charge file from a previous run, or the user added one,
! then initialize the input charges accordingly.
        inquire (file = 'CHARGES', exist = chargefile)
        if (ifixcharge .eq. 1 .and. .not. chargefile) then
         write (*,*) ' ifixcharge = 1, but there is no CHARGES file! '
         stop
        end if
        if (chargefile .and. itheory .ne. 0) then
         write (*,*) '  '
         write (*,*) ' We are reading from a charge file. '
         open (unit = 12, file = 'CHARGES', status = 'old')
         read (12,100) num_atoms, fromwhichfile
         write (*,101) fromwhichfile
         if (num_atoms .ne. natoms) then
          write (*,*) ' The charge file that you are using must not belong '
          write (*,*) ' to the basis file that you are now calculating. '
          write (*,*) ' The number of atoms differs between the two. '
          write (*,*) ' Check the CHARGES file and start over! '
          stop
         end if
         do iatom = 1, natoms
          in1 = imass(iatom)
          read (12,*) (Qin(issh,iatom), issh = 1, nssh(in1))
         end do
         close (unit = 12)
        end if
 
! Write out the input charges.
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' The Input Charges Are: '
        write (*,200)
        write (*,201)
        write (*,200)
        do iatom = 1, natoms
         in1 = imass(iatom)
         write (*,202) iatom, symbol(iatom), nssh(in1),   &
     &                 (Qin(issh,iatom), issh = 1, nssh(in1))
        end do
        write (*,200)
 
! Format Statements
! ===========================================================================
100     format (2x, i5, 2x, a40, 2x, i2)
101     format (2x, ' Charge file corresponds to basisfile = ', a40)
200     format (2x, 70('='))
201     format (2x, ' Atom # ', 2x, ' Type ', 2x, ' Shells ', 1x,' Charges ')
202     format (3x, i5, 7x, a2, 5x, i2, 4x, 8(1x, f5.2))
 
 
        return
        end
