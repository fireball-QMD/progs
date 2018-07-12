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


! readinfo.f90
! Program Description
! ===========================================================================
!       This routine reads the info.dat file.
!
! ===========================================================================
! Code originally written by Juergen Fritsch and Otto F. Sankey
 
! Code rewritten by:
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
        subroutine readinfo ()
        use charges
        use dimensions
        use interactions
        use configuration 
        use integrals
        use options, only : verbosity, inputxyz

        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Output

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ispec
        integer jspec
        integer iatom
        integer iline
        integer issh
        integer imu
        integer ins
        integer temp_nsup
        integer nzx_temp
 
        integer, dimension (50) :: temp_nsu
        character (len=2), dimension (50) :: temp_symbol
 
        real, dimension (:,:), allocatable :: cutoff  ! cutoff radius in bohr
        real :: nucz
 
        character (len=25) :: signature
        character (len=2)  :: symbolA_temp
 
! jel-grid
        character (len=1)  :: charx
        integer nchar
!        character (len=25), dimension (:,:), allocatable :: wavefxn
!        character (len=25), dimension (:,:), allocatable :: napot
! end jel-grid
 
        logical skip_it
        logical new_one
! zdenka chromcova - test that info.dat is ok
		integer foundspecies
 		
 		foundspecies=0
! Procedure
! ===========================================================================
! Initialize rcutoff array - we loop over all of it elsewhere
		
        write (*,*) '  '
!       write (*,100)
        write (*,*) ' Now we are reading from the info.dat file. '
        write (*,*) '  '

 
! Open the data file.
        open (unit = 12, file = trim(fdataLocation)//'/info.dat', status = 'old')
 
        if (verbosity .ge. 3)  write (*,*) ' You are using the database created by: '
        read (12,101) signature
        if (verbosity .ge. 3) write (*,101) signature
        read (12,*) nspecies
        if (verbosity .ge. 3) write (*,*) '  '
        if (verbosity .ge. 3) write (*,*) ' Number of species in database = ', nspecies

! Allocate nzx /
        allocate (nzx (nspecies))
        allocate (symbolA (nspecies)) 
        allocate (etotatom (nspecies)) 
        allocate (smass (nspecies)) 
        allocate (rc_PP (nspecies)) 
        allocate (rcutoff (nspecies, nsh_max)) 
        rcutoff = 0.0d0



        if (inputxyz .eq. 1) then

          open  (unit = 69, file = basisfile, status = 'old')
          read (69, *) natoms
          read(69,*)
! Loop over the number of atoms
          temp_nsup = 0
          do iatom = 1, natoms
           read (69,*) symbolA_temp
           new_one = .true.
           do imu = 1, temp_nsup
            if (trim(symbolA_temp) .eq. temp_symbol(imu)) new_one = .false.
           end do
           if (new_one) then
            temp_nsup = temp_nsup + 1
            temp_symbol(temp_nsup) = trim(symbolA_temp)
           end if
          end do
          close (unit = 69)

        else
! It's handy to not read all the files if you only want some of them
!        open (unit = 41, file = 'script.input', status = 'old')
!        read (41,*) basisfile
!        close (unit = 41)
          open  (unit = 69, file = basisfile, status = 'old')
          read (69, *) natoms

! Loop over the number of atoms
          temp_nsup = 0
          do iatom = 1, natoms
           read (69,*) nucz
           new_one = .true.
           do imu = 1, temp_nsup
            if (nucz .eq. temp_nsu(imu)) new_one = .false.
           end do
           if (new_one) then
            temp_nsup = temp_nsup + 1
            temp_nsu(temp_nsup) = nucz
           end if
          end do
          close (unit = 69)
        end if
        if (verbosity .ge. 3) write (*,*) ' We will only read these atomic indexes: '
        if (verbosity .ge. 3) write (*,*) temp_nsu(1:temp_nsup)
        if (verbosity .ge. 3) write (*,*)

        rewind (unit = 12)
        read (12,*)
        read (12,*)
 
! Make sure that nspec is not greater than nspec_max
        if (temp_nsup .gt. nspec_max) then
         write (*,*) ' nspecies = ', temp_nsup,' nspec_max = ', nspec_max
         write (*,*) ' Sorry -- redimension nspec_max in MODULES/dimensions.f90'
         stop
        end if

! allocate interactions
        allocate (cl_PP (0:nsh_max - 1, temp_nsup))
        allocate (nssh (temp_nsup))
        allocate (lssh (nsh_max, temp_nsup))
        allocate (nsshPP (temp_nsup))
        allocate (lsshPP (nsh_max, temp_nsup))
        allocate (Qneutral (nsh_max, temp_nsup))

! allocate local arrays
        allocate (cutoff (nsh_max, temp_nsup))
        allocate (wavefxn (nsh_max, temp_nsup))
        allocate (napot (0:nsh_max, temp_nsup))

! Now read in the data
        ispec = 1
        nsup = 0
        do jspec = 1, nspecies
         read (12,*)
         read (12,*)
         read (12,102) symbolA_temp
         read (12,*) nzx_temp
         if (verbosity .ge. 3) write(*,102) symbolA_temp
         if (verbosity .ge. 3) write (*,*) nzx_temp
         skip_it = .true.
         do ins = 1, temp_nsup
          if (inputxyz .eq. 0 .and. temp_nsu(ins) .eq. nzx_temp) skip_it = .false.
          if (inputxyz .eq. 1 .and. temp_symbol(ins) .eq. trim(symbolA_temp)) skip_it = .false.
         end do

! Just skip this atom
         if (skip_it) then
 
! nsu and nsup list the atoms the are skipped, by position in
! the info.dat file.  This is the format that the rest of fireball
! uses, since it is simpler to use in read1c, but we read in the atoms
! we want, since they are nicer to input.
! This means that the code in this subroutine is not nice to look at.
          nsup = nsup + 1
          nsu(nsup) = jspec
          do iline = 1, 12
           read (12,*)
          end do
         else
          foundspecies=foundspecies+1
          symbolA(ispec) = symbolA_temp
          nzx(ispec) = nzx_temp

          read (12,*) smass(ispec)
          read (12,*) nssh(ispec)

! Check and make sure that the number of shells is not greater than
! the maximum - nsh_max
          if (nssh(ispec) .gt. nsh_max) then
           write (*,*) ' nssh(ispec) = ', nssh(ispec),' nsh_max = ', nsh_max
           write (*,*) ' Sorry -- redimension nsh_max in MODULES/dimensions.f90'
           stop
          end if
          if (nssh(ispec) .gt. 8) then
           write (*,*) ' nssh(ispec) = ', nssh(ispec)
           write (*,*) ' Sorry -- Currently the basis set cannot be larger '
           write (*,*) ' than s, s*, p, p*, d, d*, f, and f* '
           stop
          end if
 
          read (12,*) (lssh(issh,ispec), issh = 1, nssh(ispec))
          read (12,*) nsshPP(ispec)
          read (12,*) (lsshPP(issh,ispec), issh = 1, nsshPP(ispec))
          read (12,*) rc_PP(ispec)
          read (12,*) (Qneutral(issh,ispec), issh = 1, nssh(ispec))

          read (12,*) (cutoff(issh,ispec), issh = 1, nssh(ispec))
          do issh = 1, nssh(ispec)
           rcutoff(ispec, issh) = cutoff(issh,ispec)*0.529177d0
          end do
 
          read (12,103) (wavefxn(issh,ispec), issh = 1, nssh(ispec))
          read (12,103) (napot(issh,ispec), issh = 0, nssh(ispec))
! jel-grid
! adjust the potential file name
          do issh = 0, nssh(ispec)
           nchar = len_trim(napot(issh,ispec))
           do imu = 1,nchar
            charx = napot(issh,ispec)(imu:imu)
            if ( lle(charx,'/') .and. lge(charx,'/') ) then 
             ins = imu
            endif
           enddo ! do imu
	   napot(issh,ispec) = trim('/basis/'//napot(issh,ispec)(ins+1:nchar))
          enddo ! do issh
! adjust the wavefunction file name
          do issh = 1, nssh(ispec)
           nchar = len_trim(wavefxn(issh,ispec))
           do imu = 1,nchar
            charx = wavefxn(issh,ispec)(imu:imu)
            if ( lle(charx,'/') .and. lge(charx,'/') ) then 
             ins = imu
            endif
           enddo ! do imu
	   wavefxn(issh,ispec) = trim('/basis/'//wavefxn(issh,ispec)(ins+1:nchar))
          enddo ! do issh 
! end jel-grid
          read (12,*) etotatom(ispec)
          read (12,*)

! Write out.
          if (verbosity .ge. 3)  then
            write (*,100)
            write (*,301) ispec
            write (*,302) symbolA(ispec)
            write (*,303) nzx(ispec)
            write (*,304) smass(ispec)
            write (*,305) nssh(ispec)
            write (*,306) (lssh(issh,ispec), issh = 1, nssh(ispec))
            write (*,307) nsshPP(ispec)
            write (*,308) (lsshPP(issh,ispec), issh = 1, nsshPP(ispec))
            write (*,314) rc_PP(ispec)
            write (*,309) (Qneutral(issh,ispec), issh = 1, nssh(ispec))
            write (*,310) (cutoff(issh,ispec), issh = 1, nssh(ispec))
            write (*,311) (wavefxn(issh,ispec), issh = 1, nssh(ispec))
            write (*,312) (napot(issh,ispec), issh = 0, nssh(ispec))
            write (*,313) etotatom(ispec)
            write (*,100)
          endif !verbosity 
! Increment ispec, since we read in data
          ispec = ispec + 1
         end if
        end do
        
        nspecies = temp_nsup

        if(foundspecies /= nspecies) then
        	write(*,*)'In Fdata/info.dat is defined only ',foundspecies,' species from ',nspecies
        	write(*,*)'Exiting. Please check the info.dat'
			stop
		end if

! Deallocate Arrays
! ===========================================================================
! jel-grid
        deallocate (cutoff)
! end jel-grid
 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
101     format (2x, a25)
102     format (2x, a2)
103     format (9(2x,a25))
301     format (2x, i2, ' - Information for this species ')
302     format (2x, a2, ' - Element ')
303     format (2x, i3, ' - Nuclear Z ')
304     format (2x, f7.3, ' - Atomic Mass ')
305     format (2x, i2, ' - Number of shells ')
306     format (2x, 8(2x,i1), ' - L; quantum number for each shell ')
307     format (2x, i2, ' - Number of shells (Pseudopotential) ')
308     format (2x, 8(2x,i1), ' - L; quantum number for each shell ')
309     format (2x, 8(2x,f5.2), ' - Occupation numbers ')
310     format (2x, 8(2x,f5.2), ' - Radial cutoffs ')
311     format (2x, 9(2x,a25), ' - Wavefunction files ')
312     format (2x, 9(2x,a25), ' - (Non)-neutral atom potentials ')
313     format (2x, f12.4, ' - Atomic energy ')
314     format (2x, f12.4, ' - Radial cutoffs (Pseudopotential) ')
        close (unit = 12)

        return
        end subroutine readinfo
 
