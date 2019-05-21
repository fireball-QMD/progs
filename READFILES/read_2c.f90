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


! read_2c.f90
! Program Description
! ===========================================================================
!       The data of a given two center integral are stored in a file of the
! form named root.nz1.nz2.dat or root_n.nz1.nz2.dat where root indicates the
! type of 2 center integral.  The numbers nz1, nz2 are the nuclear charges
! of the two atoms for which the two center table was created.  The index n
! numbers the different subtypes of a given interaction (e.g. neutral atom
! potential and the anglular contributions).  By calling the subroutine,
! information about the interaction type has to be provided.  The routine will
! generate the appropriate name (root) and will open the respective files.
!
! Here's how the interactions are defined:
!
!          interaction    subtype
!        -----------------------------------------------------------
!               1            0            overlap
!               2            0            vna   ontopl
!               2            isorp        vna   ontopl shell isorp
!               3            0            vna   ontopr
!               3            isorp        vna   ontopr shell isorp
!               4            0            vna   atom
!               4            isorp        vna_  atom  shell isorp
!               5            0            non-local
!               6            0, 4         xc ontop
!               7            0, 4         xc atom-atom
!               8            0, 4         xc correction
!               9            0            z-dipole
!               10           0            y-dipole
!               11           0            x-dipole
!               12           0            coulomb
!               13           0            kinetic
!               14           0            extended-hubbard
!
!       The different subtypes for the xc matrix elements contain the neutral
! atom case, and the plus/minus delta Q cases which are needed to determine
! derivatives with respect to the charges.
!                                ! created in this routine
!
! The 2-center interactions are stored in an array.  They are all stored in
! the same array (this may change).  Here's how the indexes are defined:
!
!      Type         interaction  subtypes     index
!    ----------------------------------------------------------------
!      overlap            1          0          1
!      vna on l           2          0..8       2..11   (2  + isorp)
!      vna on r           3          0..8       12..21
!      vna atm            4          0..8       22..31
!      vnl                5          0          32
!      vxc on             6          0..4       33..37
!      vxc atm            7          0..4       38..42
!      xccorr             8          0..4       43..47
!      z-dipole           9          0          48
!      y-dipole          10          0          49
!      x-dipole          11          0          50
!      coulomb           12          0          51
!      kinetic           13          0          52
!      extended-hubbard  14          0          53
!      density_ontopl    15          1..8       54..63
!      density_ontopr    16          1..8
!      density_atom      17          1..8 
!      dnuxc_ol          18          1..8       
!      dnuxc_or          19          1..8
!      sph dens ontopl   20          1..4
!      sph dens_ontopr   21          1..4
!      sph den_atom      22          1..4
!      sph overlap       23          0
!
! Why does vna on l (for example) go from 0 to 9.  Well, 0 means neutral atom,
! and 1 to 9 are for the shells.  9 shells!  Holy smoke.  That's a lot. But
! here is why its 9 shells.  We include s, p, d, f and s*, p*, d*, f*.
! --> But isn't that a waste of memory which makes fireball run more slowly?
!
! ===========================================================================
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
        subroutine read_2c (interaction, nspecies, itheory, ioff2c, nzx)
        use dimensions
        use integrals
        use interactions
        use constants_fireball
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: interaction
        integer, intent (in) :: ioff2c
        integer, intent (in) :: itheory
        integer, intent (in) :: nspecies
 
        integer, intent (in), dimension (nspecies) :: nzx
        
! Output is in module integrals
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer in1
        integer in2
        integer initype
        integer iounit
        integer isorp
        integer issh
        integer isub2c
        integer itype
        integer maxtype
        integer npseudo
        integer num_nonzero
        integer numz
        integer nzx1
        integer nzx2
 
        real rc1
        real rc2
        real zmax
        real zmin
 
! We include up to 5 non-local (L values) of the pseudopotential. Usually,
! you will have 2 (L=0, 1 and sometimes 2 (D)).
        real, dimension (nsh_max) :: cl_pseudo
 
        character (len = 200) append_string
        character (len = 200) extension
        character (len = 200) filename
        character (len = 200) root
        character (len = 200) root_isorp
 
        external append_string

 
! Procedure
! ===========================================================================
! Set up isub2c. 
        isub2c = 0
        if (interaction .eq. 6) isub2c = 4
        if (interaction .eq. 7) isub2c = 4
        if (interaction .eq. 8) isub2c = 4
 
! Allocate and initialize arrays
        if(interaction .eq. 1) then
         allocate (xintegral_2c (ME2c_max, nfofx, interactions2c_max,        &
     &                           nspecies, nspecies))
         allocate (z2cmax (interactions2c_max, nspecies, nspecies))
         allocate (numz2c (interactions2c_max, nspecies, nspecies))
         xintegral_2c = 0.0d0
         if (superspline) then
          allocate (splineint_2c(1:4, ME2c_max, nfofx, interactions2c_max,   &
     &                           nspecies, nspecies))
          splineint_2c = 0.0d0
         end if
        end if
        iounit = 71

! Here are the roots for the file names for all interactions
        if (interaction .eq. 1)  root = trim(fdataLocation)//'/overlap'
        if (interaction .eq. 2)  root = trim(fdataLocation)//'/vna_ontopl'
        if (interaction .eq. 3)  root = trim(fdataLocation)//'/vna_ontopr'
        if (interaction .eq. 4)  root = trim(fdataLocation)//'/vna_atom  '
        if (interaction .eq. 5)  root = trim(fdataLocation)//'/vnl       '
        if (interaction .eq. 6)  root = trim(fdataLocation)//'/xc_ontop'
        if (interaction .eq. 7)  root = trim(fdataLocation)//'/xc_atom   '
        if (interaction .eq. 8)  root = trim(fdataLocation)//'/xc_corr   '
        if (interaction .eq. 9)  root = trim(fdataLocation)//'/dipole_z  '
        if (interaction .eq. 10) root = trim(fdataLocation)//'/dipole_y  '
        if (interaction .eq. 11) root = trim(fdataLocation)//'/dipole_x  '
        if (interaction .eq. 12) root = trim(fdataLocation)//'/coulomb   '
        if (interaction .eq. 13) root = trim(fdataLocation)//'/kinetic   '
        if (interaction .eq. 14) root = trim(fdataLocation)//'/nuxc      '
        if (interaction .eq. 15) root = trim(fdataLocation)//'/den_ontopl'
        if (interaction .eq. 16) root = trim(fdataLocation)//'/den_ontopr'
        if (interaction .eq. 17) root = trim(fdataLocation)//'/den_atom  ' 
        if (interaction .eq. 18) root = trim(fdataLocation)//'/dnuxc_ol  '
        if (interaction .eq. 19) root = trim(fdataLocation)//'/dnuxc_or  '
        if (interaction .eq. 20) root = trim(fdataLocation)//'/denS_ontopl'
        if (interaction .eq. 21) root = trim(fdataLocation)//'/denS_ontopr'
        if (interaction .eq. 22) root = trim(fdataLocation)//'/denS_atom  ' 
        if (interaction .eq. 23) root = trim(fdataLocation)//'/overlapS  ' 

! Now generate the file name of the file to be opened.  Loop over all cases of
! the interaction (e.g. different charges of the xc-stuff)




! Loop over atoms in1
        do in1 = 1, nspecies
 
! Loop over atoms in2
         do in2 = 1, nspecies
 
! Initialize initype and maxtype
          initype = 0
          if(interaction .ge. 15 ) initype = 1
          if(interaction .eq. 23 ) initype = 0
          maxtype = 0
          if (interaction .eq. 2) isub2c = nssh(in1)
          if (interaction .eq. 3) isub2c = nssh(in2)
          if (interaction .eq. 4) isub2c = nssh(in2)
          if (interaction .eq. 15) isub2c = nssh(in1)
          if (interaction .eq. 16) isub2c = nssh(in2)
          if (interaction .eq. 17) isub2c = nssh(in2)
          if (interaction .eq. 18) isub2c = nssh(in1)
          if (interaction .eq. 19) isub2c = nssh(in2)
          if (interaction .eq. 20) isub2c = nssh(in1)
          if (interaction .eq. 21) isub2c = nssh(in2)
          if (interaction .eq. 22) isub2c = nssh(in2)

          if (itheory .eq. 1) maxtype = isub2c
        ! Harris case for average density
      if (interaction .ge. 15 .and. interaction .le. 22) maxtype = isub2c 
          do isorp = initype, maxtype

! Append the number of subtypes on root if there is more than one file for a
! given pair of atoms and a given interaction.  The result will be root_isorp
           root_isorp = root
           if (isub2c .ge. 1) then    
            write (extension,'(''_'',i2.2)') isorp
            root_isorp = append_string (root,extension)
           end if

! Append the nuclear charges at root, the result will be root_isorp.nz1.nz2.dat
           nzx1 = nzx(in1)
           nzx2 = nzx(in2)
           write (extension,'(''.'',i2.2,''.'',i2.2)') nzx1, nzx2
           filename = append_string (root_isorp, extension)
           write (extension,'(''.dat'')')
           filename = append_string (filename, extension)
!           if (isorp .eq. initype) write (*,'('' Opening data file: '',a100)') filename
           open (unit = iounit, file = filename, status = 'old')
           call readheader_2c (interaction, iounit, nsh_max, numz, rc1, rc2, &
     &                         zmin, zmax, npseudo, cl_pseudo)
           if (numz .gt. nfofx) then
            write (*,*) ' numz = ', numz, ' in read_2c.f90'
            write (*,*) ' nfofx = ',nfofx
            write (*,*) ' Fix this parameter and recompile! '
            stop
           end if
 
! For the Non-local pseudopotential ONLY.
           if (interaction .eq. 5) then
            do issh = 1, nsshPP(in2) 
             cl_PP(issh,in2) = cl_pseudo(issh)
            end do
           end if
 
! Here are the data file characteristics: number of points and the grid range
           itype = ind2c(interaction,isorp)
           z2cmax(itype,in1,in2) = zmax
           numz2c(itype,in1,in2) = numz
 
! The array index_max2c(in1,in2) is the number of non-vanishing matrix elements
! for a general 2-center integral that are stored in the field
! twocint(index_max2c(in1,in2),numz,j2x,in1,in2) at numz different bond charge
! distances.
           num_nonzero = index_max2c(in1,in2)
 
! For pseudopotential.
           if (interaction .eq. 5) num_nonzero = index_maxPP(in1,in2)
 
! For the vna_atom and xc_atoms cases, the number of interactions is equal
! to the number of shells for in1,in1 pair not in1, in2 pair.  This is because
! the wavefunctions are both located only on in1!
           if (interaction .eq. 4 .or. interaction .eq. 7)                   &
     &      num_nonzero = index_max2c(in1,in1)
! JIMM
           if (interaction .eq. 10)  num_nonzero = index_max2cDipY(in1,in2)
           if (interaction .eq. 11)  num_nonzero = index_max2cDipX(in1,in2)
!
           if (interaction .eq. 17)  num_nonzero = index_max2c(in1,in1)
           if (interaction .eq. 22)  num_nonzero = index_maxS(in1,in1)

! Special case for coulomb (short-range) part.  It's special because it is not
! a matrix element.  This in interaction = 12
           if (interaction .eq. 12 .or. interaction .eq. 14)                 &
     &      num_nonzero = nssh(in1)*nssh(in2)
! Special case for spehrical density approximation.  It's special because 
! it has different size of a matrix element.  
! This in interaction = 20,21,22
           if (interaction .eq. 20 .or. interaction .eq. 21)                 &
     &       num_nonzero = index_maxS(in1,in2)
           if (interaction .eq. 23) num_nonzero = index_maxS(in1,in2)
          
           call readdata_2c (interaction, iounit, num_nonzero, numz, zmax,   &
     &                       itype, in1, in2, ioff2c)

           close (unit = iounit)
          end do
         end do ! nspecies
        end do ! nspecies

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ========================================================================


        return
        end subroutine read_2c
