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


! doscentros.f90
! Program Description
! ===========================================================================
!      This subroutine calculates the (two-center) matrix elements (mu,nu).
! There used to be five different routines that did this for all of the
! two-center interactions - doscentros.f, dosxcatm.f, dosxcontop.f,
! dosenatm.f, and dosenontop.f.  These have now all been reduced to one
! routine in order to make Fireball more lean.
!
!      This routine also calculates the derivative with respect to the
! position of the orbital of the BRA.
!
! ===========================================================================
! Original code written by Jose Ortega.
!
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
        subroutine test2c (interaction, isub, in1, in2, in3)
        use dimensions
        use interactions
        use integrals
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: interaction
        integer, intent (in) :: isub
 
! The variable interaction is the type of two-center integral that is needed.
! The variable isub is the subtype of interation
! Here is information from read2c
!
!           interaction    subtype
!               1            isorp = 0, 8    average density 
!               1            0               overlap
!               2            0               vna_ontopl
!               2            isorp = 1, 8    vna_ontopl shell isorp
!               3            0               vna_ontopr
!               3            isorp = 1, 8    vna_ontopr shell isorp
!               4            0               vna_atom
!               4            isorp = 1, 8    vna_atom  shell isorp
!               5            0               non-local
!               6            0, 4            xc_ontop
!               7            0, 4            xc_atom
!               8            0, 4            xc_correction
!               9            0               z-dipole
!               10           0               y-dipole
!               11           0               x-dipole
!               12           0               coulomb
!               13           0               kinetic
!               14           0               extended-hubbard
!               15           isorp = 1, 8    density_ontopl
!               16           isorp = 1, 8    density_ontopr
!               17           isorp = 1, 8    density_atom  
!  spherical density (only doscentrosS)
!               18           isorp = 1, 4    sph dens ontopl
!               19           isorp = 1, 4    sph dens_ontopr
!               20           isorp = 1, 4    sph den_atom
!               21           isorp = 0       sph overlap
 
! Derivatives are computed if and only if iforce = 1.

        integer, intent (in) :: in1
        integer, intent (in) :: in2
        integer, intent (in) :: in3

 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer inu
        integer index
        integer iforce
        integer igrid
        integer numz
        integer num_nonzero
        integer itype
 
! -slist = output list of matrix elements
! -dslist = output list of derivatives of matrix elements
        real, dimension (ME2c_max) :: dslist
        real, dimension (ME2c_max) :: slist
        real zmax
        real distance 
        real dx

        character(80)   :: title
        character(2)   :: tit1
        character(2)   :: tit2

        logical switch

! Procedure
! ===========================================================================
! For the atom case, in3 = in1, but for everything else in3 = in2.
! For the ontop case, in2 = in1 (left) or in3 (right).
! Initialize sm, scam and sx to zero.
! < n1 | n2 | n3 >
        iforce = 0
        switch = .true.
        if(interaction .eq. 2) switch = .false.
        if(interaction .eq. 15) switch = .false.
        if(interaction .eq. 18) switch = .false.

        write (*,*) ind2c(3,:)
        write (*,*) ind2c(4,:)
        write (*,*) ind2c(5,:)
        itype = ind2c(interaction,isub)
        write (*,*) itype, interaction, isub
        zmax = z2cmax(itype,in1,in2)
        numz = numz2c(itype,in1,in2)
        num_nonzero = index_max2c(in1,in2)
        dx = zmax / numz
        distance = 0.0d0
        write (tit1,'(i2.2)') interaction
        write (tit2,'(i2.2)') isub
        title = 'fdata_'//tit1//'_'//tit2//'.dat'
        open (unit=10, file=title,status='unknown')
! This subroutine calls the subroutine intrp1d as needed to find the value of
! the matrix elements for any given atomic separation distance.
! -slist = output list of matrix elements
! -dslist = output list of derivatives of matrix elements
        do igrid = 1,numz 
         do index = 1, index_max2c(in1,in3)
          if ( switch ) then
           call interpolate_1d (interaction, isub, in1, in2, index, iforce,   &
     &                         distance, slist(index), dslist(index))
          else
           call interpolate_1d (interaction, isub, in1, in3, index, iforce,   &
     &                         distance, slist(index), dslist(index))
          end if
         end do
         distance = distance + dx
          do index=1,num_nonzero
            write (10,100, advance="no") slist(index)
          enddo 
          write(10,*) '  '
        end do

        close (10)
 
! Format Statements
! ===========================================================================
100 format(f15.9) 
        return
      end subroutine test2c
