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


! readphi.f90
! Program Description
! ===========================================================================
!       This reads in the dynammat.vib file
!
! ===========================================================================
! Original code from Otto F. Sankey with modification by Alex A. Demkov
! and Jose Ortega
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
        subroutine readphi (natoms, ratom, nstepi, nstepf)
        use dynamo
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (inout)     :: nstepi
        integer, intent (inout)     :: nstepf
        integer, intent (in)        :: natoms
        real,  dimension (3, natoms), intent (in) :: ratom

! Ouput
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom,ix

        logical readsup   ! list of free atoms ??

! Procedure
! ===========================================================================
        write (*,*) ' You are doing a dynamical matrix calculations'
        write (*,*) ' of the vibrational spectrum.'
        write (*,*) ' Read the parameters from dyn.optional'

! open input file
        open (file ='dyn.optional', unit=17, status='old')

! read elementary displacement
        read (17,*) u0disp
        write (*,*) ' Elementary  displacement = ', u0disp

! read vector of valid directions 
        read (17,*) (ndvec(ix), ix = 1, 3) 
        ndx = ndvec(1) + ndvec(2) + ndvec(3) 
        write (*,100) ndx, (ndvec(ix), ix = 1, 3) 
        if (ndx .le. 0 .or. ndx .gt. 3 ) then
         write (*,*) ' bad nspace !!'
         stop
        end if

! read name of output file 
        read (17,*) filephi
        write (*,200) filephi
  
! read flag if list of atoms or all atoms will be moved
        read (17,*) readsup

        if (.not. readsup) then 
         natoms_dm = natoms
         allocate (jatoms_dm(natoms))

! assign number
         do iatom = 1, natoms
          jatoms_dm(iatom) = iatom
         end do
        else

! read number of moved atoms 
         read (17,*) natoms_dm
         allocate (jatoms_dm(natoms))

! read the atom list 
         do iatom = 1, natoms_dm
          read (17,*) jatoms_dm(iatom)
         end do
        end if

! close input file
        close (unit = 17)

! allocate 
        allocate ( phidm(ndx*natoms,ndx*natoms) )
        allocate ( u0vec (ndx*natoms) )
        allocate ( ftot1 (3,natoms) )
        allocate ( ratom0 (3,natoms) )
      
! save initial atomic positions
        ratom0 (:,:) = ratom (:,:)
  
  
! set number of time steps
        nstepi = 1
        nstepf = 2*ndx*natoms_dm
        write (*,300) nstepf
        write (*,400)




! set flag
! remember we have two displacement ('+' and '-')for each elementary 
! displacement to avoid certain problems with symmetries 
! (e.g. diatomic molecule)
  ldynamo = .true.
  ltime = 1 

! Format Statements
! ===========================================================================
100    format ('   nspace = ', i3, '      x,y,z directions :',3i2)
200    format ('   Dynamical matrix will be written into ',a30,' file')
300    format ('   Total number of steps : ',i4)
400    format ('  ===========================================================')

  return
end subroutine readphi
