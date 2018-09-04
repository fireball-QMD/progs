! copyright info:
!
!                             @Copyright 2006
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Academy of Sciences of the Czech Republic - Pavel Jelinek

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
! C.Gonzalez Pascual, UAM, Spain

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


! writeout_dosng.f90
!If dosng = 1 there must be a dosng.optional file with the following structure:
!---------------------------------------
!etaL   !eta factor for the Lorentzian
!E_ini,S,step !Initial energy, number of energies and step
!dosats !number of atoms over which to project the DOS
! dosat(k) !list of atoms over which we project the DOS
!The output file consists of S rows and each row is of the form:
!E, D,
!with E an energy value and D the DOS at that energy.
!#########
!If dosng = 2, this program outputs files giving the Lowdin coefficients of 
!of user-requested eigenstates. The states.optional file in this case has the
!following structure
!-----------------------------------------
!Ns !Number of eigenstates to be printed. There will be a file for each eigenstate
!number of eigenstate1
!number of eigenstate2
!...
!number of eigenstate Ns
!Each output file consist of Natoms rows, where Natoms is the total number of orbitals. 
!In the file state_x.out, each row is of the form:
!n c
!where n is the atom number and c is the coefficient (squared)
!of the x-th state over the n-th atom (using Löwdin orbitals).
!#########
!If iwrtdos = 3 then we calculate both the DOS and the projection of the
!electronic states.
!#############################################################################
! Program Description
! ===========================================================================
! This routine writes the calculates and writes the DOS of the atoms we want
! ===========================================================================
! Code written by:
! ==========================================================================

        subroutine writeout_dosng( )

        use dimensions
        use interactions
        use neighbor_map
        use module_dos
        use configuration
        use charges
        use density
        use outputs

        implicit none

! Argument Declaration and Description
! Input

! Output

! Local Parameters and Data Declaration
! ===========================================================================
       real, parameter :: pi=3.141592653589793230d0

! Local Variable Declaration and Description
! ==========================================================================
        integer         :: iatom
        integer         :: jatom
        integer         :: iorb
        integer         :: jorb
        real            :: energy
        real            :: E
        integer k       !loop variable
        integer alpha   !loop variable
        integer B       !number of eigenvalues
        integer S       !number of points in the energy grid
        integer ii      !index for the loop over the grid points
        integer jj      !index for the loop over all the eigenvalues
        real    subtot  !temporary variable for the sum of coefficients
        real    tot     !temporary variable storing values of DOS
        real    etaL     !width of the lorentzian
        real    step    !separation between points in the energy grid
        integer, dimension(natoms) :: dosat       !array storing the atoms to
                                                  !be taken into account
        integer, dimension(natoms) :: dosatsp     !array storing the species
                                                  !of the atoms in DOS
                                                  !calculation
        integer, dimension(numorb_max) :: dosor   !array storing the orbitals to
                                                  !be taken into account
        integer dosats  !total number of atoms in the DOS calculation
        integer dosorbs !total number of orbitals in the DOS calculation
        integer in1     !index storing atomic species for (dosng=1)
        real    E_ini   !initial value of energy in the grid (dosng=1)
        integer Nstates !Number of molecular states to be printed out (dosng=2)
        integer, dimension(:), allocatable :: states !dosng=2
        character(len=60) :: windex !dosng=2
! Procedure
! ===========================================================================
!Initialize stuff
!FIRST ORDER OF BUSSINESS: Read info from dosng.optional
!Format for dosng.optional file

write(*,*)  'HI FROM DOSNG!'
write(*,*)  'This subroutine computes the Atom-projected Density of States without computing the Green function'

if ((iwrtdosng .eq. 1) .or. (iwrtdosng .eq. 3)) then
open(unit = 172, file = 'dosng.optional', status = 'old')
read(172,*) etaL   !eta factor for the Lorentzian
read(172,*) E_ini,S,step !Initial energy, number of energies and step
read(172,*) dosats !number of atoms over which to project the DOS
do k = 1,dosats    

read(172,*) dosat(k) !list of atoms over which we project the DOS
dosatsp(k) = imass(dosat(k))  !with this we store the species of the k-th atom
! in the DOSNG.OPTIONAL file

end do !end do do k = 1,dosats
close(172)

!##### TEST VALUES OF E_ini, S and step

write(*,*) 'DONE READING INPUT FILE'

!E_ini = E_KS(1,1)

!step=0.5

!S = (E_KS(norbitals,1)-E_ini)/step


!##### END OF TEST

! Open file dosng
open(unit = 173, file = 'dosng.out', status = 'unknown')
! OJO!!! ONLY FOR ICLUSTER = 1 FOR THE TIME BEING
!Eks: array with all the eingenvalues. Call B to the number of eigenvalues, which is the dimension of Eks.

write(173,*) '--------NEW STEP--------, '

do ii = 1,S !loop over grid points
tot=0.0d0
E=E_ini+ii*step
do iorb = 1,norbitals !loop over eigenvalues

subtot=0.0d0

do k = 1,dosats
iatom = dosat(k)
in1 = dosatsp(k)
do jorb = 1,num_orb(in1)
alpha = degelec(iatom)+jorb
subtot=subtot+dngcof(alpha,iorb,1)**2
!this is the DOS in terms of Lowdin coefficients. CAREFUL! We need to use
!the TRANSPOSE MATRIX! 
end do !end do iorb = 1,num_orb(in1)
end do ! end do k = 1,dosats ( over orbitals to consider in the DOS we are calculating)

tot = tot+subtot*etaL*(1/(etaL**2+(E-E_KS(iorb,1))**2))*(1/pi) 
!Here we multiply by the Lorentzian centered around E_KS(iorb)
!The second index of E_KS is ikpoint. For the time being this is
!only for icluster = 1

end do !end do iorb = 1,norbitals
write(173,*) E,tot
end do !end do ii = 1,S
close(173)


!###########END OF DOSNG = 1
!######################################3
!########### NOW DOSNG = 2: 
else if ((iwrtdosng .eq. 2) .or. (iwrtdosng .eq. 3)) then   

open(unit = 172, file = 'states.optional', status = 'old')
read(172,*) Nstates   !number of molecular states to print out
allocate(states(Nstates))
do k = 1,Nstates
read(172,*) states(k) !list of molecular states to print out
end do !end do do k = 1,Nstates
close(172)

do k = 1,Nstates
iorb = states(k)
write(windex,'(i0)')iorb
open(unit = 174 , file = 'state_'//trim(windex)//'.out', status = 'unknown')
write(174,*) '--------NEW STEP--------, '
do iatom = 1,natoms
in1=imass(iatom)
tot = 0
do jorb = 1,num_orb(in1)
alpha = degelec(iatom)+jorb
tot = tot+dngcof(alpha,iorb,1)**2
end do !end do jorb = 1,num_orb(in1)
write(174,*) iatom, tot
end do !end do iatom = 1,natoms

close(174)
end do !end do k = 1,Nstates
end if 

!OLD############3
! Format Statements
! ===========================================================================
100     format (3f16.9, 2x, 2i2) 
200     format (2f16.6, i5) 

      return
  
      end subroutine writeout_dosng
