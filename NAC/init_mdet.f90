! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
! Universidad Autonoma de Madrid - Jose Ortega

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

 
! init_mdet.f90
! Program Description
! ===========================================================================
!       This routine gives the initial state for molecular dynamics with
!       electronic transitions (MDET) (nonadiabatic calculation)
!
! ===========================================================================
! Code written by:
! Jose Ortega Mateo
! Departmento de Fisica Teorica de la Materia Condensada
! Universidad Autonoma de Madrid
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine init_mdet (natoms)
!       use dimensions
!       use constants_fireball
        use kpoints
        use interactions
        use density
        use charges
        use nonadiabatic
        use options
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: natoms
 
!       real, intent(inout), dimension (norbitals, nkpoints) :: eigen_k
 
! Output
!       integer, intent(out), dimension (norbitals, nkpoints) :: ioccupy_k
!       real, intent(out), dimension (norbitals, nkpoints) :: foccupy
 
! Local Parameters
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ikpoint
        integer iband
        integer jband
        integer nswitchs
        integer is
        integer ia
        integer iele
        integer jele
        integer istp
        integer nfermi
        integer, dimension (:), allocatable :: iocc
 
        real qztot
        real qcharge
        real norm

        complex a0
        complex a1
        complex caux
        complex cnorm

        real, parameter :: tol = 1.0d-4

        logical read_switchs
 
! Procedure
! ===========================================================================
! Allocate bbnkre blowre and eigen_k
        allocate (bbnkre (norbitals, norbitals, nkpoints))
        allocate (blowre (norbitals, norbitals, nkpoints))
        allocate (eigen_k (norbitals, nkpoints))
! Initialize some things
        a0 = cmplx(0.0d0,0.0d0)
        a1 = cmplx(1.0d0,0.0d0)

! Add in the qstate to the total charge
! HAO--possible fix
         qztot = ztot 
! JOM-info : number of states involved in transitions
         open (unit = 22, file = 'mdet.input', status = 'old' )
         write (*,*) ' Reading from mdet.input file! '
! JOM-info number of switches
         read(22,*) nele
         allocate (map_ks (nele) )
         allocate (map_proj (nele))
         allocate (iocc (nele))

         do iele = 1, nele
           read (22,*) iband, iocc(iele)
           map_ks(iele) = iband
         end do
         close (22)

!---------------------------------------------------------------------

        allocate (gks (3, natoms, nele, nele))
        allocate (dnac (nele, nele))
! jel-nac        
        if (imdet .eq. 2) then
          allocate (ratom_opt(3,natoms))
          allocate (dnac_opt(nele,nele))
        endif
! JOM-q or norbitals_new

         allocate ( foccupy_na ( norbitals, nkpoints) )
         allocate ( ioccupy_na ( norbitals, nkpoints) )
!---------------------------------------------------------------------
! JOM-info:  ioccupy_na = 0, 1 or 2 depending on the occupation of the
! state
!---------------------------------------------------------------------
         allocate ( c_na (nele, nele, nkpoints) )
! Initialize the occupation numbers and foccupy to zero.
        ioccupy_na = 0
        foccupy_na = 0.0d0
 
        nfermi = int(qztot) / 2
!       write(*,*)'nfermi,qztot',nfermi,qztot,int(qztot),2*nfermi

        do ikpoint = 1, nkpoints
         do iband = 1, nfermi
          foccupy_na (iband,ikpoint) = 1.0d0
          ioccupy_na (iband, ikpoint) = 2
         end do
        end do
        if (int(qztot) .gt. 2*nfermi) then
         do ikpoint = 1, nkpoints
          ioccupy_na(nfermi+1,ikpoint) = 1
          foccupy_na (iband,ikpoint) = 0.5d0
         end do
        end if

! check
        qcharge = 0.0d0
        do ikpoint = 1, nkpoints
         do iband = 1, norbitals
          if (ioccupy_na(iband,ikpoint) .ne. 0) then
        qcharge = qcharge + 2.0d0*foccupy_na(iband,ikpoint)*weight_k(ikpoint)
          end if
         end do
        end do
        if (abs(qcharge - qztot) .gt. tol) then
         write (*,*) '          qcharge = ', qcharge
         write (*,*) '          qztot = ', qztot
         write (*,*) 'must stop in subroutine init_mdet 1'
         stop
        end if


! Now, initialize the electronic states (c_na)
        
!       we must define the set of eigenstates for which we are going to
!       follow their time-evolution via the c_na(t)
!       for example:
        
        c_na = a0
        do ikpoint = 1, nkpoints
         do iele = 1, nele
        c_na (iele, iele, ikpoint) = a1
!       write(*,*)'c_na',iband,c_na (iband, iband, ikpoint)
         end do 
        end do 
 
! check Normalization !
        do ikpoint = 1, nkpoints
         do iele = 1, nele
         cnorm = a0
          do jele = 1, nele
          caux =c_na(iele,jele,ikpoint)
          cnorm = cnorm + caux*conjg(caux)
          end do
         norm = cabs (cnorm)
         write(*,*)'Norm of initial states',iele, norm
         end do
        end do
        
! change occupations using mdet.input values  (iocc(iele))        
        do ikpoint = 1, nkpoints
         do iele = 1, nele
          iband = map_ks(iele)
          foccupy_na(iband,ikpoint) = iocc(iele)*0.5d0
          ioccupy_na(iband,ikpoint) = iocc(iele)
         end do
        end do
        
! check
        qcharge = 0.0d0
        do ikpoint = 1, nkpoints
         do iband = 1, norbitals
          if (ioccupy_na(iband,ikpoint) .ne. 0) then
        qcharge = qcharge + 2.0d0*foccupy_na(iband,ikpoint)*weight_k(ikpoint)
          end if
         end do
        end do
        if (abs(qcharge - qztot) .gt. tol) then
         write (*,*) '          qcharge = ', qcharge
         write (*,*) '          qztot = ', qztot
         write (*,*) 'must stop in subroutine init_mdet 2'
         stop
        end if        
        
        
        
        
!--------------------------------------------------------------------
! write
         do iele = 1, nele
          write(*,*)'map_ks',iele,map_ks(iele)
         end do
        do ikpoint = 1, nkpoints
         do iband = 1, norbitals
          write(*,*)'foccupy',iband,ioccupy_na(iband,ikpoint),foccupy_na(iband,ikpoint)
         end do
        end do
        do ikpoint = 1, nkpoints
         do iele = 1, nele
          do jele = 1, nele
        write(*,*)'c_na',iele,jele,c_na (iele, jele, ikpoint)
          end do 
         end do 
        end do 

! Format Statements
! ===========================================================================
 
        return
        end
