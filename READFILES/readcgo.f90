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
! Ohio State University - Dave Drabold
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


! readcgo.f90
! Program Description
! ===========================================================================
!       This routine reads in the information from the cgopt.optional file.
!
! ===========================================================================
! Code rewritten by:
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine readcgo (natoms, iforce)
  
        use optimization 
        use fragments
        use interactions

        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms

! Output
        integer, intent (inout) :: iforce
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ix

! Procedure
! ===========================================================================
        write (*,*) '  '
        write (*,100)
        write (*,*) ' You have chosen conjugate gradient method'
        write (*,*) ' to find structure with minimal energy.'
        write (*,*) ' Read information from the cgopt.optional file. '
        write (*,*) '  '
           
        open (unit = 77, file = 'cgopt.optional', status = 'old')

! Initialize conjugate gradient parameters
        write (*,*) ' In the event that you are searching for an energy '
        write (*,*) ' minimum, you can also choose energy and force '
        write (*,*) ' tolerances - energytol and forcetol. These tolerances '
        write (*,*) ' will stop execution after they are achieved. Input a '
        write (*,*) ' number for all of the following, even if option is not '
        write (*,*) ' being used! '
                 
! Read multiple constraints for the steps in search of min 
        read(77,*) cg_drmax

! Read the scale to reduce the search step if e1 < e2 
        read(77,*) cg_dummy

! Read tolerance of the convergence (difference of the consequent values) 
        read(77,*) energy_tol

! Read tolerance of convergence for forces
        read(77,*) force_tol

! Read maximum number of CG steps 
        read(77,*) cg_maxstep

! Read maximum number of steps of internal CG loop  
        read(77,*) cg_minint

! Read switch controling quenching refinement of the CG results
        read(77,*) icg2md

! write out the resume of the variables
        write(*,*) ' Parameter of CG optimization : '
        write(*,*) ' ======================================'
        write(*,*) '  drmax              = ',cg_drmax
        write(*,*) '  dummy              = ',cg_dummy
        write(*,*) '  etot_tol           = ',energy_tol
        write(*,*) '  for_tol            = ',force_tol
        write(*,*) '  max. CG steps      = ',cg_maxstep
        write(*,*) '  min. int. steps    = ',cg_minint
        write(*,*) '  switch to MD       = ',icg2md

! Reset counter of steps
        cg_iter = 0

        write(*,*)'   '
        write(*,*)'  ============= End cgtol.input ============='
        write(*,*)'   '
        close (unit = 77)

! Set fixed atoms (mass.gt.1000)
        mask = 1.0d0
        freeParamsCount = 3*natoms
        if (numfrags .ne. 0 ) then 
         do iatom = 1, natoms
          do ix = 1,3
           if (fragxyz(ix,iatom) .eq. 1)then
           	 mask(ix,iatom) = 0.0d0
           	 freeParamsCount = freeParamsCount - 1
           endif
          end do
         end do

         write(*,*) '  CG_OPTIMIZATION:    LIST OF THE FIXED ATOMS '
         write(*,*) ' ============================================='
         do iatom = 1, natoms
          if (fraggots(iatom) .eq. 1) then 
           write (*,200) iatom, ' FIXED ',(fragxyz(ix,iatom), ix=1,3)
          else
           write (*,201) iatom, ' FREE '
          end if
         end do
         write (*,*) ' ============================================='
        else 

         write(*,*) ' ==============================================='
         write(*,*) '  CG_OPTIMIZATION:  ALL ATOMS ARE FREE TO MOVE  '
         write(*,*) ' ==============================================='

        end if ! 
! Allocate auxillary vectors of cg algorithm 
        allocate (x0 (3, natoms))
        allocate (f0 (3, natoms))
        allocate (Q0 (nsh_max, natoms))
        allocate (g (3, natoms))
        allocate (h (3, natoms))
      
        cgiter = 0
        initer = 0

! We need forces for first step
        iforce = 1

! Set icg_status
        istatus = 0
       
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
200     format ('iatom = ', i3, '  status = ',a8,3i3) 
201     format ('iatom = ', i3, '  status = ',a8) 

        return
        end subroutine readcgo
