! Copyright info:
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


! initcDFT.f90
! Program Description
! ===========================================================================
! ===========================================================================
! Code written by:
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine initcDFT ()

   use scf
   use interactions
   use charges 
   use MD
   use options
   use kpoints
   use density 
   use nonadiabatic
   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

! Output

! Local Parameters and Data Declaration
! ===========================================================================


! Local Variable Declaration and Description
! ===========================================================================
        real, parameter :: tol = 1.0d-4
        real diffq
        real qztot
        real qcharge
        real norm
        integer ikpoint
        integer iband
        integer jband
        integer nfermi
        integer read_eh
        integer i
        integer iele
        integer jele

! Procedure
! ===========================================================================
! Allocate bbnkre blowre and eigen_k

         open (unit = 3005, file = 'en_cDFT.dat', status = 'unknown')
         open (unit = 4001, file = 'occ_gs.dat', status = 'unknown')
         open (unit = 4002, file = 'occ_gs_check.dat', status = 'unknown')
         open (unit = 4003, file = 'occ_es_smear.dat', status = 'unknown')
         open (unit = 4004, file = 'occ_na.dat', status = 'unknown')

         allocate (bbnkre (norbitals, norbitals, nkpoints))
         allocate (blowre (norbitals, norbitals, nkpoints))
         allocate (eigen_k (norbitals, nkpoints))

! First we do ground state scf calculations we need reference wf
         cDFT_active = .false.
         itime_step_g = 0
         call    scf_loop (0)
         qztot = ztot
         write (*,*) '  ---- Initialize cDFT -----'

! Read from the input file of cDFT 
         open (unit = 22, file = 'cDFT.optional', status = 'old' )
         write (*,*) ' Reading from cDFT.optional file! '

! Read number of wf stored for projection
         read (22,*) read_eh
         read (22,*) gs_scf 
         read (22,*) tempfe_eh 
         read (22,*) n_hist

! Weights used for construction reference function from history
         allocate (wf_weight (n_hist))
         do i=1,n_hist 
           read(22,*) wf_weight(i)    
         end do

! Number of states to be invoked in tracking
         read(22,*) nele
         allocate (map_ks (nele) )
         allocate (iocc (nele))
         allocate (hist_fix (nele))
allocate (map_ks_o (nele) )

! number of states involved in transitions
         do iele = 1, nele
           read (22,*) iband, iocc(iele)
           map_ks(iele) = iband
           map_ks_o(iele) = iband
         end do
         close (22)
!         if (read_eh .eq. 1) call init_excitation()
         allocate (blowre_hist(norbitals,nele,n_hist))       
!         allocate (Wmu_glob (nele))
         allocate (Wmu_glob (norbitals))

! store initial history of wf (in this case all stored wf are the same) 
         do iband = 1, norbitals  
            do iele = 1,nele
              do i = 1,n_hist
                blowre_hist(iband,iele,i)=blowre(iband,map_ks(iele),1)
              end do
            end do
         end do   
         flag_proj = 0
         flag_es = 0
         allocate (sgn_wf (nele))
         allocate ( foccupy_na ( norbitals, nkpoints) )
         allocate ( ioccupy_na ( norbitals, nkpoints) )
         allocate ( foccupy_wr_gs ( norbitals, nkpoints) )   
         allocate ( foccupy_wr_proj  ( norbitals, nkpoints) )
         allocate ( foccupy_wr_gs_check ( norbitals, nkpoints) )  

         allocate ( foccupy_na_o ( norbitals, nkpoints) )
         allocate ( ioccupy_na_o ( norbitals, nkpoints) )
         ioccupy_na = 0
         foccupy_na = 0.0d0

         nfermi = int(qztot) / 2
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


! change occupations using mdet.input values  (iocc(iele))        
         do ikpoint = 1, nkpoints
          do iele = 1, nele
           iband = map_ks(iele)
           foccupy_na(iband,ikpoint) = iocc(iele)*0.5d0
           ioccupy_na(iband,ikpoint) = iocc(iele)
          end do
         end do

         foccupy_na_o(:,:) = foccupy_na(:,:)
         ioccupy_na_o(:,:) = ioccupy_na(:,:)

         do iband = 1, norbitals 
          if (foccupy_na(iband,1) .gt. 0.0 ) then
             ioccupy_na(iband,1) = 1
          else
             ioccupy_na(iband,1) = 0
          end if  
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

! ===========================================================================
200     format (' Band n = ', i4, ' k-points: ioccupy = ', i2)
201     format (' Band n = ', i4, ' foccupy = ', 2f12.8 )
   return
 end subroutine initcDFT

