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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

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

   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

! Output

! Local Parameters and Data Declaration
! ===========================================================================


! Local Variable Declaration and Description
! ===========================================================================
      real diffq

! Format Statements
      write (*,*) '  ---- Initialize cDFT -----'
      open (unit = 33, file = 'cDFT.optional', status = 'old')
      read (33,*) id_hole
      read (33,*) id_elec
      read (33,*) occup_elec
      close (unit = 33)
      occup_hole = 1.0d0 - occup_elec
      write (*,200) id_hole, occup_hole
      write (*,300) id_elec, occup_elec

! allocate arrays
      allocate (wf_hole(norbitals, nkpoints))
      allocate (wf_elec(norbitals, nkpoints))

      cDFT_active = .false.
      
! only 2 MD steps (ground & excited state) & no forces
!      nstepf = 3
!      iforce = 0

! ===========================================================================
200     format (2x, ' Hole @ state : ',i8,' with occupancy ',f8.4)
300     format (2x, ' Electron @ state : ',i8,' with occupancy ',f8.4)

   return
 end subroutine initcDFT

