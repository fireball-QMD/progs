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
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! scf_loop_ks.f90
! Program Description
! ===========================================================================
!       This routine performs Kohn-Sham SCF loop
!
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine scf_loop_ks (itime_step)
 
        use options 
        use scf
        use energy
        
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
       integer, intent (in) :: itime_step

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
     
! Procedure
! ===========================================================================


! ===========================================================================
!                          Start loop over scf cycles 
! ===========================================================================
! Loop over self-consistent cycles if necessary.
! max_scf_iterations is the maximum number of times it will try to iterate.
         if (itheory .ne. 0) scf_achieved = .false.
         do Kscf = 1, max_scf_iterations
           write (*,*) '  '
           write (*,*) '  '
           write (*,*) ' ****************************************************** '
           if (itheory .eq. 0) then
            write (*,*) ' Begin time step = ', itime_step 
           else 
            write (*,*) ' Begin time step = ', itime_step, ' scf step = ', Kscf
           end if
           write (*,*) ' ****************************************************** '
           write (*,*) '  '

! ASSEMBLE HAMILTONIAN
          call assemble_h_ks ()

! ===========================================================================
! Exact Diagonalization: solve for the band-structure energy - call diag_k
! ===========================================================================

! Do k-loop diagonalization
           write (*,*)  ' Call diag_k'
           call diag_k_ks ()

! ===========================================================================
!               compute the density matrices 
! ===========================================================================
! Compute the density matrices. The results rho and cape are computed.
          call build_rho_ks (itime_step) 
          
! ===========================================================================
! At this point we are only calculating the total band-structure energy
! All the rest of the stuff that was written out here is moved down below
! until a SCF solution is obtained.
          if (itheory .eq. 0) exit

          if (scf_achieved) then
            write (*,*) '  '
            write (*,*) ' ------- THE TOTAL BAND STRUCTURE ENERGY ------- '
            write (*,*) '  '
            write (*,400) ebs, itime_step, Kscf
            write (*,*) ' ----------------------------------------------- '
            exit
          else if (.not. scf_achieved .and. Kscf .eq. max_scf_iterations) then
            write (*,*) '  '
            write (*,*) ' ------- THE TOTAL BAND STRUCTURE ENERGY ------- '
            write (*,*) '  '
            write (*,401) ebs, itime_step, Kscf
            write (*,*) ' ----------------------------------------------- '
            exit
          else if (.not. scf_achieved) then
            write (*,*) '  '
            write (*,*) ' ------- THE TOTAL BAND STRUCTURE ENERGY ------- '
            write (*,*) '  '
            write (*,402) ebs, itime_step, Kscf
            write (*,*) ' ----------------------------------------------- '
          end if
          
         end do  ! end do Kscf
 
! ===========================================================================
!                          End loop over scf cycles 
! ===========================================================================

! Deallocate Arrays
! ===========================================================================

 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
400     format (' EBS = ', f20.6, ' Time step = ', i5, ' SCF step = ',       &
     &          i4, ' SCF success! ')
401     format (' EBS = ', f20.6, ' Time step = ', i5, ' SCF step = ',       &
     &          i4, ' SCF failure! ')
402     format (' EBS = ', f20.6, ' Time step = ', i5, ' SCF step = ', i4)


        return
        end subroutine scf_loop_ks
 
