! copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
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

! MCsolar.f90
! Program Description
! ===========================================================================
!       This routine assemble forces for MDET/McWeda
!
! JOM : adapted from getforces_mcweda to also calculate 
! the gradient of the Hamiltonian
! G < mu | H | nu > and overlap contributions for the calculation
! of the nonadiabatic couplings
! ===========================================================================

! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine MCsolar (itime_step)

        use options
        use outputs
        use mpi_main
        use configuration
        use forces
        use nonadiabatic

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: itime_step

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: dimax = 0.1d0
      
! Local Variable Declaration and Description
! ===========================================================================
        integer ielec
        integer jelec
        integer iatom
        integer ix
        integer in1
        real tmp
        real, dimension (3,natoms) :: dr
        real xrand

! Procedure
! ===========================================================================
! Initialize seed for random_number
        call random_seed
        
! calc dnac F/m*dij    
        do ielec = 1, nele
          do jelec = 1, nele
            tmp = 0.0d0
            do iatom = 1, natoms
              do ix = 1,3
                tmp = tmp + ftot(ix,iatom)/xmass(iatom)*gks(ix,iatom,ielec,jelec)
              enddo ! do ix
            enddo  ! do iatom
            dnac(ielec,jelec) = tmp
          enddo  ! do jelec
        enddo ! do ielec

        write(*,*) dnac(1,1),dnac(1,2)
        write(*,*) dnac(2,1),dnac(2,2)

! first step
        if (itime_step .eq. 1) then
! save optimal configuration
          do iatom = 1, natoms
            do ix = 1,3
             ratom_opt(ix,iatom) = ratom(ix,iatom)
            enddo !do ix
          enddo ! do iatom
          do ielec = 1,nele
            do jelec = 1,nele
              dnac_opt(ielec,jelec) = dnac(ielec,jelec)
            enddo ! do jelec
          enddo ! do ielec    
! dump *.xyz         
          if (iwrtxyz .eq. 1) then 
            open (unit = 18, file = 'answer.xyz', status = 'unknown') 
            write (18,*) natoms
            write (18,200) dnac(1,2)   
            do iatom = 1, natoms
              write (18,701) symbol(iatom), ratom_opt(:,iatom) + ximage(:,iatom)
            enddo
            close (unit = 18)
          endif ! if (iwrtxyz)  
        else
          if (dnac_opt(1,2) .lt. dnac(1,2)) then
            ratom_opt = ratom
            dnac_opt = dnac
! dump *.xyz         
            if (iwrtxyz .eq. 1) then 
              open (unit = 18, file = 'answer.xyz', status = 'unknown') 
              write (18,*) natoms
              write (18,200) dnac(1,2)  
              do iatom = 1, natoms
                write (18,701) symbol(iatom), ratom_opt(:,iatom) + ximage(:,iatom)
              enddo
              close (unit = 18)
            endif ! if (iwrtxyz)                              
          endif ! if (dnac_pot)
        endif ! if (itime_step)

! generate random displacement vector dr       
        do iatom = 1, natoms
          do ix = 1,3
! Random numbers for Monte-Carlo
            call random_number(xrand)  
            dr(ix,iatom) = dimax*(xrand-0.5d0)
          enddo ! do ix
        enddo ! do iatom
        
! Random numbers for Monte-Carlo
        call random_number(xrand)
        dr = dr*xrand
        
        do iatom = 1, natoms
          do ix = 1,3
            ratom(ix,iatom) = ratom(ix,iatom) + dr(ix,iatom)
          enddo ! do ix
        enddo ! do iatom
        
! debug dump        
        write (218,*) natoms
        write (218,200) dnac(1,2)  
        do iatom = 1, natoms
          write (218,701) symbol(iatom), ratom(:,iatom) + ximage(:,iatom)
        enddo
                
! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))
200     format ('# NAC =',f8.4 )
701 format (2x, a2, 3(2x,f12.6))

        return
        end subroutine MCsolar

