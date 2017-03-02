! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
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
! University of Regensburg - Juergen Fritsch

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! hoppings.f90
! Program Description
! ===========================================================================
! This routine writes the hoppings for several distances for the STM 
! calculation.  
! It generates the file tip_sample_aux y tip_sample_log_aux.   
! We active this option (iwrthop=1) when we use an script that runs the 
! program for each distance and saves the result for all the values in a same
! file (tip_sample) we will need in the STM calculations.         
! ===========================================================================
! Code written by:
! Cesar Gonzalez Pascual
! Departamento de Fisica Teorica de la Materia Condensada
! Facultad de Ciencias. Universidad Autonoma de Madrid
! E28049 Madird SPAIN
! FAX +34-91 3974950
! Office Telephone +34-91 3978648
! ==========================================================================

       subroutine writeout_hop

        use density
        use charges
        use dimensions
        use interactions
        use neighbor_map
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input

! Output

! Local Parameters and Data Declaration
! ===========================================================================
        real, parameter :: ha=27.21160d0

! Local Variable Declaration and Description
! ===========================================================================

        integer dg1
        integer dg2
        integer in1
        integer in2
        integer iatom
        integer jatom
        integer iorb
        integer jorb
        integer index
        integer issh1
        integer issh2
        integer l1
        integer l2
        integer imu
        integer inu
        integer num_nonzero
        real, dimension (:), allocatable :: hlist

! Procedure
! ===========================================================================

        do iorb = 1,norbitals
         write(49,'(99E16.5)') (real(hamk(iorb,jorb))/ha,jorb=1,norbitals)
        end do

        open (unit = 28, file = 'tip_sample_aux', status = 'unknown')
        open (unit = 27, file = 'tip_sample_log_aux', status = 'unknown')
        iatom = 1
        in1 = imass(iatom) 
        dg1 = degelec(iatom)
        jatom = 2
        in2 = imass(jatom) 
        dg2 = degelec(jatom)
        write(*,*) iatom, in1, dg1, jatom, in2, dg2, num_orb(in1), num_orb(in2)

! interaction spd - spd
        if(num_orb(in1) .eq. 9 .and. num_orb(in2) .eq. 9) then
          write(28,231) 'AAAAAAA' , &
             real(hamk(dg1+1,dg2+1))/ha,real(hamk(dg1+1,dg2+3))/ha,real(hamk(dg1+1,dg2+7))/ha,   &
             real(hamk(dg1+3,dg2+1))/ha,real(hamk(dg1+3,dg2+3))/ha,real(hamk(dg1+2,dg2+2))/ha,   &
             real(hamk(dg1+3,dg2+7))/ha,real(hamk(dg1+4,dg2+8))/ha,real(hamk(dg1+7,dg2+1))/ha,   &
             real(hamk(dg1+7,dg2+3))/ha,real(hamk(dg1+6,dg2+2))/ha,real(hamk(dg1+7,dg2+7))/ha,   &
             real(hamk(dg1+6,dg2+6))/ha,real(hamk(dg1+5,dg2+5))/ha
          write(27,231) 'AAAAAAA' , &
             log10(abs(real(hamk(dg1+1,dg2+1))/ha)),log10(abs(real(hamk(dg1+1,dg2+3))/ha)),   &
             log10(abs(real(hamk(dg1+1,dg2+7))/ha)),log10(abs(real(hamk(dg1+3,dg2+1))/ha)),   &
             log10(abs(real(hamk(dg1+3,dg2+3))/ha)),log10(abs(real(hamk(dg1+2,dg2+2))/ha)),   &
             log10(abs(real(hamk(dg1+3,dg2+7))/ha)),log10(abs(real(hamk(dg1+4,dg2+8))/ha)),   &
             log10(abs(real(hamk(dg1+7,dg2+1))/ha)),log10(abs(real(hamk(dg1+7,dg2+3))/ha)),   &
             log10(abs(real(hamk(dg1+6,dg2+2))/ha)),log10(abs(real(hamk(dg1+7,dg2+7))/ha)),   &
             log10(abs(real(hamk(dg1+6,dg2+6))/ha)),log10(abs(real(hamk(dg1+5,dg2+5))/ha))

! interaction spd - spdd*
        elseif(num_orb(in1) .eq. 9 .and. num_orb(in2) .eq. 14) then

          write(28,231) 'AAAAAAA' , &
             real(hamk(dg1+1,dg2+1))/ha,real(hamk(dg1+1,dg2+3))/ha,real(hamk(dg1+1,dg2+7))/ha,   &
             real(hamk(dg1+3,dg2+1))/ha,real(hamk(dg1+3,dg2+3))/ha,real(hamk(dg1+2,dg2+2))/ha,   &
             real(hamk(dg1+3,dg2+7))/ha,real(hamk(dg1+4,dg2+8))/ha,real(hamk(dg1+7,dg2+1))/ha,   &
             real(hamk(dg1+7,dg2+3))/ha,real(hamk(dg1+6,dg2+2))/ha,real(hamk(dg1+7,dg2+7))/ha,   &
             real(hamk(dg1+6,dg2+6))/ha,real(hamk(dg1+5,dg2+5))/ha,real(hamk(dg1+1,dg2+12))/ha,  &
             real(hamk(dg1+3,dg2+12))/ha,real(hamk(dg1+4,dg2+13))/ha,real(hamk(dg1+7,dg2+12))/ha, &
             real(hamk(dg1+6,dg2+11))/ha,real(hamk(dg1+5,dg2+9))/ha
          write(27,231) 'AAAAAAA' , &
             log10(abs(real(hamk(dg1+1,dg2+1))/ha)),log10(abs(real(hamk(dg1+1,dg2+3))/ha)),   &
             log10(abs(real(hamk(dg1+1,dg2+7))/ha)),log10(abs(real(hamk(dg1+3,dg2+1))/ha)),   &
             log10(abs(real(hamk(dg1+3,dg2+3))/ha)),log10(abs(real(hamk(dg1+2,dg2+2))/ha)),   &
             log10(abs(real(hamk(dg1+3,dg2+7))/ha)),log10(abs(real(hamk(dg1+4,dg2+8))/ha)),   &
             log10(abs(real(hamk(dg1+7,dg2+1))/ha)),log10(abs(real(hamk(dg1+7,dg2+3))/ha)),   &
             log10(abs(real(hamk(dg1+6,dg2+2))/ha)),log10(abs(real(hamk(dg1+7,dg2+7))/ha)),   &
             log10(abs(real(hamk(dg1+6,dg2+6))/ha)),log10(abs(real(hamk(dg1+5,dg2+5))/ha))

! interaction spd - sp
        elseif(num_orb(in1) .eq. 9 .and. num_orb(in2) .eq. 4) then
          write(28,231) 'AAAAAAA' , &
             real(hamk(dg1+1,dg2+1))/ha,real(hamk(dg1+1,dg2+3))/ha,real(hamk(dg1+3,dg2+1))/ha,   &
             real(hamk(dg1+3,dg2+3))/ha,real(hamk(dg1+2,dg2+2))/ha,real(hamk(dg1+7,dg2+1))/ha,   &
             real(hamk(dg1+7,dg2+3))/ha,real(hamk(dg1+6,dg2+2))/ha
          write(27,231) 'AAAAAAA' , &
             log10(abs(real(hamk(dg1+1,dg2+1))/ha)),log10(abs(real(hamk(dg1+1,dg2+3))/ha)),   &
             log10(abs(real(hamk(dg1+3,dg2+1))/ha)),log10(abs(real(hamk(dg1+3,dg2+3))/ha)),   &
             log10(abs(real(hamk(dg1+2,dg2+2))/ha)),log10(abs(real(hamk(dg1+7,dg2+1))/ha)),   &
             log10(abs(real(hamk(dg1+7,dg2+3))/ha)),log10(abs(real(hamk(dg1+6,dg2+2))/ha))

! interaction spd - sps*p*
        elseif(num_orb(in1) .eq. 9 .and. num_orb(in2) .eq. 8) then
          write(28,231) 'AAAAAAA' , &
             real(hamk(dg1+1,dg2+1))/ha,real(hamk(dg1+1,dg2+3))/ha,real(hamk(dg1+3,dg2+1))/ha,   &
             real(hamk(dg1+3,dg2+3))/ha,real(hamk(dg1+2,dg2+2))/ha,real(hamk(dg1+7,dg2+1))/ha,   &
             real(hamk(dg1+7,dg2+3))/ha,real(hamk(dg1+6,dg2+2))/ha,                              &
             real(hamk(dg1+1,dg2+5))/ha,real(hamk(dg1+1,dg2+7))/ha,real(hamk(dg1+3,dg2+5))/ha,   &
             real(hamk(dg1+3,dg2+7))/ha,real(hamk(dg1+2,dg2+6))/ha,real(hamk(dg1+7,dg2+5))/ha,   &
             real(hamk(dg1+7,dg2+7))/ha,real(hamk(dg1+6,dg2+6))/ha
          write(27,231) 'AAAAAAA' , &
             log10(abs(real(hamk(dg1+1,dg2+1))/ha)),log10(abs(real(hamk(dg1+1,dg2+3))/ha)),   &
             log10(abs(real(hamk(dg1+3,dg2+1))/ha)),log10(abs(real(hamk(dg1+3,dg2+3))/ha)),   &
             log10(abs(real(hamk(dg1+2,dg2+2))/ha)),log10(abs(real(hamk(dg1+7,dg2+1))/ha)),   &
             log10(abs(real(hamk(dg1+7,dg2+3))/ha)),log10(abs(real(hamk(dg1+6,dg2+2))/ha)),   &
             log10(abs(real(hamk(dg1+1,dg2+5))/ha)),log10(abs(real(hamk(dg1+1,dg2+7))/ha)),   &
             log10(abs(real(hamk(dg1+3,dg2+5))/ha)),log10(abs(real(hamk(dg1+3,dg2+7))/ha)),   &
             log10(abs(real(hamk(dg1+2,dg2+6))/ha)),log10(abs(real(hamk(dg1+7,dg2+5))/ha)),   &
             log10(abs(real(hamk(dg1+7,dg2+7))/ha)),log10(abs(real(hamk(dg1+6,dg2+6))/ha))

! interaction sp - spd
        elseif(num_orb(in1) .eq. 4 .and. num_orb(in2) .eq. 9) then
          write(28,231) 'AAAAAAA' , &
             real(hamk(dg1+1,dg2+1))/ha,real(hamk(dg1+1,dg2+3))/ha,real(hamk(dg1+1,dg2+7))/ha,   &
             real(hamk(dg1+3,dg2+1))/ha,real(hamk(dg1+3,dg2+3))/ha,real(hamk(dg1+2,dg2+2))/ha,   &
             real(hamk(dg1+3,dg2+7))/ha,real(hamk(dg1+2,dg2+6))/ha
          write(27,231) 'AAAAAAA' , &
             log10(abs(real(hamk(dg1+1,dg2+1))/ha)),log10(abs(real(hamk(dg1+1,dg2+3))/ha)),   &
             log10(abs(real(hamk(dg1+1,dg2+7))/ha)),log10(abs(real(hamk(dg1+3,dg2+1))/ha)),   &
             log10(abs(real(hamk(dg1+3,dg2+3))/ha)),log10(abs(real(hamk(dg1+2,dg2+2))/ha)),   &
             log10(abs(real(hamk(dg1+3,dg2+7))/ha)),log10(abs(real(hamk(dg1+2,dg2+6))/ha))

!interaction sp - sp
        elseif(num_orb(in1) .eq. 4 .and. num_orb(in2) .eq. 4) then
          write(28,231) 'AAAAAAA' , &
             real(hamk(dg1+1,dg2+1))/ha,real(hamk(dg1+1,dg2+3))/ha,   &
             real(hamk(dg1+3,dg2+1))/ha,real(hamk(dg1+3,dg2+3))/ha,real(hamk(dg1+2,dg2+2))/ha
          write(27,231) 'AAAAAAA' , &
             log10(abs(real(hamk(dg1+1,dg2+1))/ha)),log10(abs(real(hamk(dg1+1,dg2+3))/ha)),   &
             log10(abs(real(hamk(dg1+3,dg2+1))/ha)),   &
             log10(abs(real(hamk(dg1+3,dg2+3))/ha)),log10(abs(real(hamk(dg1+2,dg2+2))/ha))
        else
          write(*,*) 'Error en numero de orbitales ****'
!          STOP
        endif

! close data files
        close(27)
        close(28)
! =========================================================================
!    write out interaction_X_Y.dat for transport calculation
! =========================================================================

! open data file
       open (unit = 30, file = 'header.dat', status = 'unknown')
       open (unit = 31, file = 'inter_aux.dat', status = 'unknown')

! create index of interactions for two_center integrals
       num_nonzero = index_max2c(in1,in2)
       iatom = 1
       in1 = imass(iatom)
       dg1 = degelec(iatom)
       jatom = 2
       in2 = imass(jatom)
       dg2 = degelec(jatom)
       allocate ( hlist(num_nonzero))

! write out header
       write (30,*) ' ********************************'
       write (30,*) ' number of steps (integer)'
       write (30,*) ' dist_min    dist_max'
       write (30,*) num_nonzero
       index = 0
       do issh1 = 1 , nssh(in1)
        l1 = lssh(issh1,in1)
        do issh2 = 1, nssh(in2)
         l2 = lssh(issh2,in2)
         do imu = -min(l1,l2), min(l1,l2)
          index = index + 1
          write (30,100) l1, l2, imu, 20.0d0, 1.0d0, 0.5d0, 1, 1, 100.0d0
         end do ! do imu
        end do ! do issh1
       end do ! do issh2

! convert non_zero elements into 1D array
       do index = 1, num_nonzero
        imu = mu(index,in1,in2)
        inu = nu(index,in1,in2)
        hlist (index) = real(hamk(dg1+imu,dg2+inu))/ha
       end do

! write out interactions      
       write (31,200) 'AAAAAAA',(hlist(index), index = 1, num_nonzero)

! close data files
       close (30)
       close (31)

! Format Statements
! ===========================================================================
100    format (3i3,3f14.6,2i2,f12.6)
200    format (a10,50E20.10)
231    format (a10,50E20.10)
    

        return

      end subroutine writeout_hop
