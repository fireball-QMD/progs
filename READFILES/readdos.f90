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

! readdos.f90
! Program Description
! ===========================================================================
! This routine reads the information for the dos calculation given in the 
! NEW 'dos.optional' input file. The file will be read if iwrtdos=1.
! dos.optional is something like this
!
!      1.0                       ! scale factor of coord   
!      1        5                ! natom_beg, natom_end for dos calculation
!        151                     ! number of energy steps
!   -10.100000  0.1              ! initial energy for dos calculation
!      0                         ! iwrttip=1 writes the file tip_e_str.inp
!     -0.51     0.01             ! ener_beg, ener_end if iwrttip=1
!      0.1                       ! eta imaginary part in the DOS calculation
! natom_beg and natom_end are the number of the initial and final atoms for the
! dos calculation. So we don't need to do the calculation for all the atoms of 
! the system.
! iwrttip=1 writes the file tip_e_str.inp we need for the STM calculation, the
! energy maximum and minimum limits are given by ener_beg an ener_end. 
  
! ===========================================================================
! Code rewritten by:
! Cesar Gonzalez Pascual
! Departamento de Fisica Teorica de la Materia Condensada
! Facultad de Ciencias. Universidad Autonoma de Madrid
! E28049 Madrid SPAIN
! FAX +34-91 3974950
! Office Telephone +34-91 3978648
!===========================================================================
!
! Program Declaration
!===========================================================================
        subroutine readdos ( )
        use dimensions
        use interactions
        use module_dos
        implicit none

! Argument Declaration and Description
! ========================================== 

! Local Variable Declaration and Description
! ==========================================================================
       integer iatom
       integer in1

! Procedure
! ===========================================================================
! Initialize some things
       norb_act = 0

! Opening the file
       open (unit = 121, file = 'dos.optional')

! Read dos.optional
       read (121,*) lattice                    ! ratio of lattice parameter 
! Read the number of atoms (the first ant the final atoms)  
! We want to have the dos
       read(121,*) natom_beg,natom_end   
       do iatom = natom_beg, natom_end
        in1 = imass(iatom)
        norb_act = norb_act + num_orb(in1)
       end do

! Read the energy stuff
       read(121,*) nener                      ! number of energies
       read(121,*) ener_beg, ener_step        ! first energy and step for DOS
       read(121,*) iwrttip                    ! 1/0 yes/no write tip_e_str.inp
       read(121,*) ener_min, ener_max         ! minimun and maximum energies 
       read(121,*) eta                        ! imaginary part green function 

       return
       end subroutine readdos
