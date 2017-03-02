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

! assemble_scissor.f90
! Program Description
! ===========================================================================
!       This routine assembles the Koopman correction potential to correctly
!	calculate the molecular gap
!
! ===========================================================================
! Code written by:
! Y.J. Dappe & Enrique Abad Gonzalez
! Departamento de Fisica Teorica de la Materia Condensada
! Facultad de Ciencias C-V
! Universidad Autonoma de Madrid
! Campus de Cantoblanco
! FAX +34-91-497-4950
! Office telephone  +34-91-497-8648
! email: enrique.abad@uam.es
!
! ===========================================================================
        subroutine assemble_scissor (natoms)

        use hartree_fock
        use interactions
        use neighbor_map
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms

! Output

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer iatom, jatom, ineigh, mbeta
        integer imu, inu, ieigen
        integer jmu, jnu
        integer irho, isigma, jrho, jsigma, katom, latom, in3, in4
        integer ineigh2, jneigh, nbeta, lbeta
        integer in1, in2
        integer inum, iremainder
        integer num_mol, mol_orbitals, nniveles, numkpts, icount
	real, dimension(:),allocatable :: deltae
        real, dimension(:,:), allocatable :: bebenkare, bebenkaim

! Procedure
! ===========================================================================
! Initialize


!        allocate(hs_mat(numorb_max,numorb_max,neigh_max,natoms))
        hs_mat = 0.0d0
        natomhf_beg = 1
        natomhf_end = 60

        num_mol = 1 ! the number of molecules we introduce in our calculations,
                    ! it'll be put in readfiles someday
!        mol_orbitals = degelec(natomhf_end) + num_orb(natomhf_end) -degelec(natomhf_beg) ! number of
        ! molecular orbitals

        mol_orbitals = 240

!Read the Lowdin coeficient coming from Koopman subroutine
        open(unit=15, file='./koopman_shift', status='old')
        open(unit=416,file='./cdcoeffs.dat', status='old')

        allocate(deltae(mol_orbitals))
        do ieigen = 1, mol_orbitals
          read (15,*) deltae(ieigen)
        end do
! Leemos los numeros de niveles y kpoints

read (416,*) nniveles
if ( mol_orbitals .ne. nniveles ) then
   write (*,*) "WARNING! The number of orbitals you've in your molecule"
   write (*,*) "does not fit with the number of orbitals you've stored!"
   STOP
end if
read  (416,*)
read  (416,*)
read  (416,*)
read  (416,*) numkpts
if ( numkpts .ne. 1 ) then
   write (*,*) "WARNING! You've stored coefficients not using ONLY"
   write (*,*) "gamma point!"
end if

! Allocateamos cosas

allocate(bebenkare(mol_orbitals,mol_orbitals))
allocate(bebenkaim(mol_orbitals,mol_orbitals))

! Leemos los coeficientes

  read  (416,*)
  do imu = 1, mol_orbitals
    read  (416,*)
    inum = int(mol_orbitals/4) ! We want int/int division
    iremainder = mol_orbitals - (inum*4)
    do inu = 1, mol_orbitals - iremainder, 4
      read  (416,102)                                               &
      &       bebenkare(inu  ,imu), bebenkaim(inu  ,imu),           &
      &       bebenkare(inu+1,imu), bebenkaim(inu+1,imu),           &
      &       bebenkare(inu+2,imu), bebenkaim(inu+2,imu),           &
      &       bebenkare(inu+3,imu), bebenkaim(inu+3,imu)
    end do
    if ( iremainder .eq. 1 ) then
             read  (416,103)                                        &
     &        bebenkare(mol_orbitals,imu),                             &
     &        bebenkaim(mol_orbitals,imu)

    else if ( iremainder .eq. 2 ) then
      read  (416,104)                                                &
      &        bebenkare(mol_orbitals-1,imu),                           &
      &        bebenkaim(mol_orbitals-1,imu),                           &
      &        bebenkare(mol_orbitals  ,imu),                           &
      &        bebenkaim(mol_orbitals  ,imu)
    else if ( iremainder .eq. 3 ) then
             read  (416,105)                                        &
     &        bebenkare(mol_orbitals-2,imu),                           &
     &        bebenkaim(mol_orbitals-2,imu),                           &
     &        bebenkare(mol_orbitals-1,imu),                           &
     &        bebenkaim(mol_orbitals-1,imu),                           &
     &        bebenkare(mol_orbitals  ,imu),                           &
     &        bebenkaim(mol_orbitals  ,imu)
    end if
  end do



        do iatom = natomhf_beg, natomhf_end ! Loop over the atoms in the central cell.
          in1 = imass(iatom)
!          do ineigh = 1, neighn(iatom)       ! NO!! The scissor hamiltonian has terms between non-neighbors
          do jatom = natomhf_beg, natomhf_end
            in2 = imass(jatom)
! I THINK IT IS NOT NECESSARY NOW CHECK!
!            mbeta = neigh_b(ineigh,iatom)
!            icount = 0
!            if(iatom.le.natomhf_end.and.iatom.ge.natomhf_beg) icount = icount+1
!            if(jatom.le.natomhf_end.and.jatom.ge.natomhf_beg) icount = icount+1
!            if(mbeta.eq.0) icount = icount+1
!            if (icount.eq.3) then
            do ineigh2 = 1, neighn(iatom)
              katom = neigh_j(ineigh2,iatom)
              nbeta = neigh_b(ineigh2,iatom)
              in3 = imass(katom)
              do jneigh = 1, neighn(jatom)
                latom = neigh_j(jneigh,jatom)
		lbeta = neigh_b(jneigh,jatom)
                in4 = imass(latom)
                icount = 0
                if(katom.le.natomhf_end.and.katom.ge.natomhf_beg) icount = icount+1
                if(latom.le.natomhf_end.and.latom.ge.natomhf_beg) icount = icount+1
		if (nbeta.eq.0.and.lbeta.eq.0) icount = icount+1
                if (icount .eq.3) then
                do imu = 1, num_orb(in1)
                   jmu = imu + degelec(iatom) !- degelec(natomhf_beg)  ! esto solo vale si natomhf_beg = 1 y
                   ! num_molec = 1, lo mismo con jnu
                   ! la solucion seria jmu = imu + degelec(iatom-natomhf_beg+1)?
                   do inu = 1, num_orb(in2)
                     jnu = inu + degelec(jatom) !- degelec(natomhf_beg)  Tal y como he reescrito hs_mat jnu es esto
		     do irho = 1, num_orb(in3)
		       jrho = irho + degelec(katom) - degelec(natomhf_beg)
		       do isigma = 1, num_orb(in4)
		         jsigma = isigma + degelec(latom) - degelec(natomhf_beg)
                         do ieigen = 1, mol_orbitals
                     hs_mat(jmu,jnu) = hs_mat(jmu,jnu) +                             &
    &  deltae(ieigen)*bebenkare(jrho,ieigen)*bebenkare(jsigma,ieigen)*               &
    &  s_mat(imu,irho,ineigh2,iatom)*s_mat(inu,isigma,jneigh,jatom)

	                 end do
		       end do
		     end do
		   end do
		end do
		end if
              end do
            end do
!            end if
          end do
        end do ! End loop over iatom.
        close(unit=15)
        close(unit=416)





! Deallocate arrays
! ===========================================================================

! Format Statements
! ===========================================================================
100     format (3f11.6)
101     format (4f11.5)
102     format (4(1x, '(', f14.11, ',', f14.11, ')'))
103     format (1(1x, '(', f14.11, ',', f14.11, ')'))
104     format (2(1x, '(', f14.11, ',', f14.11, ')'))
105     format (3(1x, '(', f14.11, ',', f14.11, ')'))

        return
        end subroutine assemble_scissor
