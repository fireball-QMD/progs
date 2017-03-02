! readbind.f90
! Program Description
! ===========================================================================
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine readbind ()

   use transport
   use interactions
   use configuration
   use dimensions
   implicit none

! Argument Declaration and Description
! ===========================================================================
! Input


! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================

   integer :: inu
   integer :: imu
   integer :: i
   integer :: j
   integer :: iatom
   integer :: jatom
   integer :: in1
   integer :: natm
   integer :: iwork
   integer :: index

   logical :: newspec

   
! Procedure
! ===========================================================================

   write (*,*) '  '
   write (*,*) '  '
   write (*,*) ' Welcome to readbind.f! '
   write (*,*) ' Read interaction.optional file '
   write (*,*) ' defining systems for transport calculation '


! open input file  
   open (unit = 17, file = 'interaction.optional', status = 'old')

! Read the basis file
   write (*,100)

!************************************************************
! Sample_1
!************************************************************

! read number of cells 
   read (17,*) sample1%ncell
! read number of atoms in the sample_1
   read (17,*) sample1%natom 
   allocate ( sample1%atom(sample1%natom) )

! read atoms of the sample_1
   read (17,*) iwork
   i = 1
   do inu = 1,iwork
      read (17,*) iatom,jatom
! loop over the interval
      do imu = iatom,jatom
! control check
         if (i .gt. natoms) then 
            write (*,*) 'ERROR, atom ',i,'exceed max. number of atoms.'
            stop
         endif
         sample1%atom(i) = imu
         i = i + 1
      enddo ! do imu
   enddo ! do inu

!************************************************************
! Tip1
!************************************************************

! read information about Tip1
! number of tip's atoms
  read(17,*) iwork
  sample1%natom_tip = iwork
!allocate arrays
  allocate ( sample1%atom_tip (iwork) )
  allocate ( sample1%t2s (iwork) )
  allocate ( sample1%ratom_tip (3,iwork) )
  allocate ( sample1%spec (nspecies))
  allocate ( sample1%ispec (iwork))

! read list of atoms of tip1
  read(17,*) sample1%atom_tip (1:iwork)

! relate the tip atoms with the general list of atoms
  sample1%norb_tip = 0
  do i = 1,sample1%natom_tip
     iatom = sample1%atom_tip(i)
     in1 = imass (iatom)
     sample1%ratom_tip(:,i) = ratom(:,iatom)
     sample1%norb_tip = sample1%norb_tip + num_orb(in1)
  enddo

! calculate number of orbitals of the sample_1  
  sample1%norb = 0
  do i = 1,sample1%natom
     iatom = sample1%atom(i)
     in1 = imass (iatom)
     sample1%norb =  sample1%norb + num_orb (in1)
  enddo

! loop over atoms in tip
  do i = 1,sample1%natom_tip
     iatom =sample1%atom_tip(i)
! loop over atoms in sample1
     do j = 1,sample1%natom
        if (sample1%atom(j) .eq. iatom) then 
           sample1%t2s(i) = j
           exit
        endif
     enddo
  enddo

! determine number of species in the tip
  iatom =sample1%atom_tip(1)
  sample1%spec(1) = imass(iatom)
  index = 1
  do i = 2, sample1%natom_tip
   iatom = sample1%atom_tip(i)
   write (*,*) iatom, imass(iatom)
   newspec = .true.
   do j = 1, index
    if (imass(iatom) .eq. sample1%spec(j)) newspec = .false.  
   enddo ! do j
   write (*,*) 'flag ',newspec
   if (newspec .eqv.  .true.) then
    index = index + 1
    sample1%spec(index) = imass(iatom)
   endif 
  end do ! do i
  sample1%nspec = index
  write (*,*) 'nspec tip1 =',sample1%nspec
  write (*,*) 'id spec tip1:', sample1%spec(1:sample1%nspec)

! map individual tip atoms into the species list of the tip1
  do i = 1, sample1%natom_tip
   iatom = sample1%atom_tip(i)
   in1 = imass(iatom)
   do j = 1,sample1%nspec
    if (in1 .eq. sample1%spec(j)) sample1%ispec(i) = j
   end do ! do j
  end do ! do i


! write informations
  write (*,*) '              Sample_1   '
  write (*,100)
  write (*,201) sample1%ncell
  write (*,200) sample1%natom
  write (*,210) sample1%norb
  write (*,220) 
  write (*,230) sample1%natom_tip
  write (*,240) sample1%norb_tip
  write (*,*) 'Atoms of Sample_1:'
  do i = 1,sample1%natom
     iatom = sample1%atom(i)
     in1 = imass(iatom)
     write (*,250) i,iatom,in1,(ratom(inu,iatom),inu=1,3)
  enddo
  write (*,*) 'Atoms of Tip_1:'
  do i = 1,sample1%natom_tip
     iatom = sample1%atom_tip(i)
     in1 = imass(iatom)
     write (*,260) i,iatom,sample1%t2s(i),in1,(ratom(inu,iatom),inu=1,3)
  enddo
!  endif

! create map of atoms into sample_1 matrix
   natm = sample1%natom
   allocate (degelec1(natm))

   degelec1(1) = 0
   do i = 2,natm
! taking the previous atom in the list
      iatom = sample1%atom(i-1)
      in1 = imass(iatom)
      degelec1(i) = degelec1(i-1) + num_orb(in1)
   enddo ! end iatom

! create map of sample atoms into tip_1 matrix
   natm = sample1%natom_tip
   allocate (pointer1(natm))

   pointer1(1) = 0
   do i = 2,natm
      iatom = sample1%atom_tip(i-1)
      in1 = imass(iatom)
      pointer1(i) = pointer1(i-1) + num_orb(in1)
   enddo ! end iatom


!************************************************************
! Sample2
!************************************************************

! now we consider number of atoms in samples equal number of atoms in 
! whole system, but maybe in future we want to change it, for this reason 
! let's introduce variable sample%natom

! read number of cells 
   read (17,*) sample2%ncell
! read number of atoms in sample_1
   read (17,*) sample2%natom 
   allocate ( sample2%atom(sample2%natom) )

! read number of interval defining atoms of the sample_1
   read (17,*) iwork
   i = 1
   do inu = 1,iwork
      read (17,*) iatom,jatom
! loop over the interval
      do imu = iatom,jatom
! set the elements to off (1)
         if (i .gt. natoms) then 
            write (*,*) 'ERROR, atom ',i,'exceed max. number of atoms.'
            stop
         endif
         sample2%atom(i) = imu
         i = i + 1
      enddo ! do imu
   enddo ! do inu

!************************************************************
! Tip2
!************************************************************

! read information about Tip2
  read(17,*) iwork
  sample2%natom_tip = iwork
  allocate ( sample2%atom_tip (iwork) )
  allocate ( sample2%t2s (iwork) )
  allocate ( sample2%ratom_tip (3,iwork) )
  allocate ( sample2%spec (nspecies))
  allocate ( sample2%ispec (iwork))

! read list of atoms included in tip1
  read(17,*) sample2%atom_tip (1:iwork)

! relate the tip atoms with list of atoms
  sample2%norb_tip = 0
  do i = 1,sample2%natom_tip
     iatom = sample2%atom_tip(i)
     in1 = imass (iatom)
     sample2%ratom_tip(:,i) = ratom(:,iatom)
     sample2%norb_tip = sample2%norb_tip + num_orb(in1)
  enddo
  
! calculate number of orbitals of sample2
  sample2%norb = 0
  do i = 1,sample2%natom
     iatom = sample2%atom(i)
     in1 = imass (iatom)
     sample2%norb =  sample2%norb + num_orb (in1)
  enddo

! loop over atoms in tip
  do i = 1,sample2%natom_tip
     iatom =sample2%atom_tip(i)
! loop over atoms in sample1
     do j = 1,sample2%natom
        if (sample2%atom(j) .eq. iatom) then 
           sample2%t2s(i) = j
           exit
        endif
     enddo
  enddo

! determine number of species in the tip2
  iatom = sample2%atom_tip(1)
  sample2%spec(1) = imass(iatom)
  index = 1
  do i = 2, sample2%natom_tip
   iatom = sample2%atom_tip(i)
   newspec = .true.
   do j = 1, index
    if (imass(iatom) .eq. sample2%spec(j)) newspec = .false.
   enddo ! do j
   if (newspec .eqv. .true.) then
    index = index + 1
    sample2%spec(index) = imass(iatom)
   endif
  end do ! do i
  sample2%nspec = index

! map individual tip atoms into the species list of the tip
  do i = 1, sample2%natom_tip
   iatom = sample2%atom_tip(i)
   in1 = imass (iatom)
   do j = 1,sample2%nspec
    if (in1 .eq. sample2%spec(j)) sample2%ispec(i) = j
   end do ! do j
  end do ! do i

  write (*,*) 'nspec tip2 =', sample2%nspec
  write (*,*) 'ispec tip2:', sample2%spec(1:sample2%nspec)

! write informations
  write (*,*) '              Sample_2   '
  write (*,100)
  write (*,201) sample2%ncell
  write (*,200) sample2%natom
  write (*,210) sample2%norb
  write (*,220) 
  write (*,230) sample2%natom_tip
  write (*,240) sample2%norb_tip
  write (*,*) 'Atoms of Sample_2:'
  do i = 1,sample2%natom
     iatom = sample2%atom(i)
     in1 = imass(iatom)
     write (*,250) i,iatom,in1,(ratom(inu,iatom),inu=1,3)
  enddo
  write (*,*) 'Atoms of Tip_2:'
  do i = 1,sample2%natom_tip
     iatom = sample2%atom_tip(i)
     in1 = imass(iatom)
     write (*,260) i,iatom,sample2%t2s(i),in1,(ratom(inu,iatom),inu=1,3)
  enddo
! endif

! create map of atoms into sample_1 matrix
   natm = sample2%natom
   allocate (degelec2(natm))

   degelec2(1) = 0
   do i = 2,natm
! taking the previous atom in the list
      iatom = sample2%atom(i-1)
      in1 = imass(iatom)
      degelec2(i) = degelec2(i-1) + num_orb(in1)
   enddo ! end iatom

! create map of sample atoms into tip_1 matrix
   natm = sample2%natom_tip
   allocate (pointer2(natm))

   pointer2(1) = 0
   do i = 2,natm
      iatom = sample2%atom_tip(i-1)
      in1 = imass(iatom)
      pointer2(i) = pointer2(i-1) + num_orb(in1)
   enddo ! end iatom

!   write (*,*) 'degelec'
!   do i=1,natoms
!      write (*,*) i,degelec1(i),degelec2(i)
!   enddo

! close input file
   close (unit = 17)

! Format Statements
! ===========================================================================
100     format (2x, 70('='))
201     format (2x, 'ncell = ',i5)
200     format (2x, 'natom = ',i5)
210     format (2x, 'num_orb = ',i5)
220     format (2x, 'Tip segment :')
230     format (4x, 'natom_tip =',i5) 
240     format (4x, 'norb_tip =',i5) 
250     format (2x, i4, i5, i3, 3f12.5)
260     format (2x, i4, i5, i4, i3, 3f12.5)

   return
   
 end subroutine readbind
