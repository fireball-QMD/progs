! experiments.f
! We put he experimental stuff here like levels, etc.
! ===============================================================
        subroutine experiments(natoms,imass,delep,delover)
        use dimensions 
        use interactions
        implicit none
        logical readl
        integer mu,nu,i,natoms,in1
        real delep(numorb_max,natoms),delover(numorb_max,natoms)
        integer, intent(in) :: imass (natoms)
 
 
! The levels experiment.
! LEVELS allows you to alter the energy levels of the on-sites.
! Thus you take levels out of resonance or change the relative ordering
! of the levels. Please dont do this except for testing.
! Also, in some cases such as alkali metals, you may wish to "remove"
! orbitals from the problem. For example, p-sates of alkali metals.
! Then you also need to reduce the overlap of those orbitals with their
! neighbors.
! Below delep changes the on site enrgy levels, and delover alter the
! over lap of that orbital with all neighbors. (Note delover does not
! reduce the overlap of that orbital with itself, because of normalization).
! Note: delep=0 means no shift of onsites, and delover=1 means keep
! the same size overlaps.
! Warning I have not set up delover for forces'
!
! Initialize first:
        do i=1,natoms
         in1=imass(i)
         do mu =1,num_orb(in1)
          delep(mu,i)=0.0
          delover(mu,i)=1.0
         enddo
        enddo
        inquire(file="LEVELS",exist=readl)
        if(readl)then
         open(unit=16,file="LEVELS",status="old")
         write(*,*) ' ******************************'
         write(*,*) ' ******************************'
         write(*,*) 'Found file LEVELS: shifting diagonal energies....'
         write(*,*) 'Very dangerous if you dont know what youre doing'
         write(*,*) ' I have NOT set delover for forces yet;'
         write(*,*) ' I have NOT set delover for forces yet;'
         write(*,*) ' I have NOT set delover for forces yet;'
         write(*,*) ' I have NOT set delover for forces yet;'
         write(*,*) ' ******************************'
         write(*,*) ' ******************************'
         do i=1,natoms
          in1=imass(i)
          read(16,*)(delep(mu,i),mu=1,num_orb(in1)),   &
     &    (delover(nu,i),nu=1,num_orb(in1))
          write(*,814)(delep(mu,i),mu=1,num_orb(in1)),   &
     &    (delover(nu,i),nu=1,num_orb(in1))
814       format(' E=',4f14.5,/,' S=',4f10.5)
         enddo
         close(unit=16)
        endif
        return
        end
 
