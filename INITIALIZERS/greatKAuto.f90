! ====================================================================
! Otto  F. Sankey
! otto.sankey@asu.edu
!
! j keith
! jbrkeith@gmail.com
! 
! Latest Update June 13, 1998.
! Revised Aug. 2, 1998.
! Revised Nov, 2005
! ====================================================================
! 
! Let me give you an example:
! Suppose you want forces in a diamond lattice.
! The special k-points you will get normally
! are (using 4,4,4 for the iq values)
! (1,1,1) (3,1,1) (1,3,1) (1,1,3) in certain units.
! But where is (-1,-1,1) etc. Well, they are not
! there since the special kpoint algorithm does
! not insert symmetry related kpoints. It only
! satisfies the sum of lattice vectros condition.
! The points give above certainly satisfy this
! condition. However, if you now use these k-points to
! calculate forces, you will not get zero force
! on an atom even at a Td position. So to fix this
! we tell this program  that we want our k points
! rotated to have a certain symmetry.
! using td.sym, the 4 kpoints above become:
!   -0.2892   -0.2892   -0.2892
!    0.2892   -0.8678    0.2892
!   -0.8678    0.2892    0.2892
!   -0.2892   -0.2892    0.8678
!   -0.2892    0.2892    0.2892
!   -0.8678   -0.2892   -0.2892
!    0.2892   -0.2892    0.8678
!   -0.2892   -0.8678    0.2892
!    0.2892   -0.2892    0.2892
!   -0.2892    0.2892    0.8678
!   -0.2892   -0.8678   -0.2892
!   -0.8678   -0.2892    0.2892
!    0.2892    0.2892   -0.2892
!   -0.2892   -0.2892   -0.8678
!   -0.2892    0.8678    0.2892
!   -0.8678    0.2892   -0.2892
!
! So this program will ask you for 2 symmetry files.
! The first one rotates the k points to generates
! a symmetry set, and the second one reduces the set by
! the given symmetry operators.
!
! So if you want a symmetry set because you
! are doing forces, then
! (First set) Insert the symmetry you want (e.g. td.sym)
! (2nd   set) Insert identity (id.sym)
! 
! So if you want energies and are tyring to reduce the
! set, then
! (First set) Insert identity (id.sym)
! (2nd   set) Insert symmetry to reduce set (eg. td.sym)
!
! Program Declaration
! ===========================================================================
        subroutine greatKAuto(lvsfile,basisfile,mpmesh,iquench, &
     &     ireducekpts,numb)
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: mpmesh(3), iquench, ireducekpts
        character(len=40), intent(in) :: lvsfile, basisfile
! Output
        integer, intent(out) :: numb

! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: nkmax=100000
        real accuracy, pi
        data accuracy/5.0d-4/
        data pi/3.1415926/
 
! Local Variable Declaration and Description
! ===========================================================================
        integer natoms,iatomtype(:),npgops,i1max,i2max,i3max, &
     &   ncase1,ncase2,ncase3,iplan2,i,j,n1,n2,n3,i1,i2,i3,iplan,L, &
     &   ix,ik,numirr,ikillminus,numbnew,ik1,ikn,icont,ik2 

        real atom_position(:,:),a1(3),a2(3),a3(3),  &
     &   acrossb(3),bcrossc(3),ccrossa(3),rlv1(3),rlv2(3),rlv3(3),  &
     &   spkpt(3),sk(3,nkmax),g(3),skirred(3,nkmax),wt(nkmax),  &
     &   weight(nkmax),wtnew(nkmax),ikill(nkmax),  &
     &   primitive_lattice(3,3),temp2(:,:),pgopscart(3,3,96),  &
     &   opvector(3,48),gmag,denom,rlv1mag,rlv2mag,  &
     &   rlv3mag,s1,s2,s3,f1,f2,f3,size1,outornot

!        common /imaxes/i1max,i2max,i3max
!        common /lattice/a1,a2,a3,rlv1,rlv2,rlv3
        allocatable iatomtype,atom_position,temp2

! Procedure
! ===========================================================================
        ncase1=0
        ncase2=0
        ncase3=0
!	write(*,*)' Very general. Also gives special k '
!	write(*,*)' points that have a certain symmetry.'
!	write(*,*)' Let me give you an example: '
!	write(*,*)' Suppose you want forces in a diamond lattice. '
!	write(*,*)' The special k-points you will get normally '
!	write(*,*)' are (using 4,4,4 for the iq values) '
!	write(*,*)' (1,1,1) (3,1,1) (1,3,1) (1,1,3) in certain units.'
!	write(*,*)' But where is (-1,-1,1) etc. Well, they are not '
!	write(*,*)' there since the special kpoint algorithm does '
!	write(*,*)' not insert symmetry related kpoints.'
!	write(*,*)' After all, it uses only the lattice vectors --'
!	write(*,*)' how would it know about the point group'
!	write(*,*)' The special k points is constructed to '
!	write(*,*)' satisfy the sum of lattice vectors condition. '
!	write(*,*)' The 4 points give above certainly satisfy this '
!	write(*,*)' condition. However, if you now use these k-points'
!	write(*,*)' to calculate forces, you will not get zero force '
!	write(*,*)' on an atom even at a Td position. So to fix this '
!	write(*,*)' we tell thias program  that we want our k points '
!	write(*,*)' rotated to have a certain symmetry. '
!	write(*,*)' using td.sym, the 4 kpoints above become: '
!	write(*,*)'   -0.2892   -0.2892   -0.2892 '
!	write(*,*)'    0.2892   -0.8678    0.2892 '
!	write(*,*)'   -0.8678    0.2892    0.2892 '
!	write(*,*)'   -0.2892   -0.2892    0.8678 '
!	write(*,*)'   -0.2892    0.2892    0.2892 '
!	write(*,*)'   -0.8678   -0.2892   -0.2892 '
!	write(*,*)'    0.2892   -0.2892    0.8678 '
!	write(*,*)'   -0.2892   -0.8678    0.2892 '
!	write(*,*)'    0.2892   -0.2892    0.2892 '
!	write(*,*)'   -0.2892    0.2892    0.8678 '
!	write(*,*)'   -0.2892   -0.8678   -0.2892 '
!	write(*,*)'   -0.8678   -0.2892    0.2892 '
!	write(*,*)'    0.2892    0.2892   -0.2892 '
!	write(*,*)'   -0.2892   -0.2892   -0.8678 '
!	write(*,*)'   -0.2892    0.8678    0.2892 ' 
!	write(*,*)'   -0.8678    0.2892   -0.2892 '
!	write(*,*)'  '
!	write(*,*)' This program will compute the symmetry of your basis'
!        write(*,*)' file and use that first to possibly rotate the '
!        write(*,*)' kpoint grid you specify and then possibly to reduce'
!	write(*,*)' the kpoint set by these same symmetry operators. '
!	write(*,*)' '
!	write(*,*)' So if you want a symmetry set because you '
!	write(*,*)' are doing forces, then '
!	write(*,*)' (First set) Insert 1 (true--use symmetry to rotate)'
!	write(*,*)' (2nd   set) Insert 0 (false--don"t use symmetry to'
!	write(*,*)' reduce)'
!        write(*,*)
!	write(*,*)' If you want energies and are tyring to reduce'
!	write(*,*)' the set, then '
!	write(*,*)' (First set) Insert 1 (true--use symmetry to rotate)'
!	write(*,*)' (2nd   set) Insert 1 (true--use symmetry to reduce)'
!	write(*,*)'  '
!	write(*,*)' The greatK.com file for a "standard" diamond kps file'
!	write(*,*)' follows:'
!	write(*,*)'cat << EOF | ./greatKAuto.x'
!	write(*,*)'diamond.lvs'
!        write(*,*)'diamond.bas'
!	write(*,*)'2,2,2             n1,n2,n3 Monkhorst-Pack integers.'
!	write(*,*)'1                 (true--use symmetry to rotate)'
!	write(*,*)'1                 1 for inversion, 0 no inversion'
!	write(*,*)'1                 (true--use symmetry to reduce)'
!	write(*,*)'1                 1 for inversion'
!	write(*,*)'1                 Remove k"s which are minus a - k.'
!	write(*,*)'EOF'

	iplan=2
	open (unit=30,file=lvsfile,status='old')
        read(30,*)primitive_lattice
        close(30)
        a1(:)=primitive_lattice(1,:)
        a2(:)=primitive_lattice(2,:)
        a3(:)=primitive_lattice(3,:)
! ====================================================================
        call gkcross(a2,a3,bcrossc)
        call gkcross(a3,a1,ccrossa)
        call gkcross(a1,a2,acrossb)
        denom=a1(1)*bcrossc(1)+a1(2)*bcrossc(2)+a1(3)*bcrossc(3)
        do i=1,3
         rlv1(i)=2.0d0*pi*bcrossc(i)/denom
         rlv2(i)=2.0d0*pi*ccrossa(i)/denom
         rlv3(i)=2.0d0*pi*acrossb(i)/denom
        end do
        rlv1mag=sqrt(rlv1(1)*rlv1(1)+rlv1(2)*rlv1(2)+rlv1(3)*rlv1(3))
        rlv2mag=sqrt(rlv2(1)*rlv2(1)+rlv2(2)*rlv2(2)+rlv2(3)*rlv2(3))
        rlv3mag=sqrt(rlv3(1)*rlv3(1)+rlv3(2)*rlv3(2)+rlv3(3)*rlv3(3))
!	    write(*,*)' The reciprocal lattice vector are:'
!        write(*,12)rlv1(1),rlv1(2),rlv1(3),rlv1mag
!        write(*,13)rlv2(1),rlv2(2),rlv2(3),rlv2mag
!        write(*,14)rlv3(1),rlv3(2),rlv3(3),rlv3mag
!12      format(' rlv1 = ',f10.5,1x,f10.5,1x,f10.5,3x,' Mag.=',f10.5)
!13      format(' rlv2 = ',f10.5,1x,f10.5,1x,f10.5,3x,' Mag.=',f10.5)
!14      format(' rlv3 = ',f10.5,1x,f10.5,1x,f10.5,3x,' Mag.=',f10.5)
! ====================================================================
!        write(*,*)' Insert filename containing atomic positions'
!	write(*,*)' e.g. junk.bas'
!        read(*,'(a30)')basisfile
!	write(*,'(a30)')basisfile
	open (unit=30,file=basisfile,status='old')
    read(30,*)natoms
    allocate(iatomtype(natoms),atom_position(3,natoms))
    read(30,*)(iatomtype(i),(atom_position(j,i),j=1,3),i=1,natoms)
!        write(*,*)(iatomtype(i),(atom_position(j,i),j=1,3),i=1,natoms)
    close(30)

!       use symmetry of atoms to reduce number of kpoints if doing 
!       optimization or if user requests it (ireducekpts)
    if(ireducekpts.eq.1) then
       allocate(temp2(3,natoms))
!         print*,'accuracy',accuracy
!         print*,'primitive_lattice',primitive_lattice
!         print*,'natoms',natoms
!         print*,'iatomtype',iatomtype
!         print*,'atom_position',atom_position
!         print*,'i',i
!         print*,'temp2',temp2
!         print*,'npgops',npgops
!         print*,'pgopscart',pgopscart
!         print*,opvector
         call find_symmetry3(accuracy,primitive_lattice,natoms, & 
     &     iatomtype,atom_position,i,temp2,npgops,pgopscart,    &
     &     opvector)
         write(6,'(a,i2,a)')'Found ',npgops,' symmetry operators'
        else
         npgops=1
         pgopscart(1,1,1)=1
         pgopscart(1,2,1)=0
         pgopscart(1,3,1)=0
         pgopscart(2,1,1)=0
         pgopscart(2,2,1)=1
         pgopscart(2,3,1)=0
         pgopscart(3,1,1)=0
         pgopscart(3,2,1)=0
         pgopscart(3,3,1)=1
         print*,'not using symmetry to reduce kpoints'
        endif

!	write(*,*)'  '
!	write(*,*)' Now insert the three Monkhorst-Pack integers:'
!	write(*,*)' (e.g. 2,2,2). Must be positive integers, 1,2,3,...'
!        read(*,*)n1,n2,n3
    write(*,*)'The Monkhorst-Pack mesh that automatic has chosen are'
        n1=mpmesh(1)
        n2=mpmesh(2)
        n3=mpmesh(3)
        write(*,*)n1,n2,n3
	if(n1.le.0) stop ' bad n1 (Must be positive integers, 1,2,3,...)'
	if(n2.le.0) stop ' bad n2 (Must be positive integers, 1,2,3,...)'
	if(n3.le.0) stop ' bad n3 (Must be positive integers, 1,2,3,...)'

! ====================================================================
! Fractions = m/n + shifts.
! Shifts:
	s1=1./2.+ 1./(2.*n1)
	s2=1./2.+ 1./(2.*n2)
	s3=1./2.+ 1./(2.*n3)
!	write(*,*)' '
!	write(*,*)' We write each k vector as f1*g1+f2*g2+f3*g3'
!	write(*,*)' The fractions f1, f2, and f3 are:'
!	write(*,*)' s1,s2,s3=',s1,s2,s3
	numb=0
	do 16 i1=1,n1
	f1=float(i1)/float(n1) - s1
	do 17 i2=1,n2
        f2=float(i2)/float(n2) - s2
	do 18 i3=1,n3
        f3=float(i3)/float(n3) - s3
!	write(*,799)f1,f2,f3
799	format(' f1, f2, f3 = ',3f9.4)
        do 19 L=1,3
        spkpt(L)=rlv1(L)*f1+rlv2(L)*f2+rlv3(L)*f3
19      continue
        numb=numb+1
        if(numb.gt.nkmax)then
        write(*,*)' numb too large !!!! numb=',numb
        write(*,*)' Dimension is only',nkmax
        write(*,*)' Please redimension!!!!'
        stop
        end if
        do 57 ix=1,3
57      sk(ix,numb)=spkpt(ix)
18	continue
17	continue
16	continue

!        do 16 i=imin,imax, 2
!        frac1=dfloat(i)/dfloat(n1)
!        do 17 j=jmin,jmax,2
!        frac2=dfloat(j)/dfloat(n2)
!        do 18 k=kmin,kmax,2
!        frac3=dfloat(k)/dfloat(n3)
!        do 19 l=1,3
!        spkpt(l)=rlv1(l)*frac1+rlv2(l)*frac2+rlv3(l)*frac3
!19      continue
!	numb=numb+1
!
!	if(numb.gt.nkmax)then
!	write(*,*)' numb too large !!!! numb=',numb
!	write(*,*)' Dimension is only',nkmax
!	write(*,*)' Please redimension!!!!'
!	stop
!	end if
!
!	do 57 ix=1,3
!57	sk(ix,numb)=spkpt(ix)
!c
!18      continue
!17      continue
!16      continue
! ====================================================================
!	write(*,*)'  '
!	write(*,*)' We found numb=',numb,' kpoints.'
!	write(*,*)'  '
!	open (unit=31,file='greatK.kpts',status='unknown')
!	write(31,*)numb
!	write(*,*)' They are (in inverse Angstrom units):'
	do 41 ik=1,numb
	weight(ik)=1.d0/dfloat(numb)
!        write(*,20)sk(1,ik),sk(2,ik),sk(3,ik),weight(ik)
!        write(31,20)sk(1,ik),sk(2,ik),sk(3,ik),weight(ik)
20      format(1x,f13.8,1x,f13.8,1x,f13.8,8x,f13.8)
41	continue
!	close(unit=31)
!        write(*,*)'  '
!        write(*,*)' We write them out to greatK.kpts'
! ====================================================================
!!$	write(*,*)'  '
!!$	write(*,*)' Later we call shells to test which shells give zero for'
!!$	write(*,*)' these special k points.'
!!$	write(*,*)' Insert i1max,i2max,i3max'
!!$	write(*,*)' e.g.    3  ,  3 , 3'
!!$	write(*,*)' These create lattice vectors from'
!!$	write(*,*)' -i1max to +i1max and similarly for 2 and 3'
!!$	write(*,*)' Recall: L=i1*a1+i2*a2+i3*a3'
!!$	write(*,*)' We read these from greatK.input.'
!!$	write(*,*)' Sorry: Now if fix them to 3,3,3'
!	read(17,*)i1max,i2max,i3max
        i1max=3
        i2max=3
        i3max=3
!        write(*,*)i1max,i2max,i3max
!	call shells(a1,a2,a3,sk,numb,weight,i1max,i2max,i3max)
! ====================================================================
! We call kgroup, which rotates them to get a fully symmetric set,
! then weeds out repeats, and optionally keeps only those in the
! irreducible wedge.
!!$	write(*,*)'  '
!!$	write(*,*)' In a moment, you will be asked to supply a'
!!$	write(*,*)' symmetry file. Insert the symmetry file name'
!!$	write(*,*)' *** which will rotate the k vectors to have ***'
!!$	write(*,*)' *** a certain symmetry. This operator will  ***'
!!$	write(*,*)' *** genrally add k vectors to the list.     ***'
!!$	write(*,*)'  '

	call kgroup(sk,weight,numb,skirred,wt,numirr,pgopscart,npgops, &
     &    a1,a2,a3,rlv1,rlv2,rlv3,i1max,i2max,i3max)
!
!	write(*,*)'  '
!	write(*,*)' Back in kspecial... numb=',numb
!	write(*,*)' Now test the k vectors with L shells'
!	call shells(a1,a2,a3,sk,numb,weight,i1max,i2max,i3max)
! Now call kgroupagain to weed out k vectors related by symmetry according
! to ANOTHER group.
!	write(*,*)'  '
!	write(*,*)' Back in kspecial... numb=',numb
!
!!$        write(*,*)'  '
!!$        write(*,*)' In a moment, you will be asked to supply a'
!!$        write(*,*)' symmetry file. Insert the symmetry file name'
!!$        write(*,*)' *** which will rotate the k vectors and     ***'
!!$        write(*,*)' *** eliminate those related by symmetry.    ***'
!!$	write(*,*)' *** This operator will genrally subtract k  ***'
!!$	write(*,*)' *** vectors from the list.                  ***'
!!$        write(*,*)'  '
        call kgroupagain(sk,weight,numb,skirred,wt,numirr,pgopscart, &
     &   npgops,a1,a2,a3,rlv1,rlv2,rlv3,i1max,i2max,i3max)
!	write(*,*)' Now test the k vectors with L shells'
!	call shells(a1,a2,a3,sk,numb,weight,i1max,i2max,i3max)
!
!        write(*,*)' '
!        write(*,*)' NEW ADDITION OFS'
!        write(*,*)' Do you wish to remove those k"s which are minus'
!        write(*,*)' another k? Insert 1 for YES!'
!        read(*,*)ikillminus
!        write(*,*)ikillminus
! No need for ikill, since we do effectively the same thing
! with the group theory file when we ask for inversion.
	ikillminus=0
        if(ikillminus.ne.1)go to 459
! OK we check for k and -k
        do ik=1,numb
        ikill(ik)=0
        wtnew(ik)=weight(ik)
        end do
!
        numbnew=0
        do ik1=1,numb
        if(ikill(ik1).eq.0)then
        numbnew=numbnew+1
        do 511 ik2=ik1+1,numb
        size1=(sk(1,ik1)+sk(1,ik2))**2+(sk(2,ik1)+sk(2,ik2))**2+  &
     &  (sk(3,ik1)+sk(3,ik2))**2
        if(size1.lt.0.00001)then
        ikill(ik2)=1
        wtnew(ik1)=wtnew(ik1)+weight(ik2)
        end if
511     continue
        end if
        end do
        ikn=0
        do 455 ik=1,numb
        if(ikill(ik).eq.1)go to 455
        ikn=ikn+1
        sk(1,ikn)=sk(1,ik)
        sk(2,ikn)=sk(2,ik)
        sk(3,ikn)=sk(3,ik)
        weight(ikn)=wtnew(ik)
455     continue
        if(ikn.ne.numbnew)then
	write(*,*)' ikn=',ikn
	write(*,*)' numbnew=',numbnew
	stop ' ikn.ne.numbnew'
	end if
        numb=numbnew
!	write(*,*)' '
!        write(*,*)' The new number of k-points is',numb
!        write(*,*)'  '
!        write(*,*)' They are ....'
!        do 6441 ik=1,numb
!        write(*,20)sk(1,ik),sk(2,ik),sk(3,ik),weight(ik)
!6441     continue
!        write(*,*)' We found numb=',numb,' kpoints.'
!        write(*,*)'  '
!
459     continue
!
!	write(*,*)'  '
!	write(*,*)' Now test the k vectors with L shells'
!	call shells(a1,a2,a3,sk,numb,weight,i1max,i2max,i3max)
! ====================================================================
! ====================================================================
!        write(*,*)'  '
!        write(*,*)' We found numb=',numb,' kpoints.'
!        write(*,*)' At this point they are in the reciprocal unit'
!        write(*,*)' cell, but not in the first BZ'
!        write(*,*)'  '
!        write(*,*)' We write them out to temp.kpts'
!        write(*,*)'  '
!        open (unit=31,file='greatK.kpts',status='unknown')
!        write(31,*)numb
!        write(*,*)' They are ....'
!        do 1441 ik=1,numb
!        weight(ik)=1.d0/dfloat(numb)
!        write(*,20)sk(1,ik),sk(2,ik),sk(3,ik),weight(ik)
!        write(31,20)sk(1,ik),sk(2,ik),sk(3,ik),weight(ik)
!1441      continue
!        write(*,*)' We found numb=',numb,' kpoints.'
!        close(unit=31)
! ====================================================================
911	continue
!	write(*,*)' =================================================='
!	write(*,*)'  '
!	write(*,*)' We now shift the above k-points so that'
!	write(*,*)' they are inside the first BZ.'
!	write(*,*)'  '
!	write(*,*)' Continue y/n=1/0'
	icont=1
!	read(*,*)icont
	if(icont.eq.1)go to 99
	write(*,*)' Sorry I must stop.'
	stop
99	continue
!	write(*,*)' OK -- I"m going on.'
!	open(unit=44,file='junk.dat',status='unknown')
! Now check to see if in the BZ
!	write(*,*)' Now check to see if in the BZ...'
! Loop over all k vectors
	do 1000 ik=1,numb
! Loop over a large number of g's. We have surely overdone it to go from
! -3 to 3.
	do 1 i1=-3,3
	do 2 i2=-3,3
	do 3 i3=-3,3
! Compute g.
	do ix=1,3
	g(ix)=i1*rlv1(ix)+i2*rlv2(ix)+i3*rlv3(ix)
	end do
	gmag=g(1)**2+g(2)**2+g(3)**2
	gmag=sqrt(gmag)
! skip if g=0
	if(gmag.lt.1.d-08)go to 100
! Now k is inside the plane of 1/2 g if k * g /gmag < gmag/2

	outornot=sk(1,ik)*g(1)+sk(2,ik)*g(2)+sk(3,ik)*g(3)
	outornot=outornot/gmag
	if(outornot.gt.gmag/2.)go to 55
! OK. We are inside gmag/2
	go to 100
! Nope. We are outside gmag/2. Shift k inside by subtracting g.
55	continue
	do ix=1,3
	sk(ix,ik)=sk(ix,ik)-g(ix)
	end do
!	write(*,122)sk(1,ik),sk(2,ik),sk(3,ik)
122	format(' k=',3f9.4)
!
100	continue
3	continue
2	continue
1	continue
! Now let's check to make sure we did OK.
! ---
! Loop over a large number of g's. We have surely overdone it to go from
! -3 to 3.
        do 551 i1=-3,3
        do 552 i2=-3,3
        do 553 i3=-3,3
! Compute g.
        do ix=1,3
        g(ix)=i1*rlv1(ix)+i2*rlv2(ix)+i3*rlv3(ix)
        end do
        gmag=g(1)**2+g(2)**2+g(3)**2
        gmag=sqrt(gmag)
! skip if g=0
        if(gmag.lt.1.d-08)go to 55100
! Now k is inside the plane of 1/2 g if k * g /gmag < gmag/2

        outornot=sk(1,ik)*g(1)+sk(2,ik)*g(2)+sk(3,ik)*g(3)
        outornot=outornot/gmag
        if(outornot.gt.gmag/2.)go to 5555
! OK. We are inside gmag/2
        go to 55100
! Nope. We are outside gmag/2. Shift k inside by subtracting g.
5555      continue
	write(*,*)' Oppps -- didnt work'
!
55100     continue
553	continue
552	continue
551	continue

!----

1000	continue
	write(*,*)' =================================================='
        write(*,*)' Now call shells to test which shells give zero for'
        write(*,*)' these special k points'
        write(*,*)'  '
        call shells(a1,a2,a3,sk,numb,weight,i1max,i2max,i3max)
	write(*,*)' =================================================='

	write(*,*)' We now have the k-points in the 1st BZ'
	write(*,*)' Write them to greatK.kpts'
        write(*,*)'  '
        open (unit=31,file='greatK.kpts',status='unknown')
	write(31,*)numb
	write(*,*)' We found numb=',numb,' kpoints.'
        write(*,*)' They are ....'
        do 441 ik=1,numb
        write(*,20)sk(1,ik),sk(2,ik),sk(3,ik),weight(ik)
        write(31,620)sk(1,ik),sk(2,ik),sk(3,ik),weight(ik)
620     format(1x,f13.8,1x,f13.8,1x,f13.8,8x,f13.8)
441     continue
!        close(unit=44)
        close(unit=31)
	write(*,*)'  '
	write(*,*)' We found numb=',numb,' kpoints.'
        write(*,*)'  '
        write(*,*)' We write them out to greatK.kpts'
        write(*,*)'  '
        end 

!****************************************************************************

      function ipermute3(i,j,k)
      implicit none
      integer ipermute3,mat(3,3),i,j,k,ndet
      mat=0
      mat(1,i)=1
      mat(2,j)=1
      mat(3,k)=1
      ipermute3=ndet(mat)
      end

!**************************************************************************

! find rotation matrices on basis vectors of reciprocal lattice,
! given the rotation matrices on basis vectors of direct lattice
! arguments:
!     iops_count (input), number of rotation matrices
!     iops(i,j,k) (input), kth rotation matrix of direct lattice
!     iops2(i,j,k) (output), kth rotation matrix of reciprocal lattice

      subroutine invsym(iops_count,iops,iops2)
      implicit none
      integer iops_count,iops(3,3,48),iops2(3,3,48),  &
     &     i,j,k,m,n,ip,jp,kp,ndet,ipermute3
      iops2=0
      do n=1,iops_count
        do i=1,3
        do j=1,3
        do k=1,3
          if(ipermute3(i,j,k).ne.1)cycle
          do ip=1,3
          do jp=1,3
          do kp=1,3
            m=ipermute3(ip,jp,kp)
            if(m.eq.0)cycle
            iops2(ip,i,n)=iops2(ip,i,n)+m*iops(jp,j,n)*iops(kp,k,n)
          enddo
          enddo
          enddo
        enddo
        enddo
        enddo
        m=ndet(iops2(1,1,n))*ndet(iops(1,1,n))
        if(m.eq.-1)then
          do i=1,3
          do j=1,3
            iops2(j,i,n)=-iops2(j,i,n)
          enddo
          enddo
        else if(m.ne.1)then
          call bomb
        endif
      enddo
      end

!***************************************************************************

      subroutine find_symmetry3(accuracy,basis,natoms,iatom_type, &
     &     atomposin,nbasis2,basis2,nops,opmatrix,opvector)

! find the space-group operators for a crystal
! each space group operator has the form {R|v}

! arguments:
!     accuracy (input), accuracy of the data
!     basis(i,j) (nput), ith cartesian coordinate of the jth basis vector 
!       of the lattice
!     natoms (input), number of atoms
!     iatom_type(i) (input), type of ith atom
!     atomposin(i,j) (input), ith cartesian coordinate of jth atom
!     nbasis2 (output), number of equivalent positions in the unit cell
!     basis2(i,j) (output), ith cartesian coordinate of the jth equivalent
!       position
!     nops (output), number of representative operators in space group
!     opmatrix(i,j,k) (output), matrix R of the kth rep
!     opvector(i,j) (output), ith cartesian coordinate of the vector v of 
!       the jth rep

! explicitly declare every variable
      implicit none
! declarations
      integer natoms
      double precision volunit,accuracy2,dtrace,   &
     &     prim_to_cart(3,3),v3(3),d,   &
     &     cart_to_prim(3,3),atompos(:,:),fract(3),v(3),v2(3),   &
     &     point(:,:),buff(3,3),vol,ops2(3,3,48),ops(3,3,48),   &
     &     unit_to_cart(3,3),accuracy,ddet3,buff2(3,3),   &
     &     cart_to_unit(3,3),vlattice(:,:),prim_to_unit(3,3),   &
     &     unit_to_prim(3,3),v1(3),   &
     &     vlength,primlength(3),unitlength(3),   &
     &     basis(3,3),basis2(3,natoms),opmatrix(3,3,48),   &
     &     opvector(3,48),atomposin(3,natoms)
      integer i,j,k,m,n,i1,i2,i3,iatom_type(natoms),ier,iatom,   &
     &     iatom2,itype2,nperiod,isgopcount,iop,itype1,ios,   &
     &     iatom1,iatom1p,iatom2p,map_atom(:),   &
     &     natoms_save,iop_order(48),   &
     &     iops_count,   &
     &     iv(4),iv2(4),lattice_count,map_lattice(:,:),iop2,   &
     &     iatom_count(:),nbasis2,nops
      logical iszero,done(:),foundone
! declare variables for which memory will be dynamically allocated
      allocatable atompos,point,vlattice, map_atom,   &
     &     map_lattice,done,iatom_count
!-----------------------------------------------------------------------------
! transformation from dimensionless to cartesian coordinates
      unit_to_cart(1:3,1:3)=basis(1:3,1:3)
! transformation from cartesian to dimensionless coordinates
      call xmatinv(unit_to_cart,cart_to_unit,ier)
!      print*,'cart_to_unit',cart_to_unit
! error flag is set if unit_to_cart is a singular matrix
      if(ier.ne.0)then
        write(6,*)'Error 004 in find_symmetry3'
        write(6,*)'Invalid input:  basis vectors of the lattice '//  &
     &       'are coplanar'
        call bomb
      endif
! get volume of unit cell
      volunit=dabs(ddet3(unit_to_cart))
!------------------------------------------------------------------------------
! number of atoms in unit cell
      natoms_save=natoms
! allocate memory for arrays that store data about atoms
      allocate(atompos(3,natoms),point(3,natoms+3),vlattice(3,natoms), &
     &     map_atom(natoms),map_lattice(natoms,natoms),done(natoms))
! put atoms inside unit cell
      do i=1,natoms
        call unitcell_fb(cart_to_unit,unit_to_cart,atomposin(1,i),  &
     &       atompos(1,i))
!        print*,'atomposin',atomposin
!        print*,'atompos',atompos
        do j=1,i-1
          call dvsub(atompos(1,i),atompos(1,j),v,3)
          if(iszero(vlength(3,v),accuracy))then
            write(6,*)'Error 011 in find_symmetry: Atoms ',j,' and ',i,  &
     &             ' are at the same position'
!            print *,atompos
            call bomb
          endif
        enddo
      enddo
!----------------------------------------------------------------------------
! initialize some variables
      foundone=.false.
      lattice_count=1
      do i=1,3
        vlattice(i,1)=0
      enddo
      do iatom=1,natoms
        map_lattice(iatom,1)=iatom
      enddo
!----------------------------------------------------------------------------
! find the type of atom which has lowest abundance
! mark each atom
      do iatom=1,natoms
        done(iatom)=.false.
      enddo
      k=natoms+1
! try each unmarked atom
      do iatom=1,natoms
      if(.not.done(iatom))then
! find, mark, and count all atoms of that type
        m=iatom_type(iatom)
        n=0
        do i=iatom,natoms
          if(iatom_type(i).eq.m)then
            done(i)=.true.
            n=n+1
          endif
        enddo
! find type of lowest abundance
        if(n.lt.k)then
          k=n
! save type of atom
          itype1=m
! save first atom of that type
          iatom1=iatom
        endif
      endif
      enddo
!------------------------------------------------------------------------------
! if the unit cell is not  primitive but contains additional
! lattice points, then there must exist lattice vectors which are not integer
! combinations of the basis vectors of the unit cell.  We will look for these
! new lattice vectors.  Each lattice vector must take us from each atom in
! the unit cell to another atom of the same type.  So, if we examine every
! vector that takes us from iatom1 to another atom of type itype1, we will have
! tried all possibilities for new lattice vectors.  We do this search using
! the atoms of the type of lowest abundance in order to make the search as
! short at possible.
! try all atoms of the same type as iatom1
!*****skip finding new lattice points************
      iatom1ploop: do iatom1p=iatom1+1,0
!      iatom1ploop: do iatom1p=iatom1+1,natoms
        if(iatom_type(iatom1p).eq.itype1)then
! the vector from iatom1 to iatom1p is a possible lattice vector: fract
          call dvsub(atompos(1,iatom1p),atompos(1,iatom1),fract,3)
! put this vector as close to the origin as possible
          call near_origin(cart_to_unit,unit_to_cart,fract,fract)
!------------------------------------------------------------------------------
! test this possible lattice vector on every atom.  We are going to add the
! tentative lattice vector to the position of each atom and then look for
! another atom of the same type which is nearest the result.
          do iatom2=1,natoms
            itype2=iatom_type(iatom2)
            map_atom(iatom2)=0
! add the lattice vector to the coordinates of the atom.  the result is
! contained in v.
            call dvadd(fract,atompos(1,iatom2),v,3)
! try every atom of the same type as iatom2
            do iatom2p=1,natoms
              if(iatom_type(iatom2p).eq.itype2)then
! vector from position of iatom2p to v
                call dvsub(v,atompos(1,iatom2p),v2,3)
! bring this vector near the origin
                call near_origin(cart_to_unit,unit_to_cart,v2,v3)
! find length of v3
                d=0
                do i=1,3
                  d=d+v3(i)**2
                enddo
                d=sqrt(d)
! is d small enough?
                if(iszero(d,accuracy))then
                  map_atom(iatom2)=iatom2p
                  goto 1
                endif
! next iatom2p
              endif
            enddo
! next iatom1p
            cycle iatom1ploop
! next iatom2
1         continue
          enddo
!-----------------------------------------------------------------------------
! we now have a mapping from every atom to another atom of the same type.
! test mappings for consistancy
! mark each atom
          do iatom=1,natoms
            done(iatom)=.false.
          enddo
! find cycles in the mapping
! try each unmarked atom
          do iatom=1,natoms
          if(.not.done(iatom))then
! mark it
            done(iatom)=.true.
! find period of cycle.   start with atom iatom
            m=1
            i=iatom
            do while(map_atom(i).ne.iatom)
              m=m+1
! map atom i onto next atom in cycle
              i=map_atom(i)
! if an atom in this cycle was already in another cycle,
! data is not consistant.  Try another tentative lattice vector
              if(done(i))cycle iatom1ploop
! mark this atom
              done(i)=.true.
            enddo
! first cycle tried
            if(iatom.eq.1)then
! cycles must have a period greater than 1.  If not, try another tentative
! lattice vector
              if(m.eq.1)cycle iatom1ploop
! save period of cycle
              nperiod=m
            else
! all cycles must have same period.  If not, try another tentative lattice
! vector
              if(m.ne.nperiod)cycle iatom1ploop
            endif
          endif
          enddo
!-----------------------------------------------------------------------------
! add to list of primitive lattice vectors which are fractionals with respect
! to the original primitive lattice
          lattice_count=lattice_count+1
! bring lattice vector into unit cell at origin
          call unitcell_fb(cart_to_unit,unit_to_cart,fract,v)
! save it
          do i=1,3
            vlattice(i,lattice_count)=v(i)
          enddo
! save mapping of lattice vector
          do iatom=1,natoms
            map_lattice(iatom,lattice_count)=map_atom(iatom)
          enddo
! try another tentative lattice vector (next iatom1p)
        endif
      enddo iatom1ploop
!-----------------------------------------------------------------------------
! adjust each new lattice vector to be exact
      do i=2,lattice_count
! multiply vector by number of new lattice vectors.  This should result
! in a lattice vector of the original primitive lattice
        do j=1,3
          fract(j)=vlattice(j,i)*lattice_count
        enddo
! vector in terms of basis vectors of the primitive lattice
        call xvmlt(cart_to_unit,fract,fract,3,3,3)
        do j=1,3
! each component should be an integer:  adjust them to be exact
          m=idnint(fract(j))
          if(.not.iszero(fract(j)-m,accuracy*lattice_count  &
     &         /unitlength(j)))then
            write(6,*)'Error 006 in find_symmetry'
            call bomb
          endif
          fract(j)=m
        enddo
! back to cartesian coordinate and the original vector
        call xvmlt(unit_to_cart,fract,fract,3,3,3)
        do j=1,3
          vlattice(j,i)=fract(j)/lattice_count
        enddo
      enddo
! copy results to output
      nbasis2=lattice_count
      basis2(1:3,1:lattice_count)=vlattice(1:3,1:lattice_count)
!------------------------------------------------------------------------------
! find the basis vectors of the lattice. Collect together the basis vectors of
! the unit cell plus the additional lattice vectors.  We must find new basis
! vectors such that each of these vectors in the collection can be expressed as
! an integer combination of them.
! get basis vectors of the unit cell
      do i=1,3
        do j=1,3
          point(j,i)=unit_to_cart(j,i)
        enddo
      enddo
! add to the list each new lattice point
      n=3
      do i=2,lattice_count
        n=n+1
        do j=1,3
          point(j,n)=vlattice(j,i)
        enddo
      enddo
! volume of new primitive unit cell
      vol=volunit/lattice_count
! try all possible triplets of vectors from the collection as the new
! primitive lattice vectors
      do i1=1,n-2
      do i2=i1+1,n-1
      do i3=i2+1,n
! tentative basis vectors of the primitive lattice
        do i=1,3
          prim_to_cart(i,1)=point(i,i1)
          prim_to_cart(i,2)=point(i,i2)
          prim_to_cart(i,3)=point(i,i3)
        enddo
! check volume of primitive unit cell
        d=dabs(ddet3(prim_to_cart))
! if the volume is not correct, try the next triplet of vectors
        if(.not.iszero(vol-d,accuracy*vol*(1/primlength(1)    &
     &       +1/primlength(2)+1/primlength(3))))goto 250
! get inverse of prim_to_cart
        call xmatinv(prim_to_cart,cart_to_prim,ier)
! if prim_to_cart is non-singular
        if(ier.eq.0)then
! check if every vector in the collection can be expressed as an integer
! combination of these tentative primitive lattice vectors.
! try each vector in collection
          do i=1,n
! express vector as linear combination of the tentative primitive lattice
! vectors
            call xvmlt(cart_to_prim,point(1,i),v,3,3,3)
! are all components integers?
            do j=1,3
              m=nint(v(j))
! if not an integer, try another triplet of vectors
              if(.not.iszero(v(j)-m,accuracy                  &
     &             /vlength(3,cart_to_prim(1,j))))goto 250
            enddo
          enddo
! every vector successful: we have found the new primitive lattice vectors
          goto 251
        endif
! next triplet of vectors (next i1,i2,i3)
250     continue
      enddo
      enddo
      enddo
! no triplet of vectors is successful.  I don't anticipate that this error is
! possible to occur.  I would sure like to hear about it if it happens.
      write(6,*)'Error 012 in find_symmetry'
      call bomb
! found basis vectors of the primitive lattice.  They are contained in
! prim_to_cart which is also the transformation from coordinates with respect
! to basis vectors of the primitive lattice to cartesian coordinates
251   continue
! improve choice: shortest possible basis vectors
      call nice_lattice3(prim_to_cart)
      accuracy2=0
      do i=1,3
        primlength(i)=vlength(3,prim_to_cart(1,i))
        if(accuracy/primlength(i).gt.accuracy2)  &
     &       accuracy2=accuracy/primlength(i)
      enddo
! get cart_to_prim:  transformation from cartesian coordinates to coorindates
! with respect to basis vectors of the primitive lattice
      call xmatinv(prim_to_cart,cart_to_prim,ier)
! get prim_to_unit:  transformation from coordinates with respect to basis
! vectors of the primitive lattice to basis vectors of the unit cell
      call xmatmlt(cart_to_unit,prim_to_cart,prim_to_unit,3,3,3,3,3,3)
! get unit_to_prim:  transformation from coordinates with respect to basis
! vectors of the unit cell to basis vectors of the primitive lattice
      call xmatmlt(cart_to_prim,unit_to_cart,unit_to_prim,3,3,3,3,3,3)
!------------------------------------------------------------------------------
! get point operators of the lattice
      call lattice_symmetry2(prim_to_cart,iops_count,ops,accuracy)
! find order of each point operator, i.e., find value of n for which R^n=E.
      do i=1,iops_count
        n=1
        buff(1:3,1:3)=ops(1:3,1:3,i)
        do while(.not.iszero(dtrace(buff,3,3)-3,10*accuracy))
          call xmatmlt(ops(1,1,i),buff,buff,3,3,3,3,3,3)
          n=n+1
          if(n.gt.6)then
            write(6,*)'Error 001 in find_symmetry3'
          endif
        enddo
        iop_order(i)=n
      enddo
!------------------------------------------------------------------------------
! map atoms into primitive unit cell
! remove marks from atoms
      do i=1,natoms_save
        done(i)=.false.
      enddo
      natoms=0
! do each unmarked atom
      do iatom=1,natoms_save
      if(.not.done(iatom))then
! count atoms in primitive unit cell
        natoms=natoms+1
! identify iatom1 in new list of atoms
        if(iatom.eq.iatom1)iatom1=natoms
! get position of atom
        do j=1,3
          atompos(j,natoms)=atompos(j,iatom)
        enddo
! get type of atom
        iatom_type(natoms)=iatom_type(iatom)
! do each lattice vector
        do i=2,lattice_count
! find atom which this lattice vector maps iatom onto
          j=map_lattice(iatom,i)
! atom should not already be marked
          if(done(j))then
            write(6,*)'Error 007 in find_symmetry'
            call bomb
          endif
! mark atom
          done(j)=.true.
        enddo
! bring atom into primitive unit cell at origin
        call unitcell_fb(cart_to_prim,prim_to_cart,atompos(1,natoms),  &
     &        atompos(1,natoms))
      endif
      enddo
!-----------------------------------------------------------------------------
! find elements of space group
! count them
      isgopcount=0
      foundone=.false.
! try each point operator
      ioploop: do iop=1,iops_count
!------------------------------------------------------------------------------
! The basis vectors of the lattice must reflect the symmetry
! of the point operator.  The point operator must bring each basis vector of
! of the primitive lattice into another vector of the primitive lattice.
! rotate the primitive lattice vectors
      call xmatmlt(ops(1,1,iop),prim_to_cart,buff,3,3,3,3,3,3)
! put result in terms of the primitive lattice vectors
      call xmatmlt(cart_to_prim,buff,ops2(1,1,iop),3,3,3,3,3,3)
! if the primitive lattice vectors have the exact symmetry required, the result
! would be all integers.
      do j=1,3
      do k=1,3
        m=idnint(ops2(k,j,iop))
        if(.not.iszero(ops2(k,j,iop)-m,accuracy/primlength(k)))then
          write(6,*)'Error 013 in find_symmetry'
          call bomb
        endif
        ops2(k,j,iop)=m
      enddo
      enddo
!------------------------------------------------------------------------------
! try this point operator on the atoms and look for a space-group element.
! Operate on atom iatom1.  Remember that iatom1 is the first atom of the type
! of lowest abundance.  v1 contains the coordinates of the atom after
! the operation.  If any elements of the space group contain this point
! operator, then they will bring this atom to another atom of the same type
      call xvmlt(ops(1,1,iop),atompos(1,iatom1),v1,3,3,3)
! try to map the rotated atom onto every other atom of the same type
      do iatom1p=1,natoms
        if(iatom_type(iatom1p).eq.itype1)then
! the space group element consists of point operation following by a translation
! (called the fractional).  find the fractional required in the space group
! element if it is to bring iatom1 to iatom1p
          call dvsub(atompos(1,iatom1p),v1,fract,3)
! make this fractional as short as possible
          call near_origin(cart_to_prim,prim_to_cart,fract,  &
     &          fract)
!------------------------------------------------------------------------------
! try this tentative space-group operator on every atom
          do iatom2=1,natoms
            itype2=iatom_type(iatom2)
            map_atom(iatom2)=0
! operate on the atom: point operation followed by a translation.
! v contains the coordinates of the atom after the operation by the space
! group element.
            call xvmlt(ops(1,1,iop),atompos(1,iatom2),v,3,3,3)
            call dvadd(fract,v,v,3)
! find atom of same type closest to v.
! try every atom of the same type as iatom2
            do iatom2p=1,natoms
              if(iatom_type(iatom2p).eq.itype2)then
! vector from position of iatom2p to v
                call dvsub(v,atompos(1,iatom2p),v2,3)
! bring this vector near the origin
                call near_origin(cart_to_prim,prim_to_cart,  &
     &                v2,v3)
! find length of v3
                d=0
                do i=1,3
                  d=d+v3(i)**2
                enddo
                d=sqrt(d)
! is d short enough?
                if(iszero(d,accuracy))then
                  map_atom(iatom2)=iatom2p
                  goto 2
                endif
! next iatom2p
              endif
            enddo
! next iatom1p
            goto 301
! next iatom2
2         continue
          enddo
!------------------------------------------------------------------------------
! test mapping for consistancy
! find cycles in the mapping
! mark each atom
          do iatom=1,natoms
            done(iatom)=.false.
          enddo
! try each unmarked atom
          do iatom=1,natoms
          if(.not.done(iatom))then
! mark it
            done(iatom)=.true.
! find atoms in cycle.  start with atom iatom
            m=1
            i=iatom
            n=0
            do while(map_atom(i).ne.iatom)
              m=m+1
! map atom i onto next atom in cycle
              i=map_atom(i)
! if an atom in this cycle was already in another cycle,
! data is not consistant.  Try another fractional
              if(done(i))goto 301
! mark this atom
              done(i)=.true.
            enddo
! period of cycle must be consistent with period of point operator
            if(mod(iop_order(iop),m).ne.0)goto 301
          endif
          enddo
! found a space group element
          goto 302
!------------------------------------------------------------------------------
! try another fractional (next iatom1p)
301       continue
        endif
      enddo
! no fractional can be found.  This point operator in not part of any space
! group element.  Try next point operator.
      cycle ioploop
!------------------------------------------------------------------------------
! we have found a space group element
! count them
302   isgopcount=isgopcount+1
! save point operator part of element
      opmatrix(1:3,1:3,isgopcount)=ops(1:3,1:3,iop)
! save translational part of element
      opvector(1:3,isgopcount)=fract(1:3)
! next point operator (iop)
      enddo ioploop
! output number of operators
      nops=isgopcount
!------------------------------------------------------------------------------
! deallocate memory for arrays that store data about atoms
      deallocate(atompos,point,vlattice,map_atom,map_lattice,done)
      end

!**************************************************************************

      subroutine lattice_symmetry2(basis,iops_count,xmat,eps)

! find the symmetry elements of a lattice, given three basis vectors
! arguments:
!     basis(i,j) (input), ith cartesian component of the jth basis vector
!     iops_count (output), number of elements in the point group
!     xmat(i,j,k) (output), rotation matrix for the kth element
!     eps (input), tolerance on lengths of vectors

      implicit none
      double precision basis(3,3),xmat(3,3,48),temp(3,3),   &
     &     basis_inv(3,3),eps
      integer nmatrices,matrices(3,3,48),ier,i,j,k,m,n,itemp(3,3),  &
     &     iops_count

! symmetry matrices
      call dlatmat2(basis,eps,nmatrices,matrices)
      call xmatinv(basis,basis_inv,ier)
      if(ier.ne.0)then
        write(6,*)'Error 001 in lattice_symmetry'
        call bomb
      endif
! multiplication table
 1    continue
      do i=1,nmatrices
        do j=1,nmatrices
          call matmlt(matrices(1,1,i),matrices(1,1,j),itemp)
          kloop: do k=1,nmatrices
            mloop: do m=1,3
              do n=1,3
                if(matrices(n,m,k).ne.itemp(n,m))exit mloop
              enddo
              if(m.eq.3)then
                exit kloop
              endif
            enddo mloop
! missing a matrix: add it to list and try making multiplication table again
            if(k.eq.nmatrices)then
              nmatrices=nmatrices+1
              if(nmatrices.gt.48)then
                write(6,*)'Error 006 in lattice_symmetry'
                call bomb
              endif
              matrices(1:3,1:3,nmatrices)=itemp(1:3,1:3)
              goto 1
            endif
          enddo kloop
        enddo
      enddo
! rotation matrices
      iops_count=nmatrices
      do i=1,nmatrices
        do j=1,3
        do k=1,3
          temp(k,j)=matrices(k,j,i)
        enddo
        enddo
        call xmatmlt(temp,basis_inv,temp,3,3,3,3,3,3)
        call xmatmlt(basis,temp,xmat(1,1,i),3,3,3,3,3,3)
      enddo
      
      end

!***************************************************************************

      subroutine nice_lattice3(basis)

! find nice lattice vectors: based on length
! argument:
!     basis(i,j) (input/output), ith cartesian component of jth lattice vector

      implicit none
      double precision basis(3,3),v(3),d(3),x,vlength,ddet
      integer i,j,ncmpd,npos,nneg
      logical foundone
! length of each basis vector
      do i=1,3
        d(i)=vlength(3,basis(1,i))
      enddo
      foundone=.true.
! repeat until no further improvement can be made
      do while(foundone)
      foundone=.false.
! try each pair of basis vectors
      do i=2,3
        do j=1,i-1
! add them together: is the result shorter than one of the two vectors?
! if so, make a substitution
          call dvadd(basis(1,i),basis(1,j),v,3)
          x=vlength(3,v)
          if(x.lt.d(i).and.ncmpd(x-d(i)).ne.0)then
            basis(1:3,i)=v(1:3)
            d(i)=x
            foundone=.true.
            cycle
          else if(x.lt.d(j).and.ncmpd(x-d(j)).ne.0)then
            basis(1:3,j)=v(1:3)
            d(j)=x
            foundone=.true.
            cycle
          endif
! if adding did not produce a shorter vector, try subtracting them
          call dvsub(basis(1,i),basis(1,j),v,3)
          x=vlength(3,v)
          if(x.lt.d(i).and.ncmpd(x-d(i)).ne.0)then
            basis(1:3,i)=v(1:3)
            d(i)=x
            foundone=.true.
            cycle
          else if(x.lt.d(j).and.ncmpd(x-d(j)).ne.0)then
            basis(1:3,j)=v(1:3)
            d(j)=x
            foundone=.true.
            cycle
          endif
        enddo
      enddo
      enddo
! done: get rid of as many minus signs as possible
      do i=1,3
        npos=0
        nneg=0
        do j=1,3
          if(ncmpd(basis(j,i)).ne.0)then
            if(basis(j,i).gt.0.0d0)then
              npos=npos+1
            else
              nneg=nneg+1
            endif
          endif
        enddo
        if(nneg.gt.npos)then
          do j=1,3
            basis(j,i)=-basis(j,i)
          enddo
        endif
      enddo
! make the determinant positive if needed
      x=ddet(basis,3,3)
      if(x.lt.0.0d0)then
        do j=1,3
          basis(j,3)=-basis(j,3)
        enddo
      endif
      end

!***************************************************************************

      subroutine xmatinv(xmatin,xmatout,ier)
      implicit none

! invert a 3 by 3 matrix

      double precision xmatin(3,3),xmatout(3,3),buffer(3,3),x
      integer indx(3),ier,n,i,j

! dimension of matrix
      n=3
! clear error flag
      ier=0
      do i=1,n
        do j=1,n
          xmatout(i,j)=0
          buffer(i,j)=xmatin(i,j)
        enddo
        xmatout(i,i)=1
      enddo
! decomposition
      call ludcmp(buffer,n,n,indx,x)
! singular matrix
      if(x.eq.0.0d0)then
        ier=1
        return
      endif
! inverse matrix
      do j=1,n
        call lubksb(buffer,n,n,indx,xmatout(1,j))
      enddo
      end

! The following routines are from Numerical Recipes

      subroutine ludcmp(a,n,np,indx,d)
      implicit none
      integer nmax,np,n
      double precision tiny
      parameter (nmax=3,tiny=1.0d-20)
      double precision a(np,np),vv(nmax),d,aamax,dum,sum
      integer indx(n),i,j,k,imax,ncmpd
      d=1
      do i=1,n
        aamax=0
        do j=1,n
          if(dabs(a(i,j)).gt.aamax)aamax=dabs(a(i,j))
        enddo
        if(ncmpd(aamax).eq.0)then
! singular matrix
          d=0
          return
        endif
        vv(i)=1/aamax
      enddo
      do j=1,n
        do i=1,j-1
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
        enddo
        aamax=0
        do i=j,n
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if(dum.ge.aamax)then
            imax=i
            aamax=dum
          endif
        enddo
        if(j.ne.imax)then
          do k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.0d0)a(j,j)=tiny
        if(j.ne.n)then
          dum=1/a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
      enddo
      end

      subroutine lubksb(a,n,np,indx,b)
      implicit none
      integer n,np
      double precision a(np,np),b(n),sum
      integer indx(n),ii,i,j,ll
      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if(ii.ne.0)then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo
        else if(sum.ne.0.0d0)then
          ii=i
        endif
        b(i)=sum
      enddo
      do i=n,1,-1
        sum=b(i)
        if(i.lt.n)then
          do j=i+1,n
            sum=sum-a(i,j)*b(j)
          enddo
        endif
        b(i)=sum/a(i,i)
      enddo
      end

!****************************************************************************

      subroutine xmatmlt(x1,x2,x3,nrow1,ncol1,ncol2,nr1,nr2,nr3)

! multiply two real matrices, x3=x1*x2
! double precision version
! arguments:
!     x1,x2 (input), first and second matrix
!     x3 (output), product x1*x2
!     nrow1 (input), number of rows in x1, also the number of rows in x3
!     ncol1 (input), number of columns in x1, also the number of
!          rows in x2
!     ncol2 (input), number of columns in x2, also the number of
!          columns in x3
!     nr1 (input), number of rows in the physical array x1
!     nr2 (input), number of rows in the physical array x2
!     nr3 (input), number of rows in the physical array x3
      
      implicit none
      integer i,j,k,nrow1,ncol1,ncol2,nr1,nr2,nr3
      double precision x1(nr1,ncol1),x2(nr2,ncol2),x3(nr3,ncol2),x(:,:)
      allocatable x
      allocate(x(nrow1,ncol2))
      do i=1,ncol2
      do j=1,nrow1
        x(j,i)=0
        do k=1,ncol1
          x(j,i)=x(j,i)+x1(j,k)*x2(k,i)
        enddo
      enddo
      enddo
      do i=1,ncol2
      do j=1,nrow1
        x3(j,i)=x(j,i)
      enddo
      enddo
      deallocate(x)
      
      end

!****************************************************************************

! determine if a real variable is zero
      function iszero(x,accuracy)
      implicit none
      double precision x,accuracy
      logical iszero
      if(dabs(x).lt.accuracy)then
        iszero=.true.
      else
        iszero=.false.
      endif
      end

!****************************************************************************

      subroutine matmlt(mat1,mat2,mat3)
      implicit none
!
!	MULTIPLY TWO 3X3 MATRICES
!	MAT3=MAT1*MAT2
!
      integer mat1(3,3),mat2(3,3),mat3(3,3),mat4(3,3),j,k,l
      do 1 j=1,3
      do 1 k=1,3
      mat4(j,k)=0
      do 1 l=1,3
1     mat4(j,k)=mat4(j,k)+mat1(j,l)*mat2(l,k)
      do 2 j=1,3
      do 2 k=1,3
2     mat3(j,k)=mat4(j,k)
      return
      end

!*****************************************************************************

      subroutine reduc2(ivec)
      implicit none
! remove any common factors from ivec
      integer ivec(4),i,factor
! common denominator cannot be zero
      if(ivec(4).eq.0)call bomb
! transfer minus sign from common denominator to components
      if(ivec(4).lt.0)then
        do i=1,4
          ivec(i)=-ivec(i)
        enddo
      endif
! remove common factors
      i=factor(4,ivec) ! FIXME : this does not make sense (?) should be i=factor*ivec(4) ?
      end

!****************************************************************************

      subroutine reduc1(n)
      implicit none
!
!	ADD AND SUBTRACT INTEGERS TO BRING EACH COMPONENT OF THE VECTOR N
!	TO A VALUE GREATER THAN OR EQUAL TO 0 BUT LESS THAN 1.
!
!	THIS IS MATHEMATICALLY EQUIVALENT TO THE FOLLOWING:
!	REDUCE  N(1)+J1*N(4),  N(2)+J2*N(4),  N(3)+J3*N(4)
!	TO SMALLEST POSSIBLE NON-NEGATIVE INTEGERS FOR ALL VALUES OF J1,J2,J3
!
      integer n(4),j
      if(n(4).eq.0)call bomb
      do 3 j=1,3
2     if(n(j)-n(4).lt.0)goto 1
      n(j)=n(j)-n(4)
      goto 2
1     if(n(j)+n(4).ge.n(4))goto 3
      n(j)=n(j)+n(4)
      goto 1
3     continue
      return
      end

!*****************************************************************************

      subroutine xvmlt(x,v1,v2,nrow,ncol,nr)

! multiply a double precision vector by a double precision matrix, v2=x*v1
! arguments:
!     x (input), matrix
!     v1 (input), vector
!     v2 (output), product x*v1
!     nrow (input), number of rows in x, also the number of rows in v2
!     ncol (input), number of columns in x, also the number of rows in v1
!     nr (input), number of rows in the physical array x
      
      implicit none
      integer nrow,ncol,nr,i,j
      double precision x(nr,ncol),v1(ncol),v2(nrow),v(:)
      allocatable v
      allocate(v(nrow))
      do i=1,nrow
        v(i)=0
        do j=1,ncol
          v(i)=v(i)+x(i,j)*v1(j)
        enddo
      enddo
      v2(1:nrow)=v(1:nrow)
      deallocate(v)
      
      end

!****************************************************************************

      function vlength(n,v)
      implicit none
      integer n
      double precision v(n),vlength,x
      integer i
      x=0
      do i=1,n
        x=x+v(i)**2
      enddo
      vlength=dsqrt(x)
      end

!***************************************************************************

      function ncmpd(x)
      implicit none
!
!	COMPARE X WITH ZERO
!	NCMP=0 IF X IS CLOSE ENOUGH TO ZERO
!	NCMP=1 OTHERWISE
!	X IS DOUBLE PRECISION
!
      integer ncmpd
      double precision x,delta
      data delta/1.d-6/
      ncmpd=0
      if(dabs(x).gt.delta)ncmpd=1
      return
      end

!***************************************************************************

      function ndet(mat)
      implicit none
!
!	FIND THE DETERMINANT OF A 3X3 MATRIX MAT
!
      integer ndet,mat(3,3)
      ndet=mat(1,1)*(mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2))  &
     & -mat(1,2)*(mat(2,1)*mat(3,3)-mat(2,3)*mat(3,1))  &
     & +mat(1,3)*(mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1))
      return
      end

!**************************************************************************

      subroutine bomb
      implicit none
      write(6,'(a)')'This program has bombed'
      write(6,'(a)')'exit'
      stop
      end

!*************************************************************************

! determinant of matrix
      function ddet3(x)
      implicit none
      double precision x(3,3),ddet3
      ddet3=x(1,1)*(x(2,2)*x(3,3)-x(2,3)*x(3,2))  &
     &      +x(1,2)*(x(2,3)*x(3,1)-x(2,1)*x(3,3))  &
     &      +x(1,3)*(x(2,1)*x(3,2)-x(2,2)*x(3,1))
      end

!*************************************************************************

! bring a point into the unit cell at the origin

      subroutine unitcell_fb(cart_to_prim,prim_to_cart,v1,v2)
      implicit none
      double precision cart_to_prim(3,3),prim_to_cart(3,3),v1(3),  &
     &     v2(3),buff(3)
      integer i,ncmpd
! change coordinates of point to linear combination of basis vectors of the
! primitive lattice
      call xmatmlt(cart_to_prim,v1,buff,3,3,1,3,3,3)
! in the unit cell at the origin, the coefficient must be greater than or
! equal to zero and less than one.
      do i=1,3
        do while(buff(i).gt.1.0d0.or.ncmpd(buff(i)-1).eq.0)
          buff(i)=buff(i)-1
        enddo
        do while(buff(i).lt.0.0d0.and.ncmpd(buff(i)).ne.0)
          buff(i)=buff(i)+1
        enddo
      enddo
! return to cartesian coordinates
      call xmatmlt(prim_to_cart,buff,v2,3,3,1,3,3,3)
      end

!***************************************************************************

      subroutine dvsub(v1,v2,v3,nrow)

! subtract two real vectors: v3=v1-v2
! double precision version
! arguments:
!     v1,v2 (input), vectors
!     v3 (output), vector v1-v2
!     nrow (input), number of rows in each vector
      
      implicit none
      integer i,nrow
      double precision v1(nrow),v2(nrow),v3(nrow)
      do i=1,nrow
        v3(i)=v1(i)-v2(i)
      enddo
      end

!****************************************************************************

! put vector near origin by adding lattice vectors to it
      subroutine near_origin(xmat1,xmat2,v1,v2)
! xmat1 transforms cartesian coordinates to coordinates based on lattice
! xmat2 transforms coordinates based on lattice to cartesian coordinates
! v1 is the input vector in cartesian coordinates
! v2 is the output vector in cartesian coordinates
      implicit none
      double precision xmat1(3,3),xmat2(3,3),v1(3),v2(3),v(3)
      integer i
      call xvmlt(xmat1,v1,v,3,3,3)
      do i=1,3
        do while(v(i).gt.0.5d0)
          v(i)=v(i)-1
        enddo
        do while(v(i).le.-0.5d0)
          v(i)=v(i)+1
        enddo
      enddo
      call xvmlt(xmat2,v,v2,3,3,3)
      end

!****************************************************************************

      subroutine dvadd(v1,v2,v3,nrow)
! add two real vectors: v3=v1+v2
! double precision version
! arguments:
!     v1,v2 (input), vectors
!     v3 (output), vector v1+v2
!     nrow (input), number of rows in each vector     
      implicit none
      integer i,nrow
      double precision v1(nrow),v2(nrow),v3(nrow)
      do i=1,nrow
        v3(i)=v1(i)+v2(i)
      enddo
      end

!****************************************************************************

      function dtrace(xmat,n,nr)
! find the trace of an irrep matrix
! the trace of the matrix xmat of dimension n
      implicit none
      integer n,i,nr
      real*8 xmat(nr,n),dtrace,x
      x=0
      do i=1,n
        x=x+xmat(i,i)
      enddo
      dtrace=x
      return
      end

!*************************************************************************

      function factor(n,numbers)
      implicit none
! remove the greatest common factor contained in n integers in numbers
      integer n,numbers(1000),min,i,j,factor
      factor=1
! find a nonzero integer
      do i=1,n
        if(numbers(i).ne.0)goto 2
      enddo
! all zeros
      return
! find the minimum absolute nonzero value among the integers
2     min=iabs(numbers(i))
      do i=2,n
        if(numbers(i).ne.0.and.iabs(numbers(i)).lt.min)   &
     &      min=iabs(numbers(i))
      enddo
! try each number from the minimum on down to 2
      do i=min,2,-1
! is i a common factor?
        do j=1,n
          if(mod(numbers(j),i).ne.0)goto 1
        enddo
! yes, divide it out
        do j=1,n
          numbers(j)=numbers(j)/i
        enddo
! save it too
        factor=factor*i
! done
        return
! try next number
1       continue
      enddo
      end

!****************************************************************************

      subroutine dlatmat2(cart,eps,nmatrices,matrices)

! find symmetry matrices for a given lattice

! arguments:
!     cart(i,j) (input), ith cartesian component of jth basis vector
!     eps (input), tolerance for length
!     nmatrices (output), number of matrices
!     matrices(i,j,k) (output), kth matrix

      implicit none

      integer nmax
      parameter(nmax=100)

      integer n,i,j,j1,j2,j3,k,m,i1,i2,i3,nshort(3),ndet,itrans(3,3),  &
     &     nmatrices,matrices(3,3,48),ichoose(3)
      double precision eps, dshort(nmax,3),ishort(3,nmax,3),  &
     &     vshort(3,nmax,3),v(3),xmax, vlength,x,dvdot,d,abc(3,3),  &
     &     cart(3,3)
      logical foundone,tried(48,48)

! some initialization
      dshort=0
      ishort=0
      vshort=0
      nshort=0
      do i=1,3
        abc(i,i)=vlength(3,cart(1,i))
      enddo
      do i=2,3
        do j=1,i-1
          call dvsub(cart(1,i),cart(1,j),v,3)
          abc(i,j)=vlength(3,v)
          abc(j,i)=abc(i,j)
        enddo
      enddo
      i=0
      foundone=.true.
! longest lattice parameter
      xmax=dmax1(abc(1,1),abc(2,2),abc(3,3))+eps
! try each shell until every vector in a shell is longer than the longest
! lattice parameter
      do while(foundone)
        i=i+1
        foundone=.false.
! find all lattice vectors in shell
        do j1=-i,i
        do j2=-i,i
        do j3=-i,i
        if(iabs(j1).eq.i.or.iabs(j2).eq.i.or.iabs(j3).eq.i)then
! length of lattice vector
          v=0
          do k=1,3
            v(k)=v(k)+j1*cart(k,1)
            v(k)=v(k)+j2*cart(k,2)
            v(k)=v(k)+j3*cart(k,3)
          enddo
          d=vlength(3,v)
! if shorter than longest lattice parameter, then do next shell too
          if(d.lt.xmax)foundone=.true.
! check each lattice parameter a,b,c
          do k=1,3
! equal to length of lattice parameter to within tolerance
            if(dabs(d-abc(k,k)).lt.eps)then
! count them
              nshort(k)=nshort(k)+1
              if(nshort(k).gt.nmax)then
                write(6,*)'Error in dlatmat:  nmax too small'
                call bomb
              endif
! length
              dshort(nshort(k),k)=d
! dimensionless coordinates
              ishort(1,nshort(k),k)=j1
              ishort(2,nshort(k),k)=j2
              ishort(3,nshort(k),k)=j3
! cartesian coordinates
              vshort(1,nshort(k),k)=v(1)
              vshort(2,nshort(k),k)=v(2)
              vshort(3,nshort(k),k)=v(3)
            endif
          enddo
! next vector in shell
        endif
        enddo
        enddo
        enddo
! next shell
      enddo

! try mappings of basis vectors onto vectors the "same" length
      nmatrices=1
      do i1=1,nshort(1)
      ichoose(1)=i1
      itrans(1:3,1)=ishort(1:3,i1,1)
      do i2=1,nshort(2)
      ichoose(2)=i2
      itrans(1:3,2)=ishort(1:3,i2,2)
      i3loop: do i3=1,nshort(3)
      ichoose(3)=i3
      itrans(1:3,3)=ishort(1:3,i3,3)
! determinant of the transformation matrix must be equal to 1
      if(iabs(ndet(itrans)).ne.1)cycle
! lengths of differences of lattice vectors must match to within tolerance
      do i=2,3
        do j=1,i-1
          call dvsub(vshort(1,ichoose(i),i),vshort(1,ichoose(j),j),v,3)
          x=vlength(3,v)
          if(dabs(x-abc(i,j)).gt.eps)cycle i3loop
        enddo
      enddo
! found a transformation:  count them and save it
! if this is the identity op, just put it into the first matrix where we
! have reserved a place for it
      if(itrans(1,1)+itrans(2,2)+itrans(3,3).eq.3)then
        matrices(1:3,1:3,1)=itrans(1:3,1:3)
      else
        nmatrices=nmatrices+1
        if(nmatrices.gt.48)then
          write(6,*)'Error in dlatmat2: more than 48 point operators'
          call bomb
        endif
        matrices(1:3,1:3,nmatrices)=itrans(1:3,1:3)
      endif
! next mapping
      enddo i3loop
      enddo
      enddo
! find any additional matrices by multiplication
      foundone=.true.
      tried=.false.
      do while(foundone)
        foundone=.false.
        do i=1,nmatrices
          do j=1,nmatrices
            if(.not.tried(i,j))then
              tried(i,j)=.true.
              call matmlt(matrices(1,1,i),matrices(1,1,j),itrans)
              kloop: do k=1,nmatrices
                mloop: do m=1,3
                do n=1,3
                  if(matrices(n,m,k).ne.itrans(n,m))exit mloop
                  if(m.eq.3.and.n.eq.3)exit kloop
                enddo
                enddo mloop
                if(k.eq.nmatrices)then
                  foundone=.true.
                  nmatrices=nmatrices+1
                  if(nmatrices.gt.48)then
                    write(6,*)'Error in dlatmat2: '//  &
     &                   'more than 48 point operators'
                    call bomb
                  endif
                  matrices(1:3,1:3,nmatrices)=itrans(1:3,1:3)
                endif
              enddo kloop
            endif
          enddo
        enddo
      enddo
      end

!*************************************************************************

! determinant of a double precision matrix (adapted from xrowop2)
!    xmat (input), square matrix
!    nd (input), dimension of square matrix
!    nnd (input), actual number of rows in array xmat

      double precision function ddet(xmat,nd,nnd)

      implicit none
      integer nd,nnd
      double precision xmat(nnd,nd),xna(:,:),x
      integer nrow,ncol,nout,j,k,l,ncmpd
      allocatable xna

      allocate(xna(nd,nd))
      xna(1:nd,1:nd)=xmat(1:nd,1:nd)
      nrow=nd
      ncol=nd
      ddet=1
!       TRANSFORM MATRIX BY ROWS
      nout=0
      do 1 j=1,nrow
5     if(j+nout.gt.ncol)call bomb
!       PUT NON-ZERO ELEMENTS ON DIAGONAL
      if(ncmpd(xna(j,j+nout)).ne.0)goto 4
      do 2 k=j+1,nrow
2     if(ncmpd(xna(k,j+nout)).ne.0)goto 3
!       IF COLUMN HAS ALL ZEROES, determinant is zero
      ddet=0
      goto 99
!       INTERCHANGE ROWS
3     ddet=-ddet
      do 6 l=1,ncol
      x=xna(j,l)
      xna(j,l)=xna(k,l)
6     xna(k,l)=x
!       NORMALIZE ELEMENTS IN ROW
4     x=xna(j,j+nout)
      ddet=ddet*x
      do 20 k=j+nout,ncol
20    xna(j,k)=xna(j,k)/x
!       REDUCE ALL OTHER ROWS
      do 7 k=1,nrow
      if(k.eq.j)goto 7
      if(ncmpd(xna(k,j+nout)).eq.0)goto 7
      do 8 l=1,ncol
      if(l.eq.j+nout)goto 8
      xna(k,l)=xna(k,l)-xna(j,l)*xna(k,j+nout)
8     continue
      xna(k,j+nout)=0.
7     continue
1     continue
 99   deallocate(xna)
      end

!***************************************************************************
