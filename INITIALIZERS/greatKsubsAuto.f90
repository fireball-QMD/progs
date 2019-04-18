! ====================================================================
       subroutine gkcross(a,b,c)
       implicit none !double precision(a-h,o-z)
       real, intent(in) :: a,b
       real, intent(out) :: c
       dimension a(3), b(3), c(3)
! computes c = a X b.
       c(1)=a(2)*b(3)-a(3)*b(2)
       c(2)=a(3)*b(1)-a(1)*b(3)
       c(3)=a(1)*b(2)-a(2)*b(1)
       return
       end
!====================================================================
!kgroup.f90 -- Read in the point group data file, and gives the
!            k vectors that we have found which are irreducible.
!=======================================================================
!for questions or futher information:
!o.f. sankey, dept of physics, ASU, 602-965-4334, bitnet: SANKEY@ASUCPS.
!=======================================================================
subroutine kgroup(sk,weight,numb,skirred,wt,numirr,pcgrp,  &
     &     npgops,a1,a2,a3,rlv1,rlv2,rlv3,i1max,i2max,i3max)
    implicit none !double precision (a-h,o-z)
! passed variables -----------------------------------------------------        
    integer, intent(in) :: npgops, i1max, i2max, i3max
    integer, intent(out) :: numirr
    integer, intent(inout) :: numb
    real, intent(in) :: a1, a2, a3, rlv1,rlv2,rlv3
    real, intent(inout) :: sk, weight, skirred, wt, pcgrp
! internal variables -----------------------------------------------------
    character(len=50) filename
    integer i,j,nkmax,igrpsz,iwrite,inversion,ngroup,igrp,jgrp, &
     &   ix,iy,istop,mu,iadd,itest,irepeat,l,ic,iwhich,icount,ik, &
     &   iwedge,irred,jitis,i1,iaddmore,nirred,i2,i3
	real kxp,kyp,kzp,wait,xcomp,ycomp,zcomp,sum1
    parameter (nkmax=100000)
	parameter (igrpsz=96)
	dimension weight(nkmax)
    dimension a1(3),a2(3),a3(3)
	dimension rlv1(3),rlv2(3),rlv3(3)
    dimension sk(3,nkmax)
	dimension skirred(3,nkmax),irred(nkmax),nirred(nkmax)
	dimension wt(nkmax),wait(nkmax),irepeat(nkmax)
	dimension iwhich(nkmax)
	dimension pcgrp(3,3,igrpsz)
!	common /lattice/a1,a2,a3,rlv1,rlv2,rlv3
!        common /imaxes/i1max,i2max,i3max
!================================================================
!	write(*,*)'  '
!!	write(*,*)' Welcome to subroutine kgroup...'
!write(*,*)'  '
!	write(*,*)' We call kgroup, which rotates them to get '
!	write(*,*)' a fully symmetric set, then weeds out repeats,'
!	write(*,*)' and optionally keeps only those in the'
!	write(*,*)' irreducible wedge.'
!
!For extended output set iwrite=1, otherwise iwrite=0
!Terse
	iwrite=0
!Extended
!	iwrite=1
!================================================================
!$$$	write(*,*)' Insert the file name for the point group operators.'
!$$$	write(*,*)' Note: the first operator must be the identity.'
!$$$	read(*,539)filename
!$$$539	format(a50)
!$$$	write(*,539)filename
!$$$c	write(*,539)filename
!$$$	open(unit=30,file=filename,status='old')
!$$$	read(30,*)ngroup
!$$$	write(*,*)' There are',ngroup,' group operators.'
!$$$	if(ngroup.gt.igrpsz) stop 'too many group operators'
!$$$	do 540 igrp=1,ngroup
!$$$	read(30,*)jgrp
!$$$	if(jgrp.ne.igrp) stop ' something wrong in group file'
!$$$c pcgrp are the group operators in cartesian coordinates.
!$$$c (x',y',z') = pcgrp * (x,y,z) where x,y,z is x*xhat+y*yhat+z*zhat.
!$$$	do 541 i=1,3
!$$$	read(30,*) pcgrp(igrp,i,1),pcgrp(igrp,i,2),pcgrp(igrp,i,3)
!$$$c
!$$$	if(iwrite.eq.1)then
!$$$	write(*,188) pcgrp(igrp,i,1),pcgrp(igrp,i,2),pcgrp(igrp,i,3)
!$$$	end if
!$$$c
!$$$188	format('  ',3f12.8)
!$$$541	continue
!$$$540	continue
!$$$	close(unit=30)
!	write(*,*)'  '
!	write(*,*)' For special k-points, we always have'
!	write(*,*)' that k ---> -k has the same energy.'
!	write(*,*)' DO YOU WISH TO INCLUDE INVERSION?'
!	write(*,*)' insert 1 for YES, 0 for NO... usually insert 1'
!	read(*,*)inversion
        inversion=1  ! always add inversion for time reversal
!	write(*,*)' inversion=',inversion
!	inversion=1
!	write(*,*)'  '

! ========================================================
!	write(*,*)'  '
        ngroup=npgops
	if(inversion.eq.1)then
	write(*,*)' DOUBLE THE GROUP .....'
	do igrp=1,ngroup
	 jgrp=igrp+ngroup
	 do ix=1,3
 	  do iy=1,3
 	   pcgrp(ix,iy,jgrp)=-pcgrp(ix,iy,igrp)
	  end do
	 end do
    end do
	ngroup=ngroup+ngroup
	write(*,*)' There are',ngroup,' group operators.'
	 if(ngroup.gt.igrpsz) stop 'too many group operators'
	end if
! ============================================================
! We must have the first group operator being the identity.
! Check that first operator is identity.
	do i=1,3
 	 if(abs(pcgrp(i,i,1)-1.0d0).gt.0.0000001)then
	 write(*,*)' the first group operator is not the identity'
	 write(*,*)' We assume that group op. # 1=identity and it is not'
	 stop ' sorry!!!'
	 end if
	 do j=1,3
	  if(j.eq.i) cycle
	  if(abs(pcgrp(i,j,1)).gt.0.0000001)then
	   write(*,*)' the first group operator is not the identity'
	   write(*,*)' We assume that group op. # 1=identity and it is not'
	   stop ' sorry!!!'
	  end if
     end do
    end do
! =================================================================
! We make ngroup*numb points by rotating the points we have by ngroup operators.
	istop=0
	write(*,*)' Now rotate the k vectors we have using the'
	write(*,*)' group operators.'
	write(*,*)' numb=',numb
	write(*,*)' ngroup=',ngroup
    do i=1,numb
	 weight(i)=weight(i)/dfloat(ngroup)
     do mu=1,ngroup
      kxp=pcgrp(1,1,mu)*sk(1,i)+pcgrp(1,2,mu)*sk(2,i)+  &
      &  pcgrp(1,3,mu)*sk(3,i)
      kyp=pcgrp(2,1,mu)*sk(1,i)+pcgrp(2,2,mu)*sk(2,i)+  &
      &  pcgrp(2,3,mu)*sk(3,i)
      kzp=pcgrp(3,1,mu)*sk(1,i)+pcgrp(3,2,mu)*sk(2,i)+  &
      &  pcgrp(3,3,mu)*sk(3,i)
	  iadd=i+numb*(mu-1)
	  write(*,*)' i,mu=',i,mu,' iadd=',iadd
	  if(iadd.gt.nkmax)then
	   istop=1
	   iaddmore=iadd
	   cycle
	  end if
	  sk(1,iadd)=kxp
	  sk(2,iadd)=kyp
	  sk(3,iadd)=kzp
	  weight(iadd)=weight(i)
     end do
    end do
! Now check to see if we went out of bounds.
	if(istop.eq.1)then
	 write(*,*)' ***************error in kgroup ******************'
	 write(*,*)' ***************error in kgroup ******************'
	 write(*,*)' ***************error in kgroup ******************'
	 write(*,*)' ***************error in kgroup ******************'
	 write(*,*)' ***************error in kgroup ******************'
	 write(*,*)' ***************error in kgroup ******************'
	 write(*,*)' We went out of bounds with k points'
	 write(*,*)' nkmax=max dim=',nkmax,' We went to',iaddmore
	 write(*,*)' Please redimension parameter(nkmax=******) in '
	 write(*,*)' kgroup,kspecial,shells.f'
	 write(*,*)' Sorry -- I must abort!!!!'
	 stop ' BYE'
	 end if
	 numb=numb*ngroup
	 write(*,*)' After the rotation we have a '
	 write(*,*)' new number of k vectors =',numb
	 if(numb.gt.nkmax) stop ' too many points'
	 itest=0
	 if(itest.eq.1)then
         call shells(a1,a2,a3,sk,numb,weight,i1max,i2max,i3max)
	end if
! ================================================================
	write(*,*)' Weed out repeats .....'
! Now some of these repeat. We look for repeats.
	do i=1,numb
	 wait(i)=weight(i)
     irepeat(i)=0
    end do 
	do i=1,numb
	 if(i.eq.1) cycle
     do L=1,i-1
	  xcomp=sk(1,L)-sk(1,i)
	  ycomp=sk(2,L)-sk(2,i)
	  zcomp=sk(3,L)-sk(3,i)
	  if(abs(xcomp)+abs(ycomp)+abs(zcomp).lt.0.00001)then
! This i is already in L
	   irepeat(i)=1
	   wait(L)=wait(L)+weight(i)
	  end if
     end do
    end do
!	write(*,*)' The new icount=',icount
	ic=0
	do i=1,numb
	 if(irepeat(i).eq.1) cycle
     ic=ic+1
	 wait(ic)=wait(i)
	 weight(ic)=wait(ic)
	 sk(1,ic)=sk(1,i)
	 sk(2,ic)=sk(2,i)
	 sk(3,ic)=sk(3,i)
    end do
	numb=ic	
	write(*,*)' After weeding out repeats we have',numb,' kpoints'
	write(*,*)' The new k points are:'
	do i=1,numb
	write(*,712)sk(1,i),sk(2,i),sk(3,i),wait(i)
712	format(' k=',3f9.4,' weight=',f9.4)
    end do
	if(itest.eq.1)then
        call shells(a1,a2,a3,sk,numb,weight,i1max,i2max,i3max)
	end if
! ================================================================
! Now we have the option of finding only those in the irreducible wedge.
!	write(*,*)'  '
!	write(*,*)' Now we have the option of finding only those'
!	write(*,*)' in the irreducible wedge.'
!	write(*,*)' Insert 1 for irr.wedge, 0 for skip this'
!	read(*,*)iwedge
!	write(*,*)' iwedge=',iwedge
!	if(iwedge.eq.1)write(*,*)' Finding those in irr. wedge.'
!	if(iwedge.ne.1)write(*,*)' Not finding k"s in irred. wedge.'
	iwedge=0
	if(iwedge.ne.1)go to 622
! Finding those in irred. wedge.
! Initialize by setting each k to an irreducible k.
	write(*,*)' number of k points=',numb
	do 123 i=1,numb
!	if(mod(i,10).eq.1)write(*,*)i
! irred(i)=1 means it is irreducible (At least initially).
	irred(i)=1
	nirred(i)=0
123	continue
!
! Now go through the list of k-vectors and check whether or not 
! it is irreducible.
	do 510 i=1,numb
	if(mod(i,10).eq.1)write(*,*)i
! ****test
	if(iwrite.eq.1)then
	write(*,*)'  '
	write(*,*)'  '
	write(*,*)'  '
	write(*,*)'  '
	write(*,*)'  '
	write(*,*)' =======================================  '
	write(*,44)i,irred(i),sk(1,i),sk(2,i),sk(3,i)
44	format(' i=',i2,' irred=',i2,' k=',3f9.4)
	write(*,*)' =======================================  '
	end if
! ***testend
	if(irred(i).eq.0)go to 510
! Now rotate kx,ky,kz to kxp,kyp,kzp
	do 600 mu=1,ngroup
!	write(*,*)' mu=',mu
	kxp=pcgrp(1,1,mu)*sk(1,i)+pcgrp(1,2,mu)*sk(2,i)+  &
     &	pcgrp(1,3,mu)*sk(3,i)
!	write(*,890)kxp,pcgrp(mu,1,1),pcgrp(mu,1,2),pcgrp(mu,1,3)
!890	format(' kxp=',f9.4,' pcgrp(mu,1,i)=',3f9.4)
	kyp=pcgrp(2,1,mu)*sk(1,i)+pcgrp(2,2,mu)*sk(2,i)+  &
     &	pcgrp(2,3,mu)*sk(3,i)
	kzp=pcgrp(3,1,mu)*sk(1,i)+pcgrp(3,2,mu)*sk(2,i)+  &
     &	pcgrp(3,3,mu)*sk(3,i)
! ****test
	if(iwrite.eq.1)then
	write(*,45)mu,kxp,kyp,kzp
45	format(' mu=',i3,' kp=',3f9.4)
	end if
! ***testend
! Now search throught the list to see which j k' is.
	jitis=0
	do 500 j=1,numb
! Check it against different zones
	do 412 i1=-1,1
	do 412 i2=-1,1
	do 412 i3=-1,1
	xcomp=sk(1,j)+i1*rlv1(1)+i2*rlv2(1)+i3*rlv3(1)
	ycomp=sk(2,j)+i1*rlv1(2)+i2*rlv2(2)+i3*rlv3(2)
	zcomp=sk(3,j)+i1*rlv1(3)+i2*rlv2(3)+i3*rlv3(3)
	if(abs(kxp-xcomp)+abs(kyp-ycomp)+abs(kzp-zcomp) &
     &	.lt. 0.00001)go to 488
412	continue
	go to 500
488	continue
	jitis=j
	go to 501
500	continue
!
501	continue
! ----
! jitis is the k(j) in which k(i) is rotated into.
	if(jitis.eq.0) stop ' bad jitis'
	if(j.lt.i)write(*,*)' Whoa.... j<i j=',j,' i=',i
	if(j.ne.i)irred(j)=0
! Now we determine how many unique k vectors vector i has been rotated into.
! Check with previous r*k(i) to see if this one has popped up before.
	if(nirred(i).eq.0)go to 42	
	do 22 L=1,nirred(i)
	if(iwhich(L).eq.j)then
! This one came up before
	go to 32
	end if
22	continue
! This one has never come up before
42	continue
	nirred(i)=nirred(i)+1
	iwhich(nirred(i))=j
32	continue
!
!
600	continue
510	continue
! ================================================================
! summary
	write(*,*)' The initial number of k vectors =',numb
!	write(*,*)' irred(i)=0/1: 0 means it is reducible.'
!	write(*,*)' If a kvector is irred, then nirred is the'
!	write(*,*)' number of kvectors which have been reduced to it.'
!	do i=1,numb
!	write(*,410)i,irred(i),nirred(i)
!410	format(' k-vector #',i4,' irred=',i2,' number of reducible=',i3)
!	end do
! Now put irreducible k vectors into skirred(3,nkmax) and weight in wt.
	icount=0
	do 245 i=1,numb
	if(irred(i).eq.0)go to 245
	icount=icount+1
	skirred(1,icount)=sk(1,i)
	skirred(2,icount)=sk(2,i)
	skirred(3,icount)=sk(3,i)
!	wt(icount)=weight*dfloat(nirred(i))
	wt(icount)=wait(i)*dfloat(nirred(i))
245	continue
	numirr=icount
	write(*,*)'  '
	write(*,*)' Finished finding irreducible k-vectors'
	write(*,*)' The irreducible k vectors are: numirr=',numirr
	sum1=0.0
	do 46 i=1,numirr
	write(*,47)skirred(1,i),skirred(2,i),skirred(3,i),wt(i)
	sum1=wt(i)+sum1
47	format(' k=',3f9.4,' weight factor=',3f9.4)
46	continue
	write(*,*)'  '
	write(*,49)sum1
49	format(' Sum of weight factors=',f14.7)
!        write(*,*)'      '
!        write(*,*)' We write them out to tempirr.kpts'
!        write(*,*)'  '
!        open (unit=31,file='tempirr.kpts',status='unknown')
!        write(31,*)numirr
!        write(*,*)' We found numirr=',numirr,' kpoints.'
!        write(*,*)' They are ....'
!        do 441 ik=1,numirr
!        write(*,620)skirred(1,ik),skirred(2,ik),skirred(3,ik),wt(ik)
!        write(31,620)skirred(1,ik),skirred(2,ik),skirred(3,ik),wt(ik)
!620     format(1x,f13.8,1x,f13.8,1x,f13.8,8x,f13.8)
!441     continue
!        close(unit=31)
        write(*,*)'  '
! Now put irreducible k's into k array
	numb=numirr
	do 555 ik=1,numb
	weight(ik)=wt(ik)
	sk(1,ik)=skirred(1,ik)
	sk(2,ik)=skirred(2,ik)
	sk(3,ik)=skirred(3,ik)
555	continue
	if(itest.eq.1)then
        call shells(a1,a2,a3,sk,numb,weight,i1max,i2max,i3max)
	end if
! ================================================================
! Skip the irreducible wedge part.
622	continue
! ================================================================
!        write(*,*)' The number of k vectors is ',numb
!	write(*,*)' successfully finished kgroup.f -------- bye.'
!
	return
	end
! kgroup.f -- Read in the point group data file, and gives the
!             k vectors that we have found which are irreducible.
! =======================================================================
! for questions or futher information:
! o.f. sankey, dept of physics, ASU, 602-965-4334, bitnet: SANKEY@ASUCPS.
! =======================================================================
subroutine kgroupagain(sk,weight,numb,skirred,wt,numirr,pcgrp,  &
     &    npgops,a1,a2,a3,rlv1,rlv2,rlv3,i1max,i2max,i3max)
    implicit none !implicit double precision (a-h,o-z)
! passed variables -----------------------------------------------------        
    integer, intent(in) :: npgops, i1max, i2max, i3max
    integer, intent(out) :: numirr
    integer, intent(inout) :: numb
    real, intent(in) :: a1, a2, a3, rlv1, rlv2, rlv3
    real, intent(inout) :: sk, weight, skirred, wt, pcgrp
! internal variables -----------------------------------------------------
    character(len=50) filename
    integer i,j,nkmax,igrpsz,iwrite,inversion,ngroup,igrp,jgrp, &
     &   ix,iy,istop,mu,iadd,itest,irepeat,l,ic,ik,iskip,numer, &
     &   iwedge,irred,jitis,i1,iaddmore,nirred,i2,i3,iwhich,icount,&
     &   ii,igoto,numbnew,jj,kitis
    real kxp,kyp,kzp,wait,xcomp,ycomp,zcomp, sum1
    parameter (nkmax=100000)
    parameter (igrpsz=96)
    dimension weight(nkmax)
    dimension kitis(igrpsz)
    dimension a1(3),a2(3),a3(3)
    dimension rlv1(3),rlv2(3),rlv3(3)
    dimension sk(3,nkmax)
    dimension skirred(3,nkmax),irred(nkmax),nirred(nkmax)
    dimension wt(nkmax),wait(nkmax),irepeat(nkmax)
    dimension iwhich(igrpsz,nkmax),numer(nkmax)
    dimension pcgrp(3,3,igrpsz)
!	common /lattice/a1,a2,a3,rlv1,rlv2,rlv3
!        common /imaxes/i1max,i2max,i3max
! ================================================================
!	write(*,*)'  '
!	write(*,*)' Welcome to subroutine kgroupagain...'
!	write(*,*)'  '
	write(*,*)' We call kgroupagain, which weeds out kvectors'
	write(*,*)' related by symmetry you input.'
! Terse
	iwrite=0
! Extended
!	iwrite=1
! ================================================================
!$$$	write(*,*)' Insert the file name for the point group operators.'
!$$$	write(*,*)' Note: the first operator must be the identity.'
!$$$	read(*,539)filename
!$$$539	format(a50)
!$$$c	write(*,539)filename
!$$$	write(*,539)filename
!$$$	open(unit=30,file=filename,status='old')
!$$$	read(30,*)ngroup
!$$$	write(*,*)' There are',ngroup,' group operators.'
!$$$	if(ngroup.gt.igrpsz)stop 'too many group operators'
!$$$	do 540 igrp=1,ngroup
!$$$	read(30,*)jgrp
!$$$	if(jgrp.ne.igrp)stop ' something wrong in group file'
!$$$c pcgrp are the group operators in cartesian coordinates.
!$$$c (x',y',z') = pcgrp * (x,y,z) where x,y,z is x*xhat+y*yhat+z*zhat.
!$$$	do 541 i=1,3
!$$$	read(30,*) pcgrp(igrp,i,1),pcgrp(igrp,i,2),pcgrp(igrp,i,3)
!$$$c
!$$$	if(iwrite.eq.1)then
!$$$	write(*,188) pcgrp(igrp,i,1),pcgrp(igrp,i,2),pcgrp(igrp,i,3)
!$$$	end if
!$$$c
!$$$188	format('  ',3f12.8)
!$$$541	continue
!$$$540	continue
!$$$	close(unit=30)
!
! 
! otto
        write(*,*)' For special k-points, we always have'
        write(*,*)' that k ---> -k has the same energy.'
        write(*,*)' DO YOU WISH TO INCLUDE INVERSION???'
        write(*,*)' insert 1 for YES, 0 for NO... usually insert 1'
!        read(*,*)inversion
        inversion=1
        write(*,*)' inversion=',inversion
!        inversion=1
!       write(*,*)'  '

! ========================================================
!       write(*,*)'  '
        ngroup=npgops
        if(inversion.eq.1)then
        write(*,*)' DOUBLE THE GROUP .....'
        do 5440 igrp=1,ngroup
        jgrp=igrp+ngroup
        do ix=1,3
        do iy=1,3
        pcgrp(ix,iy,jgrp)=-pcgrp(ix,iy,igrp)
        end do
        end do
5440    continue
        ngroup=ngroup+ngroup
        write(*,*)' There are',ngroup,' group operators.'
        if(ngroup.gt.igrpsz)stop 'too many group operators'
        end if

! ===========================================================c
! New test. We first get rid of any k vectors which are related by
! reciprocal lattice vectors.
! ============================================================
! We must have the first group operator being the identity.
! Check that first operator is identity.
	do 833 i=1,3
	if(abs(pcgrp(i,i,1)-1.0d0).gt.0.0000001)then
	write(*,*)' the first group operator is not the identity'
	write(*,*)' We assume that group op. # 1=identity and it is not'
	stop ' sorry!!!'
	end if
	do 834 j=1,3
	if(j.eq.i)go to 834
	if(abs(pcgrp(i,j,1)).gt.0.0000001)then
	write(*,*)' the first group operator is not the identity'
	write(*,*)' We assume that group op. # 1=identity and it is not'
	stop ' sorry!!!'
	end if
834	continue
833	continue
! =================================================================
! ================================================================
! Now we have the option of finding only those in the irreducible wedge.
!	write(*,*)'  '
!	write(*,*)' We have the option weeding out those related'
!	write(*,*)' by symmetry.'
!	write(*,*)'  '
!	write(*,*)' Insert 1 for weed out, 0 for skip this'
!	read(*,*)iwedge
	iwedge=1
!	write(*,*)' iwedge=',iwedge
	if(iwedge.eq.1)write(*,*)' Finding those in irr. wedge.'
	if(iwedge.ne.1)write(*,*)' Not finding ks in irred. wedge.'
	if(iwedge.ne.1)go to 622
! Finding those in irred. wedge.
! Initialize by setting each k to an irreducible k.
	write(*,*)' number of k points=',numb
	do 123 i=1,numb
!	if(mod(i,10).eq.1)write(*,*)i
! irred(i)=1 means it is irreducible (At least initially).
	irred(i)=1
	nirred(i)=0
123	continue
!	write(*,*)' A:'
        do 808 i=1,numb
        wait(i)=weight(i)
808     irepeat(i)=0
!	write(*,*)' B:'
! sum weights
	sum1=0.0d0
	do i=1,numb
	sum1=sum1+weight(i)
	end do
!	write(*,*)' C:'
	write(*,*)' Initial sum of weights=',sum1
	if(abs(sum1-1.0d0).gt.0.0001) stop ' bad sum of weights'
! ===========================================================c
! New test. We first get rid of any k vectors which are related by
! reciprocal lattice vectors.
	numbnew=numb
	iskip=1
	if(iskip.eq.1)go to 1111
	write(*,*)'  '
	write(*,*)' numb=',numb,' First we remove those k'
	write(*,*)' vectors which are connected by RLV"s'
	write(*,*)'  '
	do 455 i=1,numb
        irred(i)=1
        nirred(i)=0
	numer(i)=1
455	continue
!	write(*,*)' D:'

! New test. We first get rid of any k vectors which are related by
! reciprocal lattice vectors.
	do 610 i=1,numb
        if(irred(i).eq.0)go to 610
! Now search throught the list to see which j k' is.
	do 912 jj=1,igrpsz
912	kitis(jj)=0
!        jitis=0
	icount=0
        do 2600 j=1,numb

! Check it against different zones
!        do 612 i1=-3,3
!        do 612 i2=-3,3
!        do 612 i3=-3,3

        do 612 i1=-2,2
        do 612 i2=-2,2
        do 612 i3=-2,2
	if(i1.eq.0.and.i2.eq.0.and.i3.eq.0)go to 612
        xcomp=sk(1,j)+i1*rlv1(1)+i2*rlv2(1)+i3*rlv3(1)
        ycomp=sk(2,j)+i1*rlv1(2)+i2*rlv2(2)+i3*rlv3(2)
        zcomp=sk(3,j)+i1*rlv1(3)+i2*rlv2(3)+i3*rlv3(3)
        if(abs(sk(1,i)-xcomp)+abs(sk(2,i)-ycomp)+   &
     &	abs(sk(3,i)-zcomp).lt. 0.001)then
	icount=icount+1
	kitis(icount)=j
	end if
612     continue
!
2600     continue
	if(icount.ne.0)then
	write(*,822)i,icount
822	format(' kvector number',i3,' has several RLV related ',' ks:icount=',i3)
	end if
	do 599 jj=1,icount
	if(kitis(jj).eq.0)go to 599
	write(*,*)' i,kitis=',i,kitis(jj)
	jitis=kitis(jj)
! check dimensions
	if(jitis.gt.nkmax)then
	write(*,*)' jitis=',jitis
	write(*,*)' nkmax (dimenison) is',nkmax
	write(*,*)' set nkmax larger in the parameter statement.'
	stop ' set nkmax larger in the parameter statement'
	end if
! jitis is the k(j) in which k(i) is rotated into.
        if(j.lt.i)write(*,*)'  RLV check! Whoa.... j<i j=',j,' i=',i
	irred(jitis)=0
	numer(i)=numer(i)+1
	weight(i) =weight(i)+weight(jitis)
599	continue
610	continue
!	write(*,*)' E:'
! ==================================
! ****error found
!	newnumb=0
	numbnew=0	
! ***test
!	write(*,*)'  **test**'
!	write(*,*)' numb=',numb
! ***testend
	do i=1,numb
!
	if(irred(i).eq.1)then
	numbnew=numbnew+1
! check dimensions
        if(numbnew.gt.nkmax)then
        write(*,*)' numbnew=',numbnew
        write(*,*)' nkmax (dimenison) is',nkmax
        write(*,*)' set nkmax(numbnew)larger parameter statement.'
        stop ' set nkmax (numbnew)larger parameter statement'
        end if
!
	sk(1,numbnew)=sk(1,i)
	sk(2,numbnew)=sk(2,i)
	sk(3,numbnew)=sk(3,i)
	weight(numbnew)=weight(i)
	end if
!
	end do
! ================================
!	write(*,*)' F:'
!	write(*,*)'  '
1111	continue
	numb=numbnew
!	write(*,*)'  '
!	write(*,*)' numb=',numb,' After we remove those k'
!	write(*,*)' vectors which were connected by RLV"s'
!	write(*,*)'  '
	write(*,*)' The new number after killing those'
	write(*,*)' conncected by reciprocal lattice vectors'
	write(*,*)' is', numb
!	write(*,*)'  '
!	write(*,*)'   They are:'
!	do ik=1,numb
!        write(*,20)sk(1,ik),sk(2,ik),sk(3,ik),weight(ik)
!20      format(1x,f13.8,1x,f13.8,1x,f13.8,8x,f13.8)
!	end do
!	write(*,*)' G:'
! sum weights
        sum1=0.0d0
        do i=1,numb
        sum1=sum1+weight(i)
        end do
!        write(*,*)'  Intermediate sum of weights=',sum1
        if(abs(sum1-1.0d0).gt.0.0001) stop ' bad sum of weights'

        do 554 i=1,numb
        irred(i)=1
        nirred(i)=0
        numer(i)=1
554     continue
! ============================================================   
	do i=1,numb
	iwhich(1,i)=i
	do ii=1,igrpsz
	iwhich(ii,i)=0
	end do
	end do
	
!
! Now go through the list of k-vectors and check whether or not 
! it is irreducible.
	do 510 i=1,numb
!	if(mod(i,10).eq.1)write(*,*)i,' Must go to ',numb
! ****test
	if(iwrite.eq.1)then
	write(*,*)'  '
	write(*,*)'  '
	write(*,*)'  '
	write(*,*)'  '
	write(*,*)'  '
	write(*,*)' =======================================  '
	write(*,44)i,irred(i),sk(1,i),sk(2,i),sk(3,i)
44	format(' i=',i2,' irred=',i2,' k=',3f9.4)
	write(*,*)' =======================================  '
	end if
! ***testend
	if(irred(i).eq.0)go to 510
! Now rotate kx,ky,kz to kxp,kyp,kzp
	do 600 mu=1,ngroup
!	write(*,*)' mu=',mu
	kxp=pcgrp(1,1,mu)*sk(1,i)+pcgrp(1,2,mu)*sk(2,i)+  &
     &	pcgrp(1,3,mu)*sk(3,i)
!	write(*,890)kxp,pcgrp(mu,1,1),pcgrp(mu,1,2),pcgrp(mu,1,3)
!890	format(' kxp=',f9.4,' pcgrp(mu,1,i)=',3f9.4)
	kyp=pcgrp(2,1,mu)*sk(1,i)+pcgrp(2,2,mu)*sk(2,i)+  &
     &	pcgrp(2,3,mu)*sk(3,i)
	kzp=pcgrp(3,1,mu)*sk(1,i)+pcgrp(3,2,mu)*sk(2,i)+  &
     &	pcgrp(3,3,mu)*sk(3,i)
! ****test
	if(iwrite.eq.1)then
	write(*,45)mu,kxp,kyp,kzp
45	format(' mu=',i3,' kp=',3f9.4)
	end if
! ***testend
! Now search throught the list to see which j k' is.
	jitis=0
	do 500 j=1,numb
! Check it against different zones
!	do 412 i1=-3,3
!	do 412 i2=-3,3
!	do 412 i3=-3,3

	do 412 i1=-2,2
	do 412 i2=-2,2
	do 412 i3=-2,2
!	do 412 i1=-1,1
!	do 412 i2=-1,1
!	do 412 i3=-1,1
	xcomp=sk(1,j)+i1*rlv1(1)+i2*rlv2(1)+i3*rlv3(1)
	ycomp=sk(2,j)+i1*rlv1(2)+i2*rlv2(2)+i3*rlv3(2)
	zcomp=sk(3,j)+i1*rlv1(3)+i2*rlv2(3)+i3*rlv3(3)
	if(abs(kxp-xcomp)+abs(kyp-ycomp)+abs(kzp-zcomp)  &
     &	.lt. 0.001)go to 488
!     1	lt. 0.00001)go to 488
412	continue
	go to 500
488	continue
	jitis=j
	if(j.gt.i)	go to 501
500	continue
!
501	continue
! otto
! ----
! jitis is the k(j) in which k(i) is rotated into.
	if(jitis.eq.0)then
	write(*,*)' Your symmetry is probably wrong.'
	stop ' bad jitis'
	end if
	if(jitis.lt.i)write(*,*)' Whoa.... j<i j=',j,' i=',i
	if(jitis.ne.i)irred(j)=0
! ****test
!	if(mu.eq.1)write(*,*)'  '
!	write(*,822)i,mu,jitis,j,irred(j)
!822	format(' i=',i3,'    mu=',i2,' jitis,j=',2i3,' irred(j)=',i2)
!	write(*,9833)(irred(ll),ll=1,numb)
!9833	format(' irred=',14i3)
! ***testend
! Now we determine how many unique k vectors vector i has been rotated into.
! Check with previous r*k(i) to see if this one has popped up before.
	if(nirred(i).eq.0)go to 42	
	do 22 L=1,nirred(i)
	if(iwhich(L,i).eq.jitis)then
! This one came up before
	go to 32
	end if
22	continue
! This one has never come up before
42	continue
	nirred(i)=nirred(i)+1
	iwhich(nirred(i),i)=jitis
32	continue
!
!
600	continue
510	continue
!	write(*,*)' '
!        write(*,*)'  '
!        do ik=1,numb
!        write(*,2001)sk(1,ik),sk(2,ik),sk(3,ik),weight(ik),
!     1	irred(ik),nirred(ik)
!2001    format(1x,f13.8,1x,f13.8,1x,f13.8,' wt',f13.8,' irr,nirr=',2i3)
!        end do
! sum weights
        sum1=0.0d0
        do i=1,numb
	wt(i)=0.0d0
        sum1=sum1+weight(i)
        end do
        write(*,*)' Final sum of weights=',sum1
        if(abs(sum1-1.0d0).gt.0.0001) stop ' bad sum of weights'
        do 2808 i=1,numb
2808    wait(i)=weight(i)
! ================================================================
! summary
!	write(*,*)'  '
!	write(*,*)'  '
!	write(*,*)' The initial number of k vectors =',numb
!	write(*,*)' irred(i)=0/1: 0 means it is reducible.'
!	write(*,*)' If a kvector is irred, then nirred is the'
!	write(*,*)' number of kvectors which have been reduced to it.'
!	do i=1,numb
!	write(*,410)i,irred(i),nirred(i)
!410	format(' k-vector #',i4,' irred=',i2,' number of reducible=',i3)
!	end do
! Now put irreducible k vectors into skirred(3,nkmax) and weight in wt.
	icount=0
	do 245 i=1,numb
	if(irred(i).eq.0)go to 245
	icount=icount+1
	skirred(1,icount)=sk(1,i)
	skirred(2,icount)=sk(2,i)
	skirred(3,icount)=sk(3,i)
!	wt(icount)=weight*dfloat(nirred(i))
	do L=1,nirred(i)
!	write(*,*)' i,L,iwhich=',i,L,iwhich(L,i)
	igoto=iwhich(L,i)
!	wt(icount)=wait(i)*dfloat(nirred(i))
	wt(icount)=wait(igoto)+wt(icount)
!	write(*,*)' igoto',igoto,' weight=',wait(igoto)
	end do
245	continue
	numirr=icount
!	write(*,*)'  '
!	write(*,*)' Finished finding irreducible k-vectors'
!	write(*,*)' The irreducible k vectors are: numirr=',numirr
	sum1=0.0
	do 46 i=1,numirr
!	write(*,47)skirred(1,i),skirred(2,i),skirred(3,i),wt(i)
	sum1=wt(i)+sum1
47	format(' k=',3f9.4,' weight factor=',3f9.4)
46	continue
!	write(*,*)'  '
	write(*,49)sum1
49	format(' Sum of weight factors=',f14.7)
!        write(*,*)'  '
! Now put irreducible k's into k array
	numb=numirr
	do 555 ik=1,numb
	weight(ik)=wt(ik)
	sk(1,ik)=skirred(1,ik)
	sk(2,ik)=skirred(2,ik)
	sk(3,ik)=skirred(3,ik)
555	continue
	itest=0
	if(itest.eq.1)then
        call shells(a1,a2,a3,sk,numb,weight,i1max,i2max,i3max)
	end if
! ================================================================
! Skip the irreducible wedge part.
622	continue
! ================================================================
!        write(*,*)' The number of k vectors is ',numb
!	write(*,*)' successfully finished kgroupagain.f -------- bye.'
!	itest=1
!	if(itest.eq.1) stop
!
! sum weights
        sum1=0.0d0
        do i=1,numb
        sum1=sum1+weight(i)
        end do
!        write(*,*)' Leaving with a sum of weights=',sum1
        if(abs(sum1-1.0d0).gt.0.0001) stop ' bad sum of weights'

	return
	end
! ===================================================================
! shells.f This program takes the three direct lattice vectors, a1 a2 a3,
! and deternies shells of lattice vectors.
! 
! It then computes sum(k) weight * cdexp(ai*K*L).
! This sum should be zero for each L in the shell, if the special points
! good enough for this shell. As the number of k points gets large, the
! size of the shell in which the sum is not zero gets large. Of course, the
! first shell, L=(0,0,0), the sum is 1.
!
! Otto F. Sankey, Dept. of Physics, Arizona State University, Tempe, AZ 85287
! 602 965-4334 sankey@edison.la.asu.edu
!
! ===================================================================
subroutine shells(a1,a2,a3,sk,numbk,weight,i1max,i2max,i3max)
    implicit none !double precision (a-h,o-z)
! passed variables --------------------------------------------------
    integer, intent(in) :: numbk, i1max, i2max, i3max
    real, intent(in) :: a1, a2, a3, weight
! internal variables ------------------------------------------------
    integer nkmax,maxshell,Lmax, &
     &   Lrot,iyes,numbshell,lattshell,nshellmax,icount,i1,i2,i3, &
     &   i,nshell,itrue,mshell,ishell,ik
    real sk,r,rvec,shellsize,size1,toll,rmag, &
     &   small,rshell,dot,rsum
! The maximum number of k points
    parameter (nkmax=50000)
! The maximum number of shells
	parameter (maxshell=2222)
! The maximum number of lattice vectors in all shells
!	parameter (Lmax=20001)
	parameter (Lmax=2001)
! The maximum number of L-vectors in a shell
	parameter (Lrot=2122)
! ===================================================================
	complex :: ai,sum1  ! used to be double precision
	dimension weight(nkmax)
        dimension a1(3),a2(3),a3(3)
	dimension r(3),rvec(3)
	dimension shellsize(maxshell),size1(Lmax),iyes(Lmax)
	dimension numbshell(maxshell),lattshell(maxshell,Lrot,3)
        dimension sk(3,nkmax)
	character(len=4) mess4
! ===================================================================
	write(*,*)'  '
	write(*,*)' Welcome to shells.f --- We will test whether'
	write(*,*)' or not, sum1(k) e(i*K*L) is zero or not for'
	write(*,*)' L values in different shells.'
	write(*,*)'  '
	if(numbk.gt.nkmax)stop ' too many k points!!'
!	i1max=30
!	i2max=30
!	i3max=0
	nshellmax=maxshell
! Now just loop over all lattice vectors and put the size of each one
! into size(icount).
	icount=0
	toll=0.001
	do 400 i1=-i1max,i1max 
	do 400 i2=-i2max,i2max
	do 400 i3=-i3max,i3max
	rmag=0.0
	icount=icount+1
	iyes(icount)=0
	do i=1,3
	r(i)=i1*a1(i)+i2*a2(i)+i3*a3(i)
	rmag=rmag+r(i)**2
	end do
	rmag=sqrt(rmag)
	size1(icount)=rmag
!	write(*,*)' rmag=',rmag
400	continue
	if(icount.gt.Lmax)then
	write(*,*)' icount=',icount
	write(*,*)' dimension=',Lmax
	write(*,*)' Sorry -- you must redimension'
	end if
!
! ===================================================================
! Now define the smallest size as the first shell, the second size as the
! second shell and so on.
	small=10000.
! Loop over all shells
	itrue=0
	do 61 mshell=1,nshellmax
!
	do 55 i=1,icount	
! Skip this one if we've already counted it.
	if(iyes(i).eq.1)go to 55
! If its the smallest so far, set small=current size.
	if(size1(i).lt.small+toll)then
	small=size1(i)
	iyes(i)=1
	end if
!
55	continue
!
	if(small.lt.10000.-1.)itrue=itrue+1
	shellsize(mshell)=small
	small=10000.
61	continue
! We found less shells than nshellmax. Set upper limit of shells to correct
! value
	nshellmax=itrue
	write(*,*)' Number of shells found=',itrue
! ===================================================================
! Now put lattice vectors into each shell.
	do 100 ishell=1,nshellmax
	rshell=shellsize(ishell)
	nshell=0
!	
        do 101 i1=-i1max,i1max
        do 101 i2=-i2max,i2max
        do 101 i3=-i3max,i3max
        rmag=0.0
        do 711 i=1,3
        r(i)=i1*a1(i)+i2*a2(i)+i3*a3(i)
711     rmag=rmag+r(i)**2
        rmag=sqrt(rmag)
	if(abs(rshell-rmag).lt.toll)then
	nshell=nshell+1
! lattshell ---> shell number ishell has nshell elements and i1,i2,i3 are the
! components.
	lattshell(ishell,nshell,1)=i1
	lattshell(ishell,nshell,2)=i2
	lattshell(ishell,nshell,3)=i3
	end if
101	continue
	numbshell(ishell)=nshell
	if(nshell.gt.Lrot) stop ' nshell too large, Lrot too small'
100	continue
!
! ===================================================================
! Now test sum k exp(i*K*L)
	ai=dcmplx(0.0d0,1.0d0)
	write(*,*)' Numbk in shells=',numbk
!	do i=1,numbk
!	write(*,*)' i=',i,' weight=',weight(i)
!	end do
!
	do 888 ishell=1,nshellmax
	sum1=0.0
!	Bigsum=0.0
! Compute the sum for each L, and take the largest one as an indication 
! of the error.
	do 45 nshell=1,numbshell(ishell)
!	nshell=1
	i1=lattshell(ishell,nshell,1)
	i2=lattshell(ishell,nshell,2)
	i3=lattshell(ishell,nshell,3)
	rvec(1)=i1*a1(1)+i2*a2(1)+i3*a3(1)
	rvec(2)=i1*a1(2)+i2*a2(2)+i3*a3(2)
	rvec(3)=i1*a1(3)+i2*a2(3)+i3*a3(3)
	do 46 ik=1,numbk
	dot=rvec(1)*sk(1,ik)+rvec(2)*sk(2,ik)+rvec(3)*sk(3,ik)
	sum1=sum1+weight(ik)*cexp(ai*dot)
!	write(*,*)' ik=',ik,' dot=',dot,' sum=',sum1
46	continue
!	rsum=cdabs(sum1)
!	if(rsum.gt.Bigsum)Bigsum=rsum
45	continue
	rsum=cabs(sum1)
	mess4='FAIL'
!	if(Bigsum.lt.1.0d-6)mess4='    '
	if(rsum.lt.1.0d-6)mess4='    '
!	write(*,22)ishell,Bigsum,mess4
	write(*,22)ishell,rsum,mess4
22	format(' shell=',i5,' sum=',1pe16.6,'  ',a4)
888	continue	
! ===================================================================
	return
	end	
