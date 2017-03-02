! copyright info:
!
!                             @Copyright 1998
!                          Fireball2000 Committee
!
! ASU - Otto F. Sankey
!       Kevin Schmidt
!       Jian Jun Dong
!       John Tomfohr
!       Gary B. Adams
!
! Motorola - Alex A. Demkov
!            Jun Wang
!
! University of Regensburg - Juergen Fritsch
!
! University of Utah - James P. Lewis
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu
 
!
! RESTRICTED RIGHTS LEGEND
!
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
! DtrescentrosG.f
!
! Adapted from Dtrescentros. This gets the same three-center matrix elements 
! and derivatives gotten by Dtrescentros but uses the gaussian fits.
!
! Program Description for Dtrescentros
! ===========================================================================
!       This routine is like trescentrosG.f, except that it calculates the
! forces.
!
!       This subroutine calculates the (three-center) matrix elements (mu,nu).
! It calculates the matrix elements of a bondcharge with a NA interaction,
! using data which was stored in - threecint(nz,nx,kbc,jxx,ithet,ind1,ind2,ind3)
! This subroutine calculates the matrix elements molecular and crystal
! coordinates of a bondcharge with the neutral atom elsewhere (i.e. neither an
! "ontop" nor an atom-neutral atom).
!
!
!           X  ^                                * (NA)
!              !                             +
! X-Z PLANE    !                          +
!              !                      +   THETA               -----> Z
!              O ---------------------------------------------O
!            1 END OF BC                                OTHER END OF BC
!               MU                                          NU
!
! The KEY for isorp for x/c 3-center. Evaluate the spherically averaged wave
! functions for all orbitals of the two atoms.
! (x y z) = charge on (1 2 3).       1=Left  2=right  3=NA.
!
! isopr = 0    Neutral  (0 0 0)
!         1             (- 0 0)
!         2             (+ 0 0)
!         3             (0 - 0)
!         4             (0 + 0)
!         5             (0 0 -)
!         6             (0 0 +)
!
!
! ymin,ymax,numy: bond charge  distances grid
! xmin,xmax,numx: neutral atom distances grid
!
! Here are the data file characteristics:
!
! ===========================================================================
! Code originally written by Jose Ortega and Juergen Fritsch
 
! Code rewritten by:
! James P. Lewis
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E. RM Dock
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
!
! Code re-rewritten for gaussians by:
! John Tomfohr
! john.tomfohr@asu.edu
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine DtrescentrosG_VNA_SH (isorp, in1, in2, indna, x, y,  &
     &                                   cost, eps, depsA, depsB, rhat, &
     &                                   sighat, bcnax, f3naXa, f3naXb, &
     &                                   f3naXc, rcutoff)
        use dimensions
        use interactions
        use gaussG
        use constants_fireball
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: in1
        integer, intent (in) :: in2
        integer, intent (in) :: indna
        integer, intent (in) :: isorp
 
        real, intent (in) :: cost
        real, intent (in) :: x
        real, intent (in) :: y

        real, intent (in), dimension (3, 3) :: eps
        real, intent (in), dimension (3, 3, 3) :: depsA
        real, intent (in), dimension (3, 3, 3) :: depsB
        real, intent (in), dimension (3) :: rhat
        real, intent (in), dimension (3) :: sighat

! JKT step. added rcutoff
        real, intent (in), dimension (0:nspec_max-1, nsh_max) :: rcutoff
 
! Output
        real, intent (out), dimension (numorb_max, numorb_max) :: bcnax
        real, intent (out), dimension (3, numorb_max, numorb_max) :: f3naXa
        real, intent (out), dimension (3, numorb_max, numorb_max) :: f3naXb
        real, intent (out), dimension (3, numorb_max, numorb_max) :: f3naXc
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer inu
        integer ix

        integer iparam
 
        real argument

        real amt
        real bmt

        integer iindex
        integer jindex
        integer l1
        integer l2

        integer ish
        integer jsh
        integer im
        integer jm

        real rmatel
        real rmatel1
        real rmatel2

        real rmately1
        real rmately2
        real rmately

        real rmatelx31
        real rmatelx32
        real rmatelx3

        real rmatelz31
        real rmatelz32
        real rmatelz3
 
        real, dimension (numorb_max, numorb_max) :: bcnam
        real, dimension (numorb_max, numorb_max) :: dpbcnam
        real, dimension (numorb_max, numorb_max) :: dxbcnam
        real, dimension (numorb_max, numorb_max) :: dybcnam
        real, dimension (3, numorb_max, numorb_max) :: f3naMa
        real, dimension (3, numorb_max, numorb_max) :: f3naMb

! JKT step
        real rc1
        real rc2
        real rc3

        real r12
        real r23
        real r31

        real ptemp12
        real ptemp23
        real ptemp31

        real fermi12
        real fermi23
        real fermi31
        real sint

        real Dfermix3
        real Dfermiz3
        real Dfermiy1
        real Dfermiy2
        real Dfermiy

! JKT step for Dtresgaussians only
        real fermi123
        real x1
        real y1
        real z1

        real x2
        real y2
        real z2

        real x3
        real y3
        real z3

! Procedure
! ===========================================================================
! Initialize bcnax and forces.
        do inu = 1, num_orb(in2)
           do imu = 1, num_orb(in1)
              bcnax(imu,inu) = 0.0d0
           f3naXa(1,imu,inu) = 0.0d0
           f3naXa(2,imu,inu) = 0.0d0
           f3naXa(3,imu,inu) = 0.0d0
           f3naXb(1,imu,inu) = 0.0d0
           f3naXb(2,imu,inu) = 0.0d0
           f3naXb(3,imu,inu) = 0.0d0
           enddo
        enddo

! JKT step. why is rcutoff indexed this way?! (0:nspec-1)?!!!
! JKT step  get rc3...
        rc3=0.
        do ish=1,nssh(indna)
        if (rcutoff(indna-1,ish).gt.rc3) then
        rc3=rcutoff(indna-1,ish)
        end if
        end do
!
! JKT step. get r12,r23,r31
        sint=sqrt(abs(1.-cost**2))
        r12=y
        r23=sqrt((y/2.-x*cost)**2+(x*sint)**2)
        r31=sqrt((y/2.+x*cost)**2+(x*sint)**2)

! JKT step. get ptemps. (pseudo-temperatures)
        ptemp12=ptemp*r12
        ptemp23=ptemp*r23
        ptemp31=ptemp*r31

! get bcnam just as in tresgaussians
        iindex=0
        do ish=1,nssh(in1)
        l1=lssh(ish,in1)
! JKT step rc1
        rc1=rcutoff(in1-1,ish)
        do im=-l1,l1
        iindex=iindex+1
        jindex=0
        do jsh=1,nssh(in2)
        l2=lssh(jsh,in2)
! JKT step rc2
        rc2=rcutoff(in2-1,jsh)
        do jm=-l2,l2
        jindex=jindex+1

! if only one of im or jm is .lt. 0 then rmatel=0. by symmetry
        if ((im.lt.0 .and. jm.ge.0) .or. (jm.lt.0 .and. im.ge.0)) then
        rmatel1=0.
        rmatel2=0.
        rmatel=0.
! JKT step added else if
        else if (r12.gt.(rc1+rc2) .or. r23.gt.(rc2+rc3) .or.            &
     &           r31.gt.(rc3+rc1)) then
        rmatel1=0.
        rmatel2=0.
        rmatel=0.
        else
        call gelementsG_VNA_SH  (isorp, in1, in2, indna, l1, l2, ish,   &
     &                           jsh, im, jm, x, y, cost, rmatel1)

        call gelementsG_VNA_SH_RNA  (isorp, in1, in2, indna, l1, l2,    &
     &                               ish, jsh, im, jm, y, rmatel2)

        rmatel = rmatel1 + rmatel2
  
        fermi12=2./(exp((r12-(rc1+rc2))/ptemp12)+1)-1.
        fermi23=2./(exp((r23-(rc2+rc3))/ptemp23)+1)-1.
        fermi31=2./(exp((r31-(rc3+rc1))/ptemp31)+1)-1.
        rmatel=rmatel*fermi12*fermi23*fermi31
        end if

        bcnam(iindex,jindex)=rmatel
! im, jm
        end do
        end do
! ish, jsh
        end do
        end do
! =============================  got bcnam  =================================

! Rotate bcnam into crystal-coordinates: bcnam => bcnax
        call rotate_fb (in1, in2, eps, bcnam, bcnax)
 
! ***************************************************************************
! Now consider the components of the different forces which is determined
! by whether or not the force is with respect to atom 3 or atom 1.

! JKT step... need these variables:
        x1 = 0.
        y1 = 0.
        z1 = y/2.
        x2 = 0.
        y2 = 0.
        z2 = y/2.
        x3 = x*sqrt(abs(1.-cost**2))
        y3 = 0.
        z3 = x*cost


        iindex=0
        do ish=1,nssh(in1)
        l1=lssh(ish,in1)
! JKT step
        rc1=rcutoff(in1-1,ish)
        do im=-l1,l1
        iindex=iindex+1
        jindex=0
        do jsh=1,nssh(in2)
        l2=lssh(jsh,in2)
! JKT step rc2
        rc2=rcutoff(in2-1,jsh)
        do jm=-l2,l2
        jindex=jindex+1

! The iparam says whether the derivatives
! are w/respect to the y (the distance between centers 1 and 2) and 
! x3 or z3 (coordinates of na potential...)
!    iparam = 1 ---> y
!    iparam = 2 ---> x3 
!    iparam = 3 ---> z3

! if only one of im or jm is .lt. 0 then everything's 0. by symmetry
        if ((im.lt.0 .and. jm.ge.0) .or. (jm.lt.0 .and. im.ge.0)) then
        rmately1=0.
        rmatelx31=0.
        rmatelz31=0.

        rmately2=0.
        rmatelx32=0.
        rmatelz32=0.

        rmately=0.
        rmatelx3=0.
        rmatelz3=0.
! JKT step added else if
        else if (r12.gt.(rc1+rc2) .or. r23.gt.(rc2+rc3) .or.            &
     &           r31.gt.(rc3+rc1)) then
        rmately1=0.
        rmatelx31=0.
        rmatelz31=0.

        rmately2=0.
        rmatelx32=0.
        rmatelz32=0.

        rmately=0.
        rmatelx3=0.
        rmatelz3=0.
        else

        iparam=1
        call DgelementsG_VNA_SH  (l1,l2,ish,jsh,im,jm,in1,in2,indna,x,  &
     &                            y,cost,rmately1,iparam,isorp)

        call DgelementsG_VNA_SH_RNA (l1,l2,ish,jsh,im,jm,in1,in2,indna, &
     &                               y,rmately2,isorp)

        rmately = rmately1 + rmately2

        iparam=2
        call DgelementsG_VNA_SH  (l1,l2,ish,jsh,im,jm,in1,in2,indna,x,  &
     &                            y,cost,rmatelx3,iparam,isorp)

        iparam=3
        call DgelementsG_VNA_SH  (l1,l2,ish,jsh,im,jm,in1,in2,indna,x,  &
     &                             y,cost,rmatelz3,iparam,isorp)

! JKT step. If you really want to do it right, uncomment the lines
! Between START and STOP. Personally, I don't think it's worth the 
! extra computing time. I doubt you will notice a difference...
!
! ***START*** 
!         fermi12=2./(exp((r12-(rc1+rc2))/ptemp12)+1)-1.
!         fermi23=2./(exp((r23-(rc2+rc3))/ptemp23)+1)-1.
!         fermi31=2./(exp((r31-(rc3+rc1))/ptemp31)+1)-1.
! 
! 	Dfermix3=fermi12*( &
!      &    fermi23* &
!      &    (-0.5*(fermi31+1.)**2* &
!      &      (x3-x1)*exp((r31-(rc3+rc1))/ptemp31)/ptemp31)/r31 &
!      &        + &
!      &    (-0.5*(fermi23+1.)**2* &
!      &      (x3-x2)*exp((r23-(rc2+rc3))/ptemp23)/ptemp23)/r23 * &
!      &    fermi31)
! 	Dfermiz3=fermi12*( &
!      &    fermi23* &
!      &    (-0.5*(fermi31+1.)**2* &
!      &      (z3-z1)*exp((r31-(rc3+rc1))/ptemp31)/ptemp31)/r31 &
!      &        + &
!      &    (-0.5*(fermi23+1.)**2* &
!      &      (z3-z2)*exp((r23-(rc2+rc3))/ptemp23)/ptemp23)/r23 * &
!      &    fermi31)
! 
! 	Dfermiy1=fermi23*( &
!      &    fermi31* &
!      &    (-0.5*(fermi12+1.)**2* &
!      &      (z1-z2)*exp((r12-(rc1+rc2))/ptemp12)/ptemp12)/r12 &
!      &        + &
!      &    (-0.5*(fermi31+1.)**2* &
!      &      (z1-z3)*exp((r31-(rc3+rc1))/ptemp31)/ptemp31)/r31 * &
!      &    fermi12) 
! 	Dfermiy2=fermi31*( &
!      &    fermi12* &
!      &    (-0.5*(fermi23+1.)**2* &
!      &      (z2-z3)*exp((r23-(rc2+rc3))/ptemp23)/ptemp23)/r23 &
!      &        + &
!      &    (-0.5*(fermi12+1.)**2* &
!      &      (z2-z1)*exp((r12-(rc1+rc2))/ptemp12)/ptemp12)/r12 * &
!      &    fermi23)
! 	Dfermiy=-0.5*Dfermiy1+0.5*Dfermiy2
! 
! 	fermi123=fermi12*fermi23*fermi31
! 
! 	rmately=rmately*fermi123 + bcnam(iindex,jindex)*Dfermiy
! 	rmatelx3=rmatelx3*fermi123 + bcnam(iindex,jindex)*Dfermix3
! 	rmatelz3=rmatelz3*fermi123 + bcnam(iindex,jindex)*Dfermiz3
! ***STOP***

        end if

        dybcnam(iindex,jindex)=rmately
        dxbcnam(iindex,jindex)=sint*rmatelx3 + cost*rmatelz3
!
        if (abs(sint).lt.1.E-4) then
        dpbcnam(iindex,jindex)= x*rmatelz3
        else
        dpbcnam(iindex,jindex)=-x*cost/sint*rmatelx3 + x*rmatelz3
        end if

! im, jm
        end do
        end do
! ish, jsh
        end do
        end do

! Now get f3naMa and f3naMb

! ***************************************************************************
! Now consider the components of the different forces which is determined
! by whether or not the force is with respect to atom 3 or atom 1.
        do ix = 1, 3

! The first piece will be the force with respect to atom 3.
         if (x .gt. 1.0d-5) then
          amt = (sighat(ix) - cost*rhat(ix))/x
         else
          amt = 0.0d0
         end if

        do inu=1,num_orb(in2)
         do imu=1,num_orb(in1)
          f3naMa(ix,imu,inu) = rhat(ix)*dxbcnam(imu,inu) +              &
     &                              amt*dpbcnam(imu,inu)
         end do
        end do

! The second piece will be the force with respect to atom 1.
         bmt = (cost*sighat(ix) - rhat(ix))/y
        do inu=1,num_orb(in2)
         do imu=1,num_orb(in1)
          f3naMb(ix,imu,inu) = - sighat(ix)*dybcnam(imu,inu) +          &
     &                bmt*dpbcnam(imu,inu) - f3naMa(ix,imu,inu)/2.
         end do
        end do
! ix
        end do
! ***************************************************************************
! Convert to Crystal Coordinates
! ***************************************************************************
! The call to rotated does the rotations to crystal coordinates of these
! force things.
!
! For example:
! Suppose we have f_A(3,mu,nu), which is d/dratm M(mu,nu) where M(mu,nu)
! is in molecular. To transform M(mu,nu) to crystal, we need Udag * M * U.
! Therefore, f_A(3,mu,nu)[CRYSTAL] = (d/dratm Udag) * M * U
!                                   + Udag * M * (d/dratm U)
!                                   + Udag * f_A * U.
!
! So, to use this baby, put in deps3c (deps/dr1, deps/dr2, deps/dratm),
! and f_A and M.
!
! NOTE: rotated works on the assumption that we are adding derivatives,
! NOT forces. So f3naMa,... etc. MUST not yet be forcelike.
! We do the - sign for forces at the end.
! ***************************************************************************
! Force on the neutral atom with respect to atom 3 (f3naMa).
        call rotated (in1, in2, eps, depsA, bcnam, f3naMa, f3naXa)
 
! Force on the neutral atom with respect to atom 1 (f3naMb).
        call rotated (in1, in2, eps, depsB, bcnam, f3naMb, f3naXb)
 
! Make things force-like.
! Now determine the force f3naXc which is found from Newton's Laws.
        do inu = 1, num_orb(in2)
         do imu = 1, num_orb(in1)
          f3naXa(:,imu,inu) = - f3naXa(:,imu,inu)
          f3naXb(:,imu,inu) = - f3naXb(:,imu,inu)
          f3naXc(:,imu,inu) = - f3naXa(:,imu,inu) - f3naXb(:,imu,inu)
         end do
        end do 

118     format(2(2x, i4), 3f15.6)

! Format Statements
! ===========================================================================
 
        return
        end
