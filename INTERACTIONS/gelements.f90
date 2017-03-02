! copyright info:
!
!                             @Copyright 2003
!                            Fireball Committee
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
 
! gelements.f90
! Program Description
! ===========================================================================
! This routine gets the matrix element - <psi|V|psi>, where V is the potential
! of the neutral atom or exchange-correlation.
! Example: l1 = 1, imu = 1 and l2 = 1, jmu = 0 means get <p_x|Vna|p_z>.
! (p_x on center 1 and p_z on center 2...)
!
! ===========================================================================
! Original code from John Tomfohr
! Code rewritten by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine gelements (in1, in2, indna, l1, l2, issh, jssh, imu, jmu, &
     &                        x, y, cost, rmatel, ish00)
        use dimensions
        use interactions
        use gauss
        use constants_fireball
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: imu
        integer, intent (in) :: in1
        integer, intent (in) :: in2
        integer, intent (in) :: indna
        integer, intent (in) :: ish00
        integer, intent (in) :: issh
        integer, intent (in) :: jmu
        integer, intent (in) :: jssh
        integer, intent (in) :: l1
        integer, intent (in) :: l2
 
        real, intent (in) :: cost
        real, intent (in) :: x
        real, intent (in) :: y

! Output
        real, intent (out) :: rmatel

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        integer maxpoly

! for sp3 maxpoly = 4 :  1,x,y,z
! for sp3d5 maxpoly = 10 : 1,x,y,z,x**2,y**2,z**2,xy,xz,yz
! f-orbitals not available yet...
        parameter(maxpoly=10)

! local

        integer ial1,ial2,ial3
        real rmu,lilm,bigM,d12,dcm3,bigQ,gM,xcm,ycm,zcm
        real al1,al2,al3
        real coef1,coef2,coef3
        real factor
        real c2,c4

        real x1,y1,z1,x2,y2,z2,x3,y3,z3
        
        real, dimension (maxpoly) :: poly1,poly2
        integer ipoly

! initialize poly1 and poly2 to zero
        do ipoly=1,maxpoly
        poly1(ipoly)=0.
        poly2(ipoly)=0.
        end do

! ===============
! Start procedure
! ===============

        d12=y

        x1=0.
        y1=0.
        z1=0.
        x2=0.
        y2=0.
        z2=y
        x3=x*sqrt(abs(1.-cost**2))
        y3=0.
        z3=y/2.+x*cost

        rmatel=0.

! alpha1
        do ial1=1,nalpha(issh,in1)
        al1=alpha(ial1,issh,in1)
        coef1=gcoefficients(ial1,issh,in1)

! alpha2
        do ial2=1,nalpha(jssh,in2)
        al2=alpha(ial2,jssh,in2)
        coef2=gcoefficients(ial2,jssh,in2)

        rmu=al1*al2/(al1+al2)
        lilm=al1+al2
        dcm3=(y/2.-al2*y/(al1+al2)+x*cost)**2+x**2*(1.-cost**2)
        dcm3=sqrt(abs(dcm3))

! alpha3
        do ial3=1,nalpha(ish00,indna)
        al3=alpha(ial3,ish00,indna)
        coef3=gcoefficients(ial3,ish00,indna)

        bigM=al1+al2+al3

! The "center of mass" is
!   (al1*r1+al2*r2+al3*r3)/(al1+al2+al3) 
! This is location of the gaussian that you get when
! the three gaussians are multiplied together.
        xcm=x*sqrt(abs(1.-cost**2))*al3/bigM
        ycm=0.
        zcm=(al2*y+al3*(y/2.+x*cost))/bigM

! bigQ is a factor that you get in front of the product gaussian...
        bigQ=exp(-rmu*d12**2-lilm*al3*dcm3**2/bigM)
! gM = integral from -inf to +inf of exp^(-M*x**2)
        gM=sqrt(pi/bigM)

! c2 = factor that appears often. It's an extra factor that comes up
! if you integrate x**2 e^(-Mx**2) versus e^(-M*x**2).
        c2=1./2./bigM
! extra factor with x**4 e^(-Mx**2)... 
        c4=3./4./bigM**2

        factor=0.

! ordering for poly:
!    1=1
!    2=x
!    3=y
!    4=z
!    5=x**2
!    6=y**2
!    7=z**2
!    8=x*y
!    9=x*z
!   10=y*z

! What is poly? Basically coefficients of polynomial factors...
! for example:
! p_z state centered on atom 2 is basically:
!    (z-z2) = (z-zcm)+(zcm-z2)
! Then poly(1)=(zcm-z2), poly(4)=1., and all other polys zero.
! another example:
! d_xy on 1st center ---> (x-x1)*(y-y1) = 
!    (x-xcm)*(y-ycm)+(ycm-y1)*(x-xcm)+(xcm-x1)*(y-ycm)+(ycm-y1)*(xcm-x1)
! So poly(1)=(ycm-y1)*(xcm-x1), poly(2)=(ycm-y1), poly(3)=(xcm-x1), 
! poly(8)=1., and all other polys zero.

! JKT gauss. note. A lot of the coordinates below are zero:
!    y1,y2,ycm = x1,x2 = z1 = 0.
! I've left the general expressions in as comments.


! ===> get poly1

! S	
        if (l1.eq.0) then

           poly1(1)=1.

! P
        else if (l1.eq.1) then

           if (imu.eq.-1) then
!	      poly1(1)=(ycm-y1) = zero
              poly1(3)=1.
           else if (imu.eq.0) then
!	      poly1(1)=(zcm-z1)
              poly1(1)=zcm
              poly1(4)=1.
           else if (imu.eq.1) then
!	      poly1(1)=(xcm-x1)
              poly1(1)=xcm
              poly1(2)=1.
           end if

! D 
        else if (l1.eq.2) then

           if (imu.eq.-2) then
!	      poly1(1)=(xcm-x1)*(ycm-y1) = zero
!	      poly1(2)=(ycm-y1) = zero
!	      poly1(3)=(xcm-x1)
              poly1(3)=xcm
              poly1(8)=1.
           else if (imu.eq.-1) then
!	      poly1(1)=(ycm-y1)*(zcm-z1) = zero
!	      poly1(3)=(zcm-z1)
              poly1(3)=zcm
!	      poly1(4)=(ycm-y1) = zero
              poly1(10)=1.
           else if (imu.eq.0) then
!	      poly1(1)=2.*(zcm-z1)**2-(xcm-x1)**2-(ycm-y1)**2
              poly1(1)=2.*zcm**2-xcm**2
!	      poly1(2)=-2.*(xcm-x1)
              poly1(2)=-2.*xcm
!	      poly1(3)=-2.*(ycm-y1) = zero
!	      poly1(4)=4.*(zcm-z1) 
              poly1(4)=4.*zcm
              poly1(5)=-1.
              poly1(6)=-1.
              poly1(7)=2.
           else if (imu.eq.1) then
!	      poly1(1)=(xcm-x1)*(zcm-z1)
              poly1(1)=xcm*zcm
!	      poly1(2)=(zcm-z1)
              poly1(2)=zcm
!	      poly1(4)=(xcm-x1)
              poly1(4)=xcm
              poly1(9)=1.
           else if (imu.eq.2) then
!	      poly1(1)=(xcm-x1)**2-(ycm-y1)**2
              poly1(1)=xcm**2
!	      poly1(2)=(xcm-x1)
              poly1(2)=xcm
!	      poly1(3)=-(ycm-y1) = zero
              poly1(5)=1.
              poly1(6)=-1.
           end if

        end if

! ===> get poly2

! S
        if (l2.eq.0) then

           poly2(1)=1.

! P
        else if (l2.eq.1) then

           if (jmu.eq.-1) then
!              poly2(1)=(ycm-y2) = zero
              poly2(3)=1.
           else if (jmu.eq.0) then
              poly2(1)=(zcm-z2)
              poly2(4)=1.
           else if (jmu.eq.1) then
!              poly2(1)=(xcm-x2)
              poly2(1)=xcm
              poly2(2)=1.
           end if

! D
        else if (l2.eq.2) then

           if (jmu.eq.-2) then
!              poly2(1)=(xcm-x2)*(ycm-y2) = zero
!              poly2(2)=(ycm-y2) = zero
              poly2(3)=(xcm-x2)
              poly2(8)=1.
           else if (jmu.eq.-1) then
!              poly2(1)=(ycm-y2)*(zcm-z2) = zero
              poly2(3)=(zcm-z2)
!              poly2(4)=(ycm-y2) = zero
              poly2(10)=1.
           else if (jmu.eq.0) then
!              poly2(1)=2.*(zcm-z2)**2-(xcm-x2)**2-(ycm-y2)**2
              poly2(1)=2.*(zcm-z2)**2-xcm**2
!              poly2(2)=-2.*(xcm-x2)
              poly2(2)=-2.*xcm
!              poly2(3)=-2.*(ycm-y2) = zero
              poly2(4)=4.*(zcm-z2)
              poly2(5)=-1.
              poly2(6)=-1.
              poly2(7)=2.
           else if (jmu.eq.1) then
!              poly2(1)=(xcm-x2)*(zcm-z2)
              poly2(1)=xcm*(zcm-z2)
              poly2(2)=(zcm-z2)
!              poly2(4)=(xcm-x2)
              poly2(4)=xcm
              poly2(9)=1.
           else if (jmu.eq.2) then
!              poly2(1)=(xcm-x2)**2-(ycm-y2)**2
              poly2(1)=xcm**2
!              poly2(2)=(xcm-x2)
              poly2(2)=xcm
!              poly2(3)=-(ycm-y2) = zero
              poly2(5)=1.
              poly2(6)=-1.
           end if


        end if

! ===> get factor 

! This is what you get when you multiply the two polynomials together
! and integrate against the product gaussian...

! SS
        if (l1.eq.0 .and. l2.eq.0) then
!	   factor=poly1(1)*poly2(1)
           factor=1.
! SP
        else if (l1.eq.0 .and. l2.eq.1) then
!	   factor=poly1(1)*poly2(1)
           factor=poly2(1)
! PS
        else if (l1.eq.1 .and. l2.eq.0) then
!	   factor=poly1(1)*poly2(1)
           factor=poly1(1)
! PP
        else if (l1.eq.1 .and. l2.eq.1) then
           do ipoly=2,4
              factor=factor+poly1(ipoly)*poly2(ipoly)
           end do
           factor=c2*factor+poly1(1)*poly2(1)

! SD and DS
        else if (l1.eq.0 .and. l2.eq.2) then
           factor=poly2(1)+c2*(poly2(5)+poly2(6)+poly2(7))
        else if (l1.eq.2 .and. l2.eq.0) then
           factor=poly1(1)+c2*(poly1(5)+poly1(6)+poly1(7))

! PD and DP
        else if (l1.eq.1 .and. l2.eq.2) then
           factor=poly1(1)*poly2(1)
           do ipoly=2,4
              factor=factor+c2*poly1(ipoly)*poly2(ipoly)
           end do
           do ipoly=5,7
              factor=factor+c4*poly1(1)*poly2(ipoly)
           end do
        else if (l1.eq.2 .and. l2.eq.1) then
           factor=poly2(1)*poly1(1) 
          do ipoly=2,4
              factor=factor+c2*poly2(ipoly)*poly1(ipoly)
           end do
           do ipoly=5,7
              factor=factor+c4*poly2(1)*poly1(ipoly)
           end do

! DD
        else if (l1.eq.2 .and. l2.eq.2) then
           factor=poly1(1)*poly2(1)
           do ipoly=5,7
              factor=factor+c2*( poly1(1)*poly2(ipoly) + &
     &                           poly1(ipoly)*poly2(1) )
           end do
           do ipoly=2,4
              factor=factor+c2*poly1(ipoly)*poly2(ipoly)
           end do
           do ipoly=5,7
              factor=factor+c4*poly1(ipoly)*poly2(ipoly)
           end do
           do ipoly=8,10
              factor=factor+c2**2*poly1(ipoly)*poly2(ipoly)
           end do
        end if

!	term(ial1,ial2,ial3)=factor*coef1*coef2*coef3

        rmatel=rmatel+factor*coef1*coef2*coef3*bigQ*gM**3

        end do
        end do
        end do

! Put in spherical harmonic type factors...
        rmatel=rmatel*gfactor(l1,imu)*gfactor(l2,jmu)

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
        return
        end
