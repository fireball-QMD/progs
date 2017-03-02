! =====================================================================
! This program gets derivatives of three-center matrix elements.
! Very similar to gauss_matel but with derivatives 
! ---> See gauss_matel subroutine for more details.
! =====================================================================

        subroutine Dgelements(l1,l2,ish,jsh,im,jm,in1,in2,indna,  &
     &                          x,y,cost,rmatel,iparam,ish00)

        use dimensions
        use interactions
        use gauss
        use constants_fireball
        implicit none

        integer maxpoly

! for sp3 maxpoly is 4:  1,x,y,z
! for sp3d5 maxpoly is 10: 1,x,y,z,x**2,y**2,z**2,xy,xz,yz
! This isn't set up to do f-orbitals yet...
        parameter(maxpoly=10)

! Input
        integer, intent (in) :: in1, in2, indna

        integer, intent (in) :: l1, l2

        real, intent (in) :: cost
        real, intent (in) :: x
        real, intent (in) :: y

        integer ish,jsh,im,jm

        integer iparam 

! output
        real rmatel

! local

        integer ial1,ial2,ial3
! d12 is the distance between centers 1 and 2
! dcm3 is the distance from the "center of mass" of centers
! 1 and 2 to center 3.
        real rmu,lilm,bigM,d12,dcm3,bigQ,gM,xcm,ycm,zcm
        real dQ
        real al1,al2,al3
        real coef1,coef2,coef3
        real factor,dfactor
        real c2,c4

        real x1,y1,z1,x2,y2,z2,x3,y3,z3
! derivatives of xcm,ycm,zcm,x,y,and z...
        real dxcm,dycm,dzcm
        real dx1,dy1,dz1,dx2,dy2,dz2,dx3,dy3,dz3

        real, dimension (maxpoly) :: poly1,poly2,dpoly1,dpoly2
        integer ipoly

! WANG gauss-XC
        integer ish00


! =====================================================================
! initialize poly1 poly2 dpoly1 and dpoly2 to zero
        do ipoly=1,maxpoly
        poly1(ipoly)=0.
        poly2(ipoly)=0.
        dpoly1(ipoly)=0.
        dpoly2(ipoly)=0.
        end do

! Coordinates...

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

! Derivatives of coordinates...
        dx1=0.
        dy1=0.
        dz1=0.
        
        dx2=0.
        dy2=0.
        dz2=0.
        
        dx3=0.
        dy3=0.
        dz3=0.
        
        if (iparam.eq.1) then
           dz1=-1./2.
           dz2=1./2.
        else if (iparam.eq.2) then
           dx3=1.
        else if (iparam.eq.3) then
           dz3=1.
        end if

! ===============
! Start procedure
! ===============

        rmatel=0.

! alpha1
        do ial1=1,nalpha(ish,in1)
        al1=alpha(ial1,ish,in1)
        coef1=gcoefficients(ial1,ish,in1)

! alpha2
        do ial2=1,nalpha(jsh,in2)
        al2=alpha(ial2,jsh,in2)
        coef2=gcoefficients(ial2,jsh,in2)

        rmu=al1*al2/(al1+al2)
        lilm=al1+al2
        dcm3=(y/2.-al2*y/(al1+al2)+x*cost)**2+x**2*(1.-cost**2)
        dcm3=sqrt(dcm3)

! alpha3
! WANG gauss-XC
! ish00 = -1 for electron density
! ish00 =  0 for NA potential
        do ial3=1,nalpha(ish00,indna)
        al3=alpha(ial3,ish00,indna)
        coef3=gcoefficients(ial3,ish00,indna)

! OLD statements 
!        do ial3=1,nalpha(0,indna)
!        al3=alpha(ial3,0,indna)
!        coef3=gcoefficients(ial3,0,indna)

        bigM=al1+al2+al3

! coordinates...
        xcm=x*sqrt(1.-cost**2)*al3/bigM
        ycm=0.
        zcm=(al2*y+al3*(y/2.+x*cost))/bigM

! derivatives of coordinates...
        dxcm=0.
        dycm=0.
        dzcm=0.

        if (iparam.eq.1) then
           dzcm=(-al1+al2)/2./bigM 
        else if (iparam.eq.2) then
           dxcm=al3/bigM
        else if (iparam.eq.3) then
           dzcm=al3/bigM
        end if

! Big Q
        bigQ=exp(-rmu*d12**2-lilm*al3*dcm3**2/bigM)
        gM=sqrt(pi/bigM)

! Derivative of bigQ...
        if (iparam.eq.1) then
!           dQ=-2.*rmu*y+2.*al3*(al1-al2)*((al1*z1+al2*z2)/lilm-z3)/bigM
! WANG gauss-XC, it seems there is some error related to dQ
            dQ=-2.*rmu*y+al3*(al1-al2)*((al1*z1+al2*z2)/lilm-z3)/bigM
        else if (iparam.eq.2) then
            dQ=2.*lilm*al3*( (al1*x1+al2*x2)/lilm - x3 )/bigM
         else
            dQ=2.*lilm*al3*( (al1*z1+al2*z2)/lilm - z3 )/bigM
         end if 
         dQ=dQ*bigQ

! factor that appears often... It's an extra factor that comes up
! if you integrate x^2 e^(-Mx^2) versus e^(-Mx^2).
        c2=1./2./bigM
! extra factor with x**4 e^(-Mx**2)...
        c4=3./4./bigM**2

        factor=0.
        dfactor=0.

! The general expressions are left as comments...


! ===> get poly1 and dpoly1

! S
        if (l1.eq.0) then

           poly1(1)=1.

! P	
        else if (l1.eq.1) then

           if (im.eq.-1) then
!              poly1(1)=(ycm-y1) = zero
!	      dpoly1(1)=dycm-dy1 = zero
              poly1(3)=1.
           else if (im.eq.0) then
!              poly1(1)=(zcm-z1)
               poly1(1)=zcm
               dpoly1(1)=dzcm-dz1
               poly1(4)=1.
           else if (im.eq.1) then
!              poly1(1)=(xcm-x1)
               poly1(1)=xcm
!	      dpoly1(1)=dxcm-dx1
               dpoly1(1)=dxcm
              poly1(2)=1.
           end if

! D
        else if (l1.eq.2) then

           if (im.eq.-2) then
!              poly1(1)=(xcm-x1)*(ycm-y1) = zero
!	      dpoly1(1)=(dxcm-dx1)*(ycm-y1)+(xcm-x1)*(dycm-dy1) = zero
!              poly1(2)=(ycm-y1) = zero
!	      dpoly1(2)=dycm-dy1 = zero
!              poly1(3)=(xcm-x1)
               poly1(3)=xcm
!	      dpoly1(3)=dxcm-dx1
               dpoly1(3)=dxcm
               poly1(8)=1.
           else if (im.eq.-1) then
!              poly1(1)=(ycm-y1)*(zcm-z1) = zero
!	      dpoly1(1)=(dycm-dy1)*(zcm-z1)+(ycm-y1)*(dzcm-dz1) = zero
!              poly1(3)=(zcm-z1)
               poly1(3)=zcm
               dpoly1(3)=dzcm-dz1
!              poly1(4)=(ycm-y1) = zero
!	      dpoly1(4)=dycm-dy1 = zero
              poly1(10)=1.
           else if (im.eq.0) then
!              poly1(1)=2.*(zcm-z1)**2-(xcm-x1)**2-(ycm-y1)**2
               poly1(1)=2.*zcm**2-xcm**2
!	      dpoly1(1)=4.*(zcm-z1)*(dzcm-dz1) - &
!     &                  2.*(xcm-x1)*(dxcm-dx1) - &
!     &                  2.*(ycm-y1)*(dycm-dy1)
              dpoly1(1)=4.*zcm*(dzcm-dz1) - 2.*xcm*dxcm
!              poly1(2)=-2.*(xcm-x1)
               poly1(2)=-2.*xcm
!	      dpoly1(2)=-2.*(dxcm-dx1)
              dpoly1(2)=-2.*dxcm
!              poly1(3)=-2.*(ycm-y1) = zero
!	      dpoly1(3)=-2.*(dycm-dy1) = zero
!              poly1(4)=4.*(zcm-z1)
               poly1(4)=4.*zcm
               dpoly1(4)=4.*(dzcm-dz1)
               poly1(5)=-1.
               poly1(6)=-1.
              poly1(7)=2.
           else if (im.eq.1) then
!              poly1(1)=(xcm-x1)*(zcm-z1)
               poly1(1)=xcm*zcm
!	      dpoly1(1)=(dxcm-dx1)*(zcm-z1)+(xcm-x1)*(dzcm-dz1)
               dpoly1(1)=dxcm*zcm + xcm*(dzcm-dz1)
!              poly1(2)=(zcm-z1)
               poly1(2)=zcm
               dpoly1(2)=(dzcm-dz1)
!              poly1(4)=(xcm-x1)
               poly1(4)=xcm
!	      dpoly1(4)=(dxcm-dx1)
               dpoly1(4)=dxcm
              poly1(9)=1.
           else if (im.eq.2) then
!              poly1(1)=(xcm-x1)**2-(ycm-y1)**2
               poly1(1)=xcm**2
!	      dpoly1(1)=2.*(xcm-x1)*(dxcm-dx1) - &
!     &                  2.*(ycm-y1)*(dycm-dy1)
               dpoly1(1)=2.*xcm*dxcm
!              poly1(2)=(xcm-x1)
               poly1(2)=xcm
!	      dpoly1(2)=dxcm-dx1
              dpoly1(2)=dxcm
!              poly1(3)=-(ycm-y1) = zero
!	      dpoly1(3)=-(dycm-dy1) = zero
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

           if (jm.eq.-1) then
!              poly2(1)=(ycm-y2) = zero
!	      dpoly2(1)=dycm-dy2 = zero
              poly2(3)=1.
           else if (jm.eq.0) then
              poly2(1)=(zcm-z2)
              dpoly2(1)=dzcm-dz2
              poly2(4)=1.
           else if (jm.eq.1) then
!              poly2(1)=(xcm-x2)
              poly2(1)=xcm
!	      dpoly2(1)=dxcm-dx2
              dpoly2(1)=dxcm
              poly2(2)=1.
           end if

! D
        else if (l2.eq.2) then

           if (jm.eq.-2) then
!              poly2(1)=(xcm-x2)*(ycm-y2) = zero
!	      dpoly2(1)=(dxcm-dx2)*(ycm-y2)+(xcm-x2)*(dycm-dy2) = zero
!              poly2(2)=(ycm-y2) = zero
!	      dpoly2(2)=(dycm-dy2) = zero
!              poly2(3)=(xcm-x2)
               poly2(3)=xcm
!	      dpoly2(3)=(dxcm-dx2)
               dpoly2(3)=dxcm
              poly2(8)=1.
           else if (jm.eq.-1) then
!              poly2(1)=(ycm-y2)*(zcm-z2) = zero
!	      dpoly2(1)=(dycm-dy2)*(zcm-z2)+(ycm-y2)*(dzcm-dz2) = zero
              poly2(3)=(zcm-z2)
               dpoly2(3)=(dzcm-dz2)
!              poly2(4)=(ycm-y2) = zero
!	      dpoly2(4)=(dycm-dy2) = zero
              poly2(10)=1.
           else if (jm.eq.0) then
!              poly2(1)=2.*(zcm-z2)**2-(xcm-x2)**2-(ycm-y2)**2
               poly2(1)=2.*(zcm-z2)**2-xcm**2
!	      dpoly2(1)=4.*(zcm-z2)*(dzcm-dz2) - &
!     &                  2.*(xcm-x2)*(dxcm-dx2) - &
!     &                  2.*(ycm-y2)*(dycm-dy2)
               dpoly2(1)=4.*(zcm-z2)*(dzcm-dz2) - 2.*xcm*dxcm
!	      poly2(2)=-2.*(xcm-x2)
               poly2(2)=-2.*xcm
!	      dpoly2(2)=-2.*(dxcm-dx2)
               dpoly2(2)=-2.*dxcm
!              poly2(3)=-2.*(ycm-y2) = zero
!	      dpoly2(3)=-2.*(dycm-dy2) = zero
              poly2(4)=4.*(zcm-z2)
               dpoly2(4)=4.*(dzcm-dz2)
              poly2(5)=-1.
              poly2(6)=-1.
              poly2(7)=2.
           else if (jm.eq.1) then
!              poly2(1)=(xcm-x2)*(zcm-z2)
               poly2(1)=xcm*(zcm-z2)
!	      dpoly2(1)=(dxcm-dx2)*(zcm-z2)+(xcm-x2)*(dzcm-dz2)
               dpoly2(1)=dxcm*(zcm-z2)+xcm*(dzcm-dz2)
              poly2(2)=(zcm-z2)
               dpoly2(2)=(dzcm-dz2)
!              poly2(4)=(xcm-x2)
               poly2(4)=xcm
!	      dpoly2(4)=(dxcm-dx2)
               dpoly2(4)=dxcm
              poly2(9)=1.
           else if (jm.eq.2) then
!              poly2(1)=(xcm-x2)**2-(ycm-y2)**2
               poly2(1)=xcm**2
!	      dpoly2(1)=2.*(xcm-x2)*(dxcm-dx2) - &
!     &                  2.*(ycm-y2)*(dycm-dy2)
               dpoly2(1)=2.*xcm*dxcm
!              poly2(2)=(xcm-x2)
               poly2(2)=xcm
!	      dpoly2(2)=(dxcm-dx2)
               dpoly2(2)=dxcm
!              poly2(3)=-(ycm-y2) = zero
!	      dpoly2(3)=-(dycm-dy2) = zero
              poly2(5)=1.
              poly2(6)=-1.
           end if

         end if

! ===> get factor and dfactor

! SS
        if (l1.eq.0 .and. l2.eq.0) then
!          factor=poly1(1)*poly2(1)
           factor=1.
! SP
        else if (l1.eq.0 .and. l2.eq.1) then
!          factor=poly1(1)*poly2(1)
           factor=poly2(1)
           dfactor=dpoly2(1)
! PS
        else if (l1.eq.1 .and. l2.eq.0) then
!          factor=poly1(1)*poly2(1)
           factor=poly1(1)
           dfactor=dpoly1(1)
! PP
        else if (l1.eq.1 .and. l2.eq.1) then
           do ipoly=2,4
              factor=factor+poly1(ipoly)*poly2(ipoly)
              dfactor=dfactor+dpoly1(ipoly)*poly2(ipoly) + &
     &                        poly1(ipoly)*dpoly2(ipoly)
           end do
           factor=c2*factor+poly1(1)*poly2(1)
           dfactor=c2*dfactor+dpoly1(1)*poly2(1) + &
     &                        poly1(1)*dpoly2(1)

! SD and DS
        else if (l1.eq.0 .and. l2.eq.2) then
           factor=poly2(1)+c2*(poly2(5)+poly2(6)+poly2(7))
           dfactor=dpoly2(1)+c2*(dpoly2(5)+dpoly2(6)+dpoly2(7))
        else if (l1.eq.2 .and. l2.eq.0) then
           factor=poly1(1)+c2*(poly1(5)+poly1(6)+poly1(7))
           dfactor=dpoly1(1)+c2*(dpoly1(5)+dpoly1(6)+dpoly1(7))

! PD and DP
        else if (l1.eq.1 .and. l2.eq.2) then
           factor=poly1(1)*poly2(1)
           dfactor=dpoly1(1)*poly2(1)+poly1(1)*dpoly2(1)
           do ipoly=2,4
              factor=factor+c2*poly1(ipoly)*poly2(ipoly)
              dfactor=dfactor+c2*( dpoly1(ipoly)*poly2(ipoly) + &
     &                             poly1(ipoly)*dpoly2(ipoly) )
           end do
           do ipoly=5,7
              factor=factor+c4*poly1(1)*poly2(ipoly)
              dfactor=dfactor+c4*( dpoly1(1)*poly2(ipoly) + &
     &                             poly1(1)*dpoly2(ipoly) )
           end do
        else if (l1.eq.2 .and. l2.eq.1) then
           factor=poly2(1)*poly1(1)
           dfactor=dpoly2(1)*poly1(1)+poly2(1)*dpoly1(1)
           do ipoly=2,4
              factor=factor+c2*poly2(ipoly)*poly1(ipoly)
              dfactor=dfactor+c2*( dpoly2(ipoly)*poly1(ipoly) + &
     &                             poly2(ipoly)*dpoly1(ipoly) )
           end do
           do ipoly=5,7
              factor=factor+c4*poly2(1)*poly1(ipoly)
              dfactor=dfactor+c4*( dpoly2(1)*poly1(ipoly) + &
     &                             poly2(1)*dpoly1(ipoly) )
           end do

! DD
        else if (l1.eq.2 .and. l2.eq.2) then
           factor=poly1(1)*poly2(1)
           dfactor=dpoly1(1)*poly2(1)+poly1(1)*dpoly2(1)
           do ipoly=5,7
              factor=factor+c2*( poly1(1)*poly2(ipoly) + &
     &                           poly1(ipoly)*poly2(1) )
              dfactor=dfactor+c2*( dpoly1(1)*poly2(ipoly) + &
     &                             poly1(1)*dpoly2(ipoly) + &
     &                             dpoly1(ipoly)*poly2(1) + &
     &                             poly1(ipoly)*dpoly2(1) )
           end do
           do ipoly=2,4
              factor=factor+c2*poly1(ipoly)*poly2(ipoly)
              dfactor=dfactor+c2*( dpoly1(ipoly)*poly2(ipoly) + &
     &                             poly1(ipoly)*dpoly2(ipoly) )
           end do
           do ipoly=5,7
              factor=factor+c4*poly1(ipoly)*poly2(ipoly)
              dfactor=dfactor+c4*( dpoly1(ipoly)*poly2(ipoly) + &
     &                             poly1(ipoly)*dpoly2(ipoly) )
           end do
           do ipoly=8,10
              factor=factor+c2**2*poly1(ipoly)*poly2(ipoly)
              dfactor=dfactor+c2**2*( dpoly1(ipoly)*poly2(ipoly) + &
     &                                poly1(ipoly)*dpoly2(ipoly) )
           end do
        end if 

        rmatel=rmatel+coef1*coef2*coef3*gM**3* &
     &                (bigQ*dfactor + factor*dQ)
        end do
        end do
        end do

        rmatel=rmatel*gfactor(l1,im)*gfactor(l2,jm)

        return
        end
! =======================================================================
