! copyright info:
!
!                             @Copyright 2005
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek
! Brigham Young University - Hao Wang

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang 
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch

!
! fireball-qmd is a free (GPLv3) open project.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


! initNH.f90
! Program Description
! ===========================================================================
!       This initializes the nose-hoover chain.
!
! ===========================================================================
! Code written by:
! J. Keith
! ===========================================================================
!
! Program Declaration
! ===========================================================================
subroutine initNH(natoms,T_want)
use noseHoover
use constants_fireball
implicit none

! Passed variables
! ===========================================================================
integer, intent(in) :: natoms
real, intent(in) :: T_want
 
! Variable Declaration and Description
! ===========================================================================

real :: tau, gauss_num, sv!, Gauss_vel
integer i
logical noseHooverChainQ

! Procedure
! ===========================================================================
inquire (file = 'NH.optional', exist = noseHooverChainQ)
if (noseHooverChainQ) then
 ! read in NH.optional
 open (unit = 120, file = 'NH.optional', status = 'old')
 read(120,*)nnos
 print *,'num of chains',nnos
 allocate(omega(nnos))
 do i=1,nnos
  read(120,*)omega(i)
 end do
 close(120)
else ! these are the defaults
 nnos=4
 allocate(omega(nnos))
 tau=500 ! default 500 fs relaxation time
 omega=2*pi/tau
end if

allocate(xi(nnos))
allocate(v_xi(nnos))
allocate(Q_i(nnos))
allocate(G_i(nnos))
xi=0.0
kT=kb*T_want
gkT=natoms*3*kT
print *,'kT',kT
Q_i(1) = gkT/(omega(1)*omega(1)) !Q is eV*fs^2=M*L^2
write(*,*)'Q_i(1)',Q_i(1)
do i=2,nnos
 Q_i(i) = kT/(omega(i)*omega(i))
 write(*,*)'Q_i(i)',Q_i(i)
end do

! initialize the thermostat velocities with a Gaussian distribution
do i = 1, nnos
 sv = sqrt( kT / Q_i(i) )
 call gaussian_distro(gauss_num)
 v_xi(i) = sv * gauss_num
end do

G_i=0.0

debug=.false.

contains
subroutine gaussian_distro(num1)
 implicit none
 real, intent(out) :: num1
 real, save :: num2
 logical, save :: gaus_stored=.false.
 real x1, x2, w, rnum1, rnum2

 call random_seed 
 if(gaus_stored) then
  num1=num2
  gaus_stored=.false.
 else
  do
   call random_number(rnum1)
   call random_number(rnum2)
   x1 = 2.0 * rnum1 - 1.0;
   x2 = 2.0 * rnum2 - 1.0;
   w = x1 * x1 + x2 * x2;
   if (w<1.0) exit
  end do
  w = sqrt( (-2.0 * log( w ) ) / w );
  num1 = x1 * w;
  num2 = x2 * w;
  gaus_stored=.true.
 endif
end subroutine


!!$function ranf()
!!$!*
!!$implicit integer(a-z)
!!$
!!$real :: ranf
!!$!*
!!$!* Period parameters
!!$parameter(N     =  624)
!!$parameter(N1    =  N+1)
!!$parameter(M     =  397)
!!$parameter(MATA  = -1727483681)
!!$!*                                    constant vector a
!!$parameter(UMASK = -2147483648)
!!$!*                                    most significant w-r bits
!!$parameter(LMASK =  2147483647)
!!$!*                                    least significant r bits
!!$!* Tempering parameters
!!$parameter(TMASKB= -1658038656)
!!$parameter(TMASKC= -272236544)
!!$!*
!!$dimension mt(0:N-1)
!!$!*                     the array for the state vector
!!$common /block/mti,mt
!!$save   /block/
!!$data   mti/N1/
!!$!*                     mti==N+1 means mt[N] is not initialized
!!$!*
!!$dimension mag01(0:1)
!!$data mag01/0, MATA/
!!$save mag01
!!$!*                        mag01(x) = x * MATA for x=0,1
!!$!*
!!$   TSHFTU(y)=ishft(y,-11)
!!$   TSHFTS(y)=ishft(y,7)
!!$   TSHFTT(y)=ishft(y,15)
!!$   TSHFTL(y)=ishft(y,-18)
!!$!*
!!$   if(mti.ge.N) then
!!$
!!$   do kk = 0 , N-M-1
!!$
!!$     y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
!!$     mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
!!$
!!$   end do
!!$
!!$   do kk = N-M , N-2
!!$
!!$     y=ior(iand(mt(kk),UMASK),iand(mt(kk+1),LMASK))
!!$     mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
!!$
!!$   end do
!!$
!!$   y=ior(iand(mt(N-1),UMASK),iand(mt(0),LMASK))
!!$   mt(N-1)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
!!$   mti = 0
!!$
!!$   endif
!!$!*
!!$   y=mt(mti)
!!$   mti=mti+1
!!$   y=ieor(y,TSHFTU(y))
!!$   y=ieor(y,iand(TSHFTS(y),TMASKB))
!!$   y=ieor(y,iand(TSHFTT(y),TMASKC))
!!$   y=ieor(y,TSHFTL(y))
!!$!*
!!$   if(y.lt.0) then
!!$
!!$     ranf=(dble(y)+2.0d0**32)/(2.0d0**32-1.0d0)
!!$
!!$   else
!!$
!!$     ranf=dble(y)/(2.0d0**32-1.0d0)
!!$
!!$   endif
!!$!*
!!$end function
!!$
!!$function Gauss_vel()
!!$! generate normal distribution--Num Rec?
!!$
!!$integer :: Iset
!!$real fac, Gset, Rsq, v1, v2, ranf
!!$real Gauss_vel
!!$save Iset, Gset
!!$!external ranf
!!$
!!$data Iset /0/
!!$
!!$if (Iset==0) then
!!$
!!$1 call random_seed()
!!$ call random_number(ranf)
!!$ v1=2.0*ranf-1.0
!!$ call random_seed()
!!$ call random_number(ranf)
!!$ v2=2.0*ranf-1.0
!!$ Rsq=v1*v1+v2*v2
!!$ if( (Rsq >= 1.0) .or. (Rsq ==0.0)) goto 1
!!$
!!$!1 v1=2.0*ranf()-1.0
!!$!  v2=2.0*ranf()-1.0
!!$!  Rsq=v1*v1+v2*v2
!!$!  if( (Rsq >= 1.0) .or. (Rsq ==0.0)) goto 1
!!$
!!$ fac = sqrt(-2.0*log(Rsq)/Rsq)
!!$ Gset = v1*fac
!!$ Gauss_vel = v2*fac
!!$ Iset = 1
!!$else
!!$ Gauss_vel=Gset
!!$ Iset=0
!!$end if
!!$end function

end subroutine
