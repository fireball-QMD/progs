! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! University of Regensburg - Juergen Fritsch
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

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! getYlm.f90
! Program Description
! ===========================================================================
! getYlm returns spherical harmonics part of wavefunctions at given point
!
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: jelinekp@fzu.cz
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine getYlm (l, rvec, psi, dpsi)

   use dimensions
   use constants_fireball

   implicit none

! Argument Declaration and Description
! ===========================================================================

   integer, intent(in)                :: l        ! l-number
!   integer, intent(in)                :: m        ! m-number
   real, dimension (3), intent(in)    :: rvec     ! r(:)
   real, dimension (5), intent(out)   :: psi      ! psi(x)
   real, dimension (3,5), intent(out) :: dpsi     ! dpsi/dx

! Local Parameters and Data Declaration
! ===========================================================================
   real, parameter :: restol = 1.0d0-5

! Local Variable Declaration and Description
! ===========================================================================

   real r
   real r2
   real x
   real y
   real z
   real fact
   real facts
   real factp
   real factd


! Procedure
! ===========================================================================

! --------------------------------------------------------------------
!                -------- Spherical harmonics -------
! --------------------------------------------------------------------
! s-orbital (l=0)
!     Y0m = sqrt(1/4*Pi)
!     d(Y1m)/dr = 0
! --------------------------------------------------------------------
! p-orbital (l=1)
!     Y1m = sqrt(3/4*Pi)/r * Y1(r)
!            m=-1:   Y1(r) = y
!            m= 0:   Y1(r) = z
!            m=+1:   Y1(r) = x
!
!     d(Y1m)/dr = sqrt(2/4*Pi)d(Y1m(r)/r)/dr
!            m=-1:   Y2(r) =
!            m= 0:   Y2(r) =
!            m=+1:   Y2(r) =
!
! --------------------------------------------------------------------
! d-orbital (l=2)
!     Y2m = sqrt(15/4*Pi)/r**2 * Y2m(r)
!     d(Y2m)/dr = sqrt(15/4*Pi)/r**2 * d(Y2m(r))/dr
!            m=-2:   Y2(r) = x*y
!            m=-1:   Y2(r) = y*z
!            m= 0:   Y2(r) = (3z**2-r**2)/sqrt(12)
!            m=+1:   Y2(r) = x*z
!            m=+2:   Y2(r) = (x**2-y**2)/2
!
! ===========================================================================

! set constants_fireball
   facts = sqrt(1.0d0/(4.0d0*pi))
   factp = sqrt(3.0d0/(4.0d0*pi))
   factd = sqrt(15.0d0/(4.0d0*pi))

   psi = 0.0d0
   dpsi = 0.0d0
   x = 1.0d0*rvec(1)
   y = 1.0d0*rvec(2)
   z = 1.0d0*rvec(3)
   r = sqrt(x**2 + y**2 + z**2)

   if (r .lt. 0.0001d0 ) then
      if (l .ne. 0) then
         psi = 0.0d0
      else
         psi(1) = facts
      endif
      return
   endif

   select case ( l )
! -------------------------------------------------------------
! s-orbital
! -------------------------------------------------------------
    case (0)
      fact = facts
! value
      psi(1) = fact
! derivative
      dpsi (1,1) = 0.0d0
      dpsi (2,1) = 0.0d0
      dpsi (3,1) = 0.0d0
! -------------------------------------------------------------
! p-orbital
! -------------------------------------------------------------
    case (1)
! what's the value of spehrical harmonics if r=0 ??
      if (r .lt. restol ) then
         psi = 0.0d0
         return
      endif

      fact = factp/r
! value
      psi(1) = y*fact
      psi(2) = z*fact
      psi(3) = x*fact
! derivative
      dpsi (1,1) = 0.0d0
      dpsi (2,1) = 0.0d0
      dpsi (3,1) = 0.0d0
      dpsi (1,2) = 0.0d0
      dpsi (2,2) = 0.0d0
      dpsi (3,2) = 0.0d0
      dpsi (1,3) = 0.0d0
      dpsi (2,3) = 0.0d0
      dpsi (3,3) = 0.0d0
! -------------------------------------------------------------
! d-orbital
! -------------------------------------------------------------
    case (2)
! what's the value of spehrical harmonics if r=0 ??
      if (r .lt. restol ) then
         psi = 0.0d0
         return
      endif
      r2 = r**2
      fact = factd/r2
! value
      psi(1) = x*y * fact
      psi(2) = y*z * fact
      psi(3) = (3.0d0*z**2 - r2)/sqrt(12.0d0) * fact
      psi(4) = x*z * fact
      psi(5) = (x**2 - y**2)/2.0d0 * fact
! derivative
      dpsi (1,1) = 0.0d0
      dpsi (2,1) = 0.0d0
      dpsi (3,1) = 0.0d0
      dpsi (1,2) = 0.0d0
      dpsi (2,2) = 0.0d0
      dpsi (3,2) = 0.0d0
      dpsi (1,3) = 0.0d0
      dpsi (2,3) = 0.0d0
      dpsi (3,3) = 0.0d0
      dpsi (1,4) = 0.0d0
      dpsi (2,4) = 0.0d0
      dpsi (3,4) = 0.0d0
      dpsi (1,5) = 0.0d0
      dpsi (2,5) = 0.0d0
      dpsi (3,5) = 0.0d0

! -------------------------------------------------------------
! control block
    case default
       write (*,*) ' Error: Wrong l-number !!!'
       stop
   end select

! Format Statements
! ===========================================================================

   return
 end subroutine getYlm
