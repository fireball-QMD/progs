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

! getpsi.f90
! Program Description
! ===========================================================================
! getpsi returns amplitude of wavefunction in given point r 
!
! ===========================================================================
! Code written by:
! ===========================================================================
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162  
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: jelinekp@fzu.cz
!
! Program Declaration
! ===========================================================================
 subroutine getpsi (ispec, issh, xin, psi, dpsi)

   use dimensions
   use constants_fireball
   use wavefunction

   implicit none

! Argument Declaration and Description
! ===========================================================================

   integer, intent(in)               :: ispec    ! ispec 
   integer, intent(in)               :: issh     ! shell
   real, intent(in)                  :: xin      ! r
   real, intent(out)                 :: psi      ! psi(x)
   real, intent(out)                 :: dpsi     ! dpsi/dx
 
! Local Parameters and Data Declaration
! ===========================================================================
! tol=tolerance (may be needed to avoid roundoff error in the calling program)
! if xin > xmax but within, say, .001% of xmax then ignore
   real, parameter :: tol = 3.0d-1
   real, parameter :: trash = 0.000001d0

! Local Variable Declaration and Description
! ===========================================================================
!
! -------------------------------------------------------------
   

   integer imid
   integer nnum
   integer ilo
   integer ihi
   
   real h
   real xmax
   real, parameter :: xmin = 0.0d0
   real xxp
   real xlo
   real xhi
   
! Superspline variables
   real a
   real b
   real wflo
   real wfhi
   real d2wflo
   real d2wfhi

 
! Procedure
! ===========================================================================

! --------------------------------------------------------------------
!                   -------- Radial part -------
! --------------------------------------------------------------------


! get number of mesh points
   nnum = mesh_wf(issh,ispec)
! get xmax of R(l)
   xmax = rmax_wf(issh,ispec)
! note : the points must be equally spaced and start at 0
   h = (xmax - xmin)/(nnum - 1)
! get distance
!   xin = sqrt(rvec(1)**2 + rvec(2)**2 + rvec(3)**2)
! check interpolation range
   if (xin .lt. xmin) then
      write (*,*) ' xin, xmin = ', xin, xmin
      write (*,*) '  error in getpsi : xin < xmin'
      stop ! Negative value is very bad - should never get there
   else if (xin .ge. xmax) then
      psi = wf(nnum,issh,ispec)
      dpsi = 0.0d0
      return
   else if (xin .lt. (xmin+h)) then
      psi = wf(1,issh,ispec)
      dpsi = 0.0d0
      return  
   else
      xxp = xin - xmin
! note : imid is the point to the right (code assumes this)
      imid = int((xxp+trash)/h) + 1
      xhi = real(imid-1)*h  
      xlo = real(imid-2)*h
      ihi = imid
      ilo = imid-1
   end if
!   write (202,*)  's xin:',xin
!   write (202,300) imid,xhi, xlo,ihi,ilo
! interpolation formulae
   a = (xhi - xxp)/h   
   b = (xxp - xlo)/h
   wfhi = wf(ihi,issh,ispec)
   wflo = wf(ilo,issh,ispec)
   d2wfhi = wf_spline(ihi,issh,ispec)
   d2wflo = wf_spline(ilo,issh,ispec)

! value
   psi = a*wflo + b*wfhi + ((a**3-a)*d2wflo + (b**3-b)*d2wfhi)*(h**2)/6.0d0
! 1-th derivative
   dpsi = (wfhi-wflo)/h - (3.0d0*a**2-1.0d0)*h*d2wflo                       &
 &        + (3.0d0*b**2-1.0d0)*h*d2wfhi
! 2-nd derivative
!   d2psi = a*d2wflo + b*d2wfhi

! Format Statements
! ===========================================================================
300 format (' s:',i6,2f14.7,2i6)
   return
 end subroutine getpsi
