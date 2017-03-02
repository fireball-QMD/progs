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

! smoother.f90
! Program Description
! ===========================================================================
!       This is the smoother function.  The smoothing function is 
! smoother(r,rbegin,rend).  We define our answer as 
! smoother(r)*exact + (1 - smoother(r))*longrange.
! We also define x=frac=(r - rbegin)/(rend - rbegin) is some equations
! The following are the requirements for a sane smoother function:
! smoother(rend)=0
! smoother(rbegin)=1
! dsmoother <= 0    (smoother never increases)
! These are the requirements for a sane dsmoother:
! dsmoother(rend)=0
! dsmoother(rbegin)=0
! The following are nice to have:
! You might like symmetry: stn(x) = 1-stn(1-x)
! The maximum value of abs(dsmoother) to be small, reduce noise in derivatives
! Quick to evaluate on the computer: not (cos(pi*x)+1)/2.
!
! There are different functional forms:
! (1-x^n)^m is the "old" method
! Which obeys the requirements for n,m>1.  It is the smoothest
! for n,m=2 (smallest possible values).  Thus, a generic forth-order polynomial
! might be a better choice, since n,m=2 is 1 - 2x^2 + x^4 (no x or x^3 terms).
! A forth-order polynomial is the "new" method, but the above conditions
! force it to be of the form:
! 1 + 0x + (n-3)x^2 + (2-2n)x^3 + nx^4
! Where n must be within [-3,3], inclusive.
! For n=1, you recover the "old" method.
! To minizize abs(dsmoother), n is set to zero.
! To make symmetric, n is set to zero
! To get second-derivatives to be zero at x=0, n is to 3
! To get second-derivatives to be zero at x=1, n is to -3
!
! Finally, if you want to get second-derivatives to be zero at x=0,1 because
! you are doing a hessian, then you need a higher-order function.
! 1 + (-10-n)x^3 + (15+3n)x^4 + (-6-3n)x^5 + nx^6
! For symmetry, and mininizing abs(dsmoother), n=0 is probably best.
! This is not implemented.
!
!
! Now a word about choosing xi
! E = Sf + L(1-f),  where S=short, L=long, f=smoothing function.
! Note that this is well-behaved (xi>0), since E is bound by S and L.
! But what about derivatives (ie. forces):
! E' = S'f + L'(1-f) + Extra term.
! Which is well-behaved, except for the extra term which makes it:
! E' = S'f + L'(1-f) + f'(S-L)
! Keeping f' small is good, which is achieved by making xi smaller.
! Keeping (S-L) small is good also, which is achieved by making xi bigger.
! Clearly some sort of balance must be achieved, when choosing xi.
! 
! ===========================================================================
! Code written by:
! Kurt R. Glaesemann
! Henry Eyring Center for Theoretical Chemistry
! Department of Chemistry
! University of Utah
! 315 S. 1400 E.
! Salt Lake City, UT 84112-0850
! FAX 801-581-4353
! Office telephone 801-585-1078
! ===========================================================================
 
! Program Declaration
! ===========================================================================
        subroutine smoother (r, rend, xi, stn, dstn)
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
        real, intent (in) :: r
        real, intent (in) :: rend
        real, intent (in) :: xi
        real, intent (out) :: stn
        real, intent (out) :: dstn
 
! Local Parameters and Data Declaration
! ===========================================================================
        integer, parameter :: npower = 2
        integer, parameter :: mpower = 2
        integer, parameter :: scaler = 0

        logical, parameter :: old_method = .true.
 
! Local Variable Declaration and Description
! ===========================================================================
        real frac
        real rbegin
        real dum
 
! Procedure
! ===========================================================================
        rbegin = xi*rend
        if (r .lt. 0.0d0) then
         write (*,*) ' r < 0 in smoother *** error! '
         stop
        else if (r .gt. rend) then
         stn = 0.0d0
         dstn = 0.0d0
        else if (r .lt. rbegin) then
         stn = 1.0d0
         dstn = 0.0d0
        else
         frac = (r - rbegin)/(rend - rbegin)

         if (old_method) then
           dum = 1.0d0 - frac**npower
           stn = dum**mpower
! We also calculate d (smoother) /dr.
           dstn = - (mpower*npower)*(dum**(mpower - 1))                       &
     &                             *(frac**(npower - 1))/(rend - rbegin)
         else  ! new method
           stn = 1 + (scaler-3)*frac**2 + (2-2*scaler)*frac**3 + scaler*frac**4
           dstn=   2*(scaler-3)*frac  + 3*(2-2*scaler)*frac**2 + scaler*frac**3
         end if
        end if
 
! Format Statements
! ===========================================================================
 
        return
        end
