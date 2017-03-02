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

! unocentros.f90
! Program Description
! ===========================================================================
!      This routine calculates the one-center exchange-correlation
! interaction.
!
! ===========================================================================
! Original Code written by Juergen Fritsch.
 
! Code rewritten by:
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
!
! Average density part written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162  
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
!
! ===========================================================================

!
! Program Declaration
! ===========================================================================
        subroutine unocentros (in1, iatom, iforce, itheory, itheory_xc,    &
     &                         exc_1c, muexc_1c, dccexc_1c, mu1xc)
     
        use dimensions
        use interactions
        use integrals
        use charges
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: iatom
        integer, intent(in) :: iforce
        integer, intent(in) :: in1
        integer, intent(in) :: itheory
        integer, intent(in) :: itheory_xc

! Output
        real, intent(out) :: exc_1c        ! XC energy term of DCC
        real, intent(out) :: muexc_1c      ! XC potential term of DCC
        real, intent(out) :: dccexc_1c     ! XC DCC term
        real, intent(out), dimension (numorb_max, numorb_max) :: mu1xc
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer imu
        integer ideriv
        integer iderivmin
        integer iderivmax
        integer issh
        integer l1, l2
        integer L1old, L2old
        integer m1, m2
        integer n1, n2
        integer inu
        integer jssh
        integer kssh
 
        real dqhelp
        real q_mu
        
        real, dimension (0:2) :: dqfact
        real, dimension (0:2) :: linfact
        real, dimension (0:2) :: quadfact
 
! jel -dq
        real, dimension (nsh_max) :: dqi
        integer iissh, jjssh
! end jel-dq
! Procedure
! ===========================================================================
! In this case we do not need to pass "isub", since we are going to perform
! a loop over it (isub = 0, 4 to get 1st order Taylor expansion).
!
! 0 = neutral, 1-2 = first derivatives
   if(itheory_xc .eq. 0 ) then  
! Initialize.
      exc_1c = 0.0d0
      dccexc_1c = 0.0d0
      muexc_1c = 0.0d0
      mu1xc = 0.0d0
 
      dqhelp = 0.0d0
      if (itheory .eq. 1) then
         if (dq(in1) .ne. 0.0d0) then
            do issh = 1, nssh(in1)
               dqhelp = dqhelp + (Qin(issh,iatom) - Qneutral(issh,in1))
            end do
 
! We have +dq and -dq. The range is 2*dq. Determine the fraction of this amount.
            dqhelp = dqhelp/dq(in1)
         end if
      end if
 
! ****************************************************************************
!
! We do a quadratic expansion:
!       f(q) = f(0)  +  f'(0)*q   +  1/2 f''(0)*q*q
!
! The derivatives are computed as:
!       f'(0)  = [ f(dq) - f(-dq) ] / 2dq
!       f''(0) = [ f(dq) - 2f(0) + f(-dq) ] / dq*dq
!
! We introduce linfac(0)   =  1.0
!              linfac(i)   =  -/+ * (1/2) *  q/dq        i=1,2
!              quadfac(0)  = -2.0 * (1/2) * (q/dq)**2
!              quadfac(i)  =  1.0 * (1/2) * (q/dq)**2    i=1,2
!
!       f(0) = f(dq=0) ; f(1) = f(-dq) ; f(2) = f(dq)
!
! With this, f(q) is given as:
!       f(q) = sum_i  (linfac(i) + quadfac(i))*f(i)
!
! ****************************************************************************
      linfact(0)  = 1.0d0
      linfact(1)  = -0.5d0*dqhelp
      linfact(2)  =  0.5d0*dqhelp
      quadfact(1) =  0.5d0*dqhelp**2
      quadfact(2) =  0.5d0*dqhelp**2
      quadfact(0) = -1.0d0*dqhelp**2
 
      dqfact(0:2) = linfact(0:2) + quadfact(0:2)
 
!       write (*,*) '  '
!       write (*,*) ' ************************************************ '
!       write (*,*) ' Whoops we set d n(exc-muxc) /dr=0 for charges. '
!       write (*,*) ' Juergen says that this term is worth skipping '
!       write (*,*) ' for now - I am not quite sure about that. '
!       write (*,*) ' ************************************************ '
 
! ideriv = 0, 1, 2 for neutral, -, and + dq.
! Caution: exc1c is Integral n*(exc-muxc) d3r.
      iderivmin = 0
      iderivmax = 0
      if (itheory .eq. 1) iderivmax = 2
      do ideriv = iderivmin, iderivmax
 
! All we are doing is y = mx + b. This is written in a complicated way.
         dccexc_1c = dccexc_1c + dqfact(ideriv)*exc1c_0(in1,ideriv)
 
! Here is the fixed method.
! The variable n1 lies in the middle of -L to + L.
         n1 = 0
         L1old = 0
         do issh = 1, nssh(in1)
            L1 = lssh(issh,in1)
 
! So how do we get in the middle for the L1 set.
!               n1 = n1 + L1old + L1 + 1
! Gets to end of L1old part-^     ^-- To get to the middle of the new L1 part.
            n1 = n1 + L1old + L1 + 1
            do M1 = - L1, L1
               imu = n1 + M1
               n2 = 0
               L2old = 0
               do jssh = 1, nssh(in1)
                  L2 = lssh(jssh,in1)
                  n2 = n2 + L2old + L2 + 1
                  do M2 = - L2, L2
                     inu = n2 + M2
                     if (M1 .eq. M2) then
                        mu1xc(inu,imu) =                                  &
     &         mu1xc(inu,imu) + dqfact(ideriv)*exc1c(in1,jssh,issh,ideriv)
                     end if
                  end do
                  L2old = L2
               end do

! We set force to aro. For a neutral atom, these atoms terms ARE ZERO.
! But the dq parts give position dependent terms which do have forces.
               if (iforce .eq. 1) then
!           set to zero above
                  if (itheory .eq. 1) then
! We set force to zero. For a neutral atom, these atoms terms ARE ZERO.
! But the dq parts give position dependent terms which do have forces.
! These forces are "supposedly" negligible. Ask JPL or JHF.
                  end if
               end if
            end do
            L1old = L1
         end do
        end do
       endif

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! A V E R A G E  D E N S I T Y :   
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Note:
! In teh case of a~full charge-transfer OLSXC method we have 
! 'exc' and '<imu|Vxc|inu>' and their 1. and 2. derivatives with respect to dq
! of each shell. 
!  
! Comment on Extended Hubbard, both Exc and Vxc have to be evaluated in Harris sense 
! (atomic charges)
! in addition we calculate Uexc_1c and Umuxc_1c (energy and potential)
! variables:
! exc_1c  ... double counting correction 
! Uexc_1c ... 
! McWEDA 
      if(itheory_xc .eq. 2) then 

! Initialize.
       exc_1c = 0.0d0
       dccexc_1c = 0.0d0
       muexc_1c = 0.0d0
       mu1xc = 0.0d0

! Evaluate charge transfer in regard of the neutral charge
! for Harris and EH dqi = 0; this makes terms correpsonding 
! charge transfer zero, that's what we want 
       dqi = 0.0d0
       if (itheory .eq. 1) then
         do issh = 1, nssh(in1)
            dqi(issh) = (Qin(issh,iatom) - Qneutral(issh,in1))
         end do
       end if

!dani.JOM [
! we want to do:  if(m2 .eq. m1 .and. l1 .eq. l2) mu1xc(imu,inu) = nuxc1c(in1,issh,jssh)
!but using get?ssh variables

       do imu = 1,num_orb(in1)
        m1   = getmssh(degelec(iatom)+imu)
        l1   = getlssh(degelec(iatom)+imu)
        issh = getissh(degelec(iatom)+imu)
        do inu = 1,num_orb(in1)
         m2   = getmssh(degelec(iatom)+inu) 
         l2   = getlssh(degelec(iatom)+inu)
         jssh = getissh(degelec(iatom)+inu)
         if( m1 .eq. m2 .and. l1 .eq. l2 ) then
          mu1xc(inu,imu) = nuxc1c(in1,jssh,issh)
          do kssh = 1,nssh(in1)
           mu1xc(inu,imu) = mu1xc(inu,imu) +        &
           &  dnuxc1c(in1,jssh,issh,kssh)*dqi(kssh)
          enddo
         endif
        end do
       end do


! Double counting correction (dccexc_1c) int <mu|exc - vxc|mu>*q
! Note: now we take vxc without charge transfer correction,
! so we don't include derivatives of vxc and we multiply it 
! by the neutral charge !!!   
      ! 'neutral' term

!dani.jel {
       do issh = 1,nssh(in1)
! energy term
        exc_1c = exc_1c + exc1c0(in1,issh,issh)*Qin(issh,iatom) 
! potential term
        muexc_1c = muexc_1c + nuxc1c(in1,issh,issh)*Qin(issh,iatom)
! DCC term
        dccexc_1c = dccexc_1c +                                               &
        & (exc1c0(in1,issh,issh) - nuxc1c(in1,issh,issh))*Qin(issh,iatom)
        do jssh = 1,nssh(in1)
         exc_1c = exc_1c + dexc1c(in1,issh,issh,jssh)*dqi(jssh)*Qin(issh,iatom) 
         muexc_1c = muexc_1c +                                                &
        &        dnuxc1c(in1,issh,issh,jssh)*dqi(jssh)*Qin(issh,iatom)
         dccexc_1c = dccexc_1c +                                              &
        &       ( dexc1c(in1,issh,issh,jssh) -                         &
                 dnuxc1c(in1,issh,issh,jssh) )*dqi(jssh)*Qin(issh,iatom)
        enddo
       enddo
      !dani.jel }
!--->double counting correction <mu|exc|mu>=<mu|exc0|mu>+<mu|exc0'*dqi|mu>     

      endif ! if( itheory_xc .eq. 2)
      

! Format Statements
! ===========================================================================
 
   return
 end subroutine unocentros
