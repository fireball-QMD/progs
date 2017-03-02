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

! mixer.f90
! Program Description
! ===========================================================================

! ===========================================================================
! Code written by:

! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine mixer (natoms, itheory, ifixcharge, iwrtcharges)

   use charges
   use scf
   use interactions
   implicit none

! Argument Declaration and Description
! ===========================================================================
! input

! input and output
   integer, intent (in) :: natoms       ! number of atoms
   integer, intent (in) :: itheory      ! itheory
   integer, intent (in) :: ifixcharge   ! fixed charges
   integer, intent (in) :: iwrtcharges  ! fixed charges
   
! output
!   real, intent(out) :: sigma !

! Local Parameters and Data Declaration
! ===========================================================================


! Local Variable Declaration and Description
! ===========================================================================

   integer iatom
   integer in1
   integer issh
   integer imix

   real dqrms
   real dqmax
   real renorm
   real zcheck
   real zouttot
 
! Procedure
! ===========================================================================

! Check to see if input charges and output charges are within tolerance.
! If they are within tolerance, then perform a motional time step. If
! they are not within tolerence, then perform another iteration to make
! charges self-consistent.

! Also, only do this step if ifixcharge not equal to 1.
   if (ifixcharge .ne. 1) then

     dqrms = 0
     dqmax = -99
     do iatom = 1, natoms
       in1 = imass(iatom)
       do issh = 1, nssh(in1)
         dqmax = amax1(abs(Qin(issh,iatom) - Qout(issh,iatom)),dqmax)
         dqrms = dqrms + (Qin(issh,iatom) - Qout(issh,iatom))**2
       end do
     end do
     dqrms = sqrt(dqrms)/(2*natoms)

     write (*,300) dqrms
     write (*,301) dqmax

! ===========================================================================
!                                   mixing
! ===========================================================================
! We only do this part if we are not doing Harris, i.e. itheory = 0.
     if (itheory .eq. 1 .or. itheory .eq. 2) then
       imix = 0
       do iatom = 1, natoms
         in1 = imass(iatom)
         do issh = 1, nssh(in1)
           imix = imix + 1
           Qinmixer(imix) = Qin(issh,iatom)
           Qoutmixer(imix) = Qout(issh,iatom)
           if (ialgmix .eq. 4) then 
            mwe(imix) = 1.0d0
            drwe(imix) = 1.0d0
           endif
         end do
       end do 
    
! Mix Qinmixer and Qoutmixer to get a NEW Qinmixer, which is the NEW input for
! the next run. So Qinmixer is the thing this baby returns. We need a mixing
! factor, bmix, which is in controls.input.
! bmix = 0.08 means:  Careful mixing (only 8 percent of the new charge).
! idmix: Defines the mixing procedure: idmix=1 means simple mixing
! For larger idmix values, the choice of bmix becomes less important.
       if (iwrtcharges .eq. 1) then
         write (*,*) '  '
         write (*,*) ' Calling charge mixer '
       end if
! call mixing procedure
! modified by honza

! call mixing procedure
       select case (ialgmix)
       case (1)
          write(*,501) 'mixing with anderson algorithm'
          call anderson (Qoutmixer, Qinmixer, bmix, sigma, Kscf, idmix,    &
               &               imix , max_scf_iterations)
       case (2)
          write(*,501) 'mixing with broyden algorithm'
          call broyden (Qoutmixer, Qinmixer, bmix, sigma, Kscf, idmix,    &
               &               imix , max_scf_iterations)
       case (3)
          write(*,501) 'mixing with louie algorithm'
          call louie (Qoutmixer, Qinmixer, bmix, sigma, Kscf, idmix,    &
               &               imix , max_scf_iterations)
       case (4)
          write(*,501) 'mixing with pulay algorithm'
          call pulay (Qoutmixer, Qinmixer, bmix, sigma, Kscf, idmix,    &
               &               imix , max_scf_iterations)
       end select !ialgmix


! end of honza
       
       if (Kscf .gt. 1) then
         if (sigma .lt. sigmaold) then
           sigmaold = sigma
         end if
       else
         sigmaold = sigma
       end if

! Check convergence of charge; sigmatol is in scf.optional
       write (*,*) '  '
       if (sigma .lt. sigmatol) scf_achieved = .true.

       if (.not. scf_achieved) then
         imix = 0
         do iatom = 1, natoms
           in1 = imass(iatom)
           do issh = 1, nssh(in1)
             imix = imix + 1
             Qin(issh,iatom) = Qinmixer(imix)
           end do
         end do
       end if

     end if ! end if ( itheory 1,2)
! Skip to here if charges are fixed
    end if ! end if (ifixcharge)

! ===========================================================================
!                                 end mixing
! ===========================================================================
    zouttot = 0
    do iatom = 1, natoms
      in1 = imass(iatom)
      do issh = 1, nssh(in1)
        zouttot = zouttot + Qin(issh,iatom)
      end do
    end do
    renorm = (zouttot - ztot)/nssh_tot
    write (*,302) renorm
    zcheck = 0.0d0
    do iatom = 1, natoms
      in1 = imass(iatom)
      do issh = 1, nssh(in1)
        Qin(issh,iatom) = Qin(issh,iatom) - renorm
        zcheck = zcheck + Qin(issh,iatom)
      end do
    end do

! write out resume 
    write (*,*) ' (Before renormalization) zouttot = ', zouttot
    write (*,*) ' (After  renormalization)  zcheck = ', zcheck
    write (*,*) ' (What it must be)           ztot = ', ztot
    write (*,*) '  '
    write (*,303) sigma, sigmatol, Kscf


! Format Statements
! ===========================================================================
300     format (2x, ' Deviation (rms) of input/output charges = ', f9.4)
301     format (2x, ' Deviation (max) of input/output charges = ', f9.4)
302     format (2x, ' Renormalization of Qin:   renorm = ', 4x, f14.8)
303     format (' =============> sigma = ', e14.7,                           &
     &          ' Must be less than', e14.7, ' SCF step = ', i3)
501     format (2x,'Your chosen mixing routine is: ',a)
    
      return
    end subroutine mixer
