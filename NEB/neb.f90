! copyright info:
!
!                             @Copyright 2005
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek
!
! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola - Jun Wang
! Ohio State University - Dave Drabold

!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! neb.f90
! Program Description
! ===========================================================================
!  we use the improved  nigged elastic band method for finding MEP
! (minimum energy path) using new definition of tanget
! for detail see Journal od Chem. Phys. vol 113, no. 22, p. 9978 (2000)
!
! ===========================================================================
! Code written by:
! K. Kosmider
! and
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
 subroutine move_neb (itime, etot, force)

   use configuration
   use neb

   implicit none

! Argument Declaration and Description
! ===========================================================================
   integer, intent(in)                              :: itime
   real, intent(in)                                 :: etot

   real, dimension(3,natoms), intent(inout)         :: force

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================

   integer iatom
   integer img
   integer iimg_neb

   real fms
   real fms_max
   real rms
   real rms_max
   real ems_max

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Auxilliary parameters

   call get_iimg_neb (itime, iimg_neb)
   write (*,*) 'iimg_neb =',iimg_neb,'itime =',itime

   if (itime .gt. 2) then
! save energy & force of the given image
    detot_neb(iimg_neb) = etot_neb(iimg_neb) - etot
    etot_neb(iimg_neb) = etot
    do iatom = 1, natoms
     ftot_neb(:,iatom,iimg_neb) = force(:,iatom)
    enddo ! do iatom

    write (*,*)  ' ftot_neb :',iimg_neb
    do iatom = 1, natoms
     write (*,300) iatom, ftot_neb(:,iatom,iimg_neb)
    enddo
    write (*,*)


! next image
    if (iimg_neb .eq. (nimg_neb-1)) then
! here we move all images along projected NEB forces

     rms_max = 0.0d0
     fms_max = 0.0d0
     ems_max = 0.0d0

! loop over free images
     do img = 2,nimg_neb-1


      write (*,*) ' ------------ neb image no.',img
      write (*,*)

! 1. calc tangent
      call get_tang (natoms, img)
! 2. project out the spring forces
      call get_Fspring (natoms, img)
! 3. project out the rectangular part of total forces
      call get_Frec (natoms, img)
! 4. evaluate NEB forces
      call get_Fneb (natoms, img, fms)
      fms_max = max(fms_max,fms)
     enddo

! loop over free images
     do img = 2,nimg_neb-1
! 5. move the image according to the NEB forces
      call move_image (img, rms)
      rms_max = max(rms_max,rms)
      ems_max = max(ems_max,abs(detot_neb(img)))
     enddo ! do img

! write out convergence parameters
     write (*,99)   iter_neb
     write (*,100)  rms_max,tol_displ_neb
     write (*,101)  ems_max,tol_etot_neb
     write (*,102)  fms_max,tol_ftot_neb

! write out the actual snapshot of images
     call writeout_image (natoms, symbol, shifter)

! increase NEB iteration
     iter_neb = iter_neb + 1

! check convergence & max_number of iterations
! A. ftot
     if ((fms_max .lt. tol_ftot_neb) .and. (ems_max .lt. tol_etot_neb)) then
      write (*,*) ' Minimal energy pathway has been achievied within '
      write (*,*) ' the convergence criteria. Bye! '
      write (*,101)  ems_max,tol_etot_neb
      write (*,102)  fms_max,tol_ftot_neb
      stop
     endif

! check if the maximal number of NEB iteration is exceeded
     if (iter_neb .gt. niter_neb_max) then

      write (*,*) '  Maximal number of neb iteration =', niter_neb_max
      write (*,*) '  has been reached without convergence.'
      write (*,*) '  The program is going to stop here ....'
      stop

     endif

    endif ! if (iimg_neb)

!
    call get_iimg_neb (itime+1, iimg_neb)

! copy new image
    do iatom = 1, natoms
     ratom(:,iatom) = ratom_neb(:,iatom,iimg_neb)
    enddo ! do iatom

! Initial image
   else if (itime .eq. 1) then

! save energy of the initial image
    etot_neb(1) = etot
! copy the final image for next step (the final state)
    do iatom = 1, natoms
     ratom(:,iatom) = ratom_neb(:,iatom,nimg_neb)
    enddo ! do iatom
! Final image (itime =2)
   else
! save energy of the final image
    etot_neb(nimg_neb) =  etot
! copy 2nd image for next step
    do iatom = 1, natoms
     ratom(:,iatom) = ratom_neb(:,iatom,2)
    enddo ! do iatom
   endif ! if (itime .gt. 2)


   return

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
99  format (2x,' NEB_RES no. iter =', i4 )
100 format (2x,' NEB_RES displ :  RES =', f16.8,'  TOL = ',f16.8)
101 format (2x,' NEB_RES etot  :  RES =', f16.8,'  TOL = ',f16.8)
102 format (2x,' NEB_RES ftot  :  RES =', f16.8,'  TOL = ',f16.8)
300 format(i3,3f8.4)

   return

 end subroutine move_neb




!############################################
! subroutine returns actual image index
!############################################
 subroutine get_iimg_neb(itime, iimg_neb)

   use neb
   implicit none

! Argument Declaration and Description
! ===========================================================================
   integer,        intent(in)      ::      itime
   integer,        intent(out)     ::      iimg_neb

! Procedure
! ===========================================================================
   if(itime .gt. nimg_neb) then
    iimg_neb = mod(itime-nimg_neb-1,nimg_neb-2) + 2
   else if (itime .gt. 0 .and. itime .le. nimg_neb ) then
    if(itime .gt. 2) then
     iimg_neb = itime - 1
    else if(itime .eq. 2) then
     iimg_neb = nimg_neb
    else if(itime .eq. 1) then
     iimg_neb = itime
    else
     write(*,*) " Error : wrong itime value"
    endif
   else
    write(*,*) " Error : wrong itime value"
   endif

! Format Statements
! ===========================================================================
   return

 end subroutine

!############################################
! subroutine calcs the tangent of the image
!############################################
 subroutine get_tang (natoms, image)

   use neb

   implicit none

! Argument Declaration and Description
! ===========================================================================
   integer, intent(in)                              :: natoms
   integer, intent(in)                              :: image


! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer ix
   real  ee1
   real  ee2
   real  ee3
   real  demax
   real  demin
   real  norm
   real  r32
   real  r21

   real, dimension (3,natoms)  :: r1
   real, dimension (3,natoms)  :: r2
   real, dimension (3,natoms)  :: r3

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Auxilliary parameters

! image energies
! i - 1
   ee1 = etot_neb(image-1)
   do iatom = 1, natoms
    r1(:,iatom) = ratom_neb(:,iatom,image-1)
   enddo

! i
   ee2 = etot_neb(image)
   do iatom = 1, natoms
    r2(:,iatom) = ratom_neb(:,iatom,image)
   enddo

! i + 1
   ee3 = etot_neb(image+1)
   do iatom = 1, natoms
    r3(:,iatom) = ratom_neb(:,iatom,image+1)
   enddo

! case i+1 < i < i-1  tau- eqn.(9) JCP 113,9978(2000)
   if ((ee1 .gt. ee2) .and. (ee2 .gt. ee3)) then

    do iatom = 1, natoms
     do ix = 1,3
      tang(ix,iatom) = r2(ix,iatom) - r1(ix,iatom)
     enddo
    enddo

! case i+1 > i > i-1 tau+
   else if ((ee1 .lt. ee2) .and. (ee2 .lt. ee3)) then

    do iatom = 1, natoms
     do ix = 1,3
      tang(ix,iatom) = r3(ix,iatom) - r2(ix,iatom)
     enddo
    enddo
! case i+1 > i < i-1  or  i+1 < i > i-1
   else

    demax = max( abs(ee3-ee2),abs(ee1-ee2))
    demin = min( abs(ee3-ee2),abs(ee1-ee2))

! case i+1 > i-1
    if (ee3 .gt. ee1) then
! loop over atoms
     do iatom = 1, natoms
      do ix = 1,3
! eqn (10)
       tang(ix,iatom) = (r3(ix,iatom) - r2(ix,iatom))*demax      &
                      + (r2(ix,iatom) - r1(ix,iatom))*demin
      enddo
     enddo
! case i+1 < i-1
    else
! loop over atoms
     do iatom = 1, natoms
      do ix = 1,3
! eqn (10)
       tang(ix,iatom) = (r3(ix,iatom) - r2(ix,iatom))*demin      &
                      + (r2(ix,iatom) - r1(ix,iatom))*demax
      enddo
     enddo ! do iatom
    endif ! if (ee3 .gt. ee1)
   endif ! if ((ee1 .gt. ee2) .and. (ee2 .gt. ee3))

! simple bisection version of tangent
!  r32 = 0.0d0
!  r21 = 0.0d0
!  do iatom = 1,natoms
!   do ix = 1,3
!    r32 = r32 + (r3(ix,iatom) - r2(ix,iatom))**2
!    r21 = r21 + (r2(ix,iatom) - r1(ix,iatom))**2
!   enddo
!  enddo
!  r32 = sqrt(r32)
!  r21 = sqrt(r21)
!  do iatom = 1,natoms
!   do ix =1,3
!    tang(ix,iatom) = (r2(ix,iatom) - r1(ix,iatom))/r21 + (r3(ix,iatom) - r2(ix,iatom))/r32
!    tang(ix,iatom) = (r2(ix,iatom) - r1(ix,iatom))/r21
!   enddo
!  enddo

! normalize tangent
   norm = 0.0d0
   do iatom = 1,natoms
    do ix = 1,3
      norm = norm + tang(ix,iatom)**2
    enddo
   enddo

   norm = sqrt(norm)

   do iatom = 1,natoms
    do ix = 1,3
      tang(ix,iatom) = tang(ix,iatom)/norm
    enddo
   enddo


! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
100 format(i3,3f8.4)
   return

 end subroutine get_tang

!############################################
! subroutine calcs spring force
!############################################
 subroutine get_Fspring (natoms, image)

   use neb

   implicit none

! Argument Declaration and Description
! ===========================================================================
   integer, intent(in)                              :: natoms
   integer, intent(in)                              :: image

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer ix

   real  r32
   real  r21
   real  dot
   real  kr12r23

   real, dimension (3,natoms)  :: r1
   real, dimension (3,natoms)  :: r2
   real, dimension (3,natoms)  :: r3

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Auxilliary parameters

! i - 1
!   write (*,*)  ' RATOM :',image-1
   do iatom = 1, natoms
    r1(:,iatom) = ratom_neb(:,iatom,image-1)
!    write (*,200) iatom, r1(:,iatom)
   enddo
! i
!   write (*,*)  ' RATOM :',image
   do iatom = 1, natoms
    r2(:,iatom) = ratom_neb(:,iatom,image)
!    write (*,200) iatom, r2(:,iatom)
   enddo
! i + 1
!   write (*,*)  ' RATOM :',image+1
   do iatom = 1, natoms
    r3(:,iatom) = ratom_neb(:,iatom,image+1)
!    write (*,200) iatom, r3(:,iatom)
   enddo

! calc norm
   r32 = 0.0d0
   r21 = 0.0d0
   do iatom = 1, natoms
    do ix = 1,3
     r32 = r32 + (r3(ix,iatom) - r2(ix,iatom))**2
     r21 = r21 + (r2(ix,iatom) - r1(ix,iatom))**2
    enddo
   enddo
   r32 = sqrt(r32)
   r21 = sqrt(r21)

   kr12r23 = k_neb*(r32-r21)

   do iatom = 1, natoms
    do ix = 1,3
! the spring force component
     Fs(ix,iatom) = kr12r23*tang(ix,iatom)
    enddo
   enddo



!   dot = 0.0d0
!   do iatom = 1, natoms
!    do ix = 1,3
!     r32 = abs(r3(ix,iatom) - r2(ix,iatom))
!     r21 = abs(r2(ix,iatom) - r1(ix,iatom))
!     dot = dot + k_neb*(r21 - r32)*tang(ix,iatom)
!    enddo
!   enddo


! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
100 format(i3,3f10.4,'  ',3f10.4)
200 format(i3,3f8.4)

   return

 end subroutine get_Fspring


!############################################
! subroutine calcs projected total forces
!############################################
 subroutine get_Frec (natoms, image)

   use neb

   implicit none

! Argument Declaration and Description
! ===========================================================================
   integer, intent(in)                              :: natoms
   integer, intent(in)                              :: image

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer ix

   real, dimension (3,natoms)  :: F
   real dot
   real dot1

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Auxilliary parameters

   do iatom = 1, natoms
    F(:,iatom) = ftot_neb(:,iatom,image)
   enddo

   dot = 0.0d0
   do iatom = 1, natoms
    do ix = 1,3
     dot = dot + F(ix,iatom)*tang(ix,iatom)
    enddo
   enddo
   do iatom = 1, natoms
    do ix = 1,3
! the perpendicular total force component
     Frec(ix,iatom) = F(ix,iatom) - dot*tang(ix,iatom)
    enddo
   enddo

! test orthogonality
   dot1 = 0.0d0
   do iatom = 1, natoms
    do ix = 1,3
     dot1 = dot1 + Frec(ix,iatom) * dot*tang(ix,iatom)
    enddo
   enddo
   if (abs(dot1) .gt. 0.000000001d0 ) then
    write (*,*) 'ERROR dot1 =',dot1
    stop
   endif


! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
100 format(i3,3f8.4,' ',3f8.4,' ',3f7.4)
200 format(i3,3f8.4)

   return
 end subroutine get_Frec

!################################################
! subroutine returns NEB forces acting on image
!################################################
 subroutine get_Fneb (natoms, image, fms)

   use neb
   use fragments

   implicit none

! Argument Declaration and Description
! ===========================================================================
   integer, intent(in)                              :: natoms
   integer, intent(in)                              :: image
   real, intent(out)				    :: fms

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer ix
   real fxx
   real fmax

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Auxilliary parameters

   fms = 0.0d0
   fmax = 0.0d0
   write (*,*) ' iatom        F_neb          F_rec          Fs'

   if(allocated(fragxyz)) then
! Loop over atoms
    do iatom = 1, natoms
     do ix = 1,3
! project out fixed fragments
      if(fragxyz(ix,iatom) .eq. 1 ) then
       fxx = 0.0d0
      else
       fxx = Fs(ix,iatom) + Frec(ix,iatom)
       fms = max(fms,abs(Frec(ix,iatom)))
       fmax = max(fmax,abs(fxx))
      endif

! the spring force component
      Fneb(ix,iatom,image) = fxx
     enddo ! do ix
     write (*,100) iatom,Fneb(:,iatom,image),Frec(:,iatom),Fs(:,iatom)
    enddo ! do iatom

   else

     do iatom = 1, natoms
     do ix = 1,3
! project out fixed fragments
      fxx = Fs(ix,iatom) + Frec(ix,iatom)
      fms = max(fms,abs(Frec(ix,iatom)))
      fmax = max(fmax,abs(fxx))
! the spring force component
      Fneb(ix,iatom,image) = fxx
     enddo ! do ix
     write (*,100) iatom,Fneb(:,iatom,image),Frec(:,iatom),Fs(:,iatom)
    enddo ! do iatom

   endif ! if (allocated)

   write (*,*) 'Image no.',image,' Max. neb force =',fmax
! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
100 format(i3,3f8.4,' ',3f8.4,' ',3f7.4)

   return
 end subroutine get_Fneb


!#################################################
! subroutine moves image according to NEB forces
!#################################################

 subroutine move_image (image, rms)

! subroutine moves image according to NEB forces  using Verlet algorithm
! and filters out velocities pointing in the oposite direction to the NEB forces

   use neb
   use charges
   use interactions
   use configuration

   implicit none

! Argument Declaration and Description
! ===========================================================================
   integer, intent(in)                              :: image
   real, intent(out)                                :: rms

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer ix

   real  vdotF
   real  norm
   real  vmax
   real, dimension (3,natoms)     :: vel
   real, dimension (3,natoms)     :: pos
   real, dimension (3,natoms)     :: uvec
   real, dimension (3,natoms)     :: f

   character(2)    :: idn
   character(20)   :: fname

! Allocate Arrays
! ===========================================================================

! Procedure
! ===========================================================================
! Auxilliary parameters

! do local copies
   do iatom = 1,natoms
! copy positon of the image
    pos (:,iatom) = ratom_neb (:,iatom,image)
! copy velocities of the image
    vel (:,iatom) = vatom_neb (:,iatom,image)
    f (:,iatom) = Fneb(:,iatom,image)
   enddo

! keep the velocities parallel to the NEB forces
! if velocities points in a oposite direction
! we set the velocity component to zero
   norm = 0.0d0
   do iatom = 1,natoms
    do ix = 1,3
     norm = norm + f(ix,iatom)**2
    enddo
   enddo

   norm = sqrt(norm)
   vdotF = 0.0d0
   do iatom = 1,natoms
    do ix = 1,3
     uvec(ix,iatom) = f(ix,iatom)/norm
     vdotF = vdotF + vel(ix,iatom)*uvec(ix,iatom)
    enddo
   enddo
   write (*,*) 'image = ',image,' vdotF = ',vdotF
   if ( vdotF .gt. 0.0d0) then
    do iatom = 1,natoms
     vel(:,iatom) = vdotF*uvec(:,iatom)
    enddo
   else
    write (*,*) '-NEB- reset velocities for image no.',image
    vel(:,:) = 0.0d0
   endif

! do Verlet algorithm to update position and velocities
   call vverlet (natoms, pos, vel, f, xmass, dt_neb, image)



! quench velocities if total energy increased
!   if (detot_neb(image) .lt. 0.0d0) then
!    write (*,*) '-NEB- quench velocities for image no.',image
!    vel(:,:) = 0.0d0
!   endif

!   write (*,*) 'image = ',image,' VELOCITIES :'
   rms = 0.0d0
   vmax = 0.0d0
! save actual position and velocities for future
   do iatom = 1,natoms
    do ix = 1,3
     rms = max(rms,abs(ratom_neb (ix,iatom,image)-pos (ix,iatom)))
     vmax = max(vmax,abs(vel(ix,iatom)))
    end do
!    write (*,100) iatom,vel(:,iatom)
! copy positon of the image
    ratom_neb (:,iatom,image) = pos (:,iatom)
! copy velocities of the image
    vatom_neb (:,iatom,image) = vel (:,iatom)
   enddo
   write (*,*) 'Image no.',image,' Max. velocity = ',vmax

! dump actual atomic configuration into the *.xyz and *.bas file
   write (idn,'(i2.2)') image
! open *.xyz file
   fname = 'answer_neb_'//idn//'.xyz'
   open (unit= 90, file=fname, status='unknown', position = 'append')
! open *.bas file
   fname = 'answer_neb_'//idn//'.bas'
   open (unit= 91, file=fname, status='unknown')

! print header of *.xyz file
   write (90,*) natoms
   write (90,200) image, iter_neb, etot_neb(image)

! print header of *.bas file
   write (91,*) natoms

! print atoms
   do iatom = 1,natoms
    write (90,800) symbol(iatom), pos(:,iatom) - shifter(:)
    write (91,900) nzx(imass(iatom)), pos(:,iatom) - shifter(:)
   end do

! close files
   close(unit=90)
   close(unit=91)

! Deallocate Arrays
! ===========================================================================

! Format Statements
! ===========================================================================
100 format (i3,3f8.4)
200 format ('# image no. ',i2,' NEB iteration no.',i3,' Energy ',f20.6)
800 format (2x, a2, 3(2x,f12.6))
900 format (2x, i2, 3(2x,f12.6))

   return

 end subroutine move_image






!############################################
! velocity Verlet algorithm for n time step
!############################################
 subroutine vverlet(natoms, pos, vel, force, mass, dt, image)

   implicit none
! Argument Declaration and Description
! ===========================================================================
   real, dimension (3, natoms), intent(inout)		:: pos	 ! atom position
   real, dimension (3, natoms), intent(inout)		:: vel	 ! velocity
   real, dimension (3, natoms), intent(inout)		:: force ! force
   real, dimension (natoms), intent(in)			:: mass	 ! atom mass
   real, intent(in)		        		:: dt	 ! time step
   integer, intent(in)		        		:: natoms ! number of atoms
   integer, intent(in)		        		:: image ! number of atoms

! Local Variable Declaration and Description
! ===========================================================================
   integer :: iatom	! iteration
   integer :: ix
   real :: dt2		! aux time step
   real :: dt22		! aux time step
   real :: rmax
   real :: drr

!   real, parameter ::  dmax = 0.1d0

! Procedure
! ===========================================================================
! define aux time steps
  dt2 = dt*0.5d0
  dt22 = dt*dt*0.5d0

  rmax = 0.0d0
! loop over atoms
  do iatom = 1, natoms
! update velocity
   vel(:,iatom) = vel(:,iatom) + dt2*force(:,iatom)/mass(iatom)
! update position
   pos(:,iatom) = pos(:,iatom) + dt*vel(:,iatom)                     &
                      + dt22*force(:,iatom)/mass(iatom)
   do ix = 1,3
    drr = (dt*vel(ix,iatom) + dt22*force(ix,iatom)/mass(iatom))
    rmax = max(abs(drr),rmax)
   enddo
! update velocity (half step)
   vel(:,iatom) = vel(:,iatom) + (dt2*force(:,iatom)/mass(iatom))
  end do ! do iatom
  write (*,*) 'Image no.',image,' Max. displ= ',rmax

! Format Statements
! ===========================================================================

  return

 end subroutine vverlet



!############################################
! subroutine dumping actual images
!############################################
 subroutine writeout_image (natoms, symbol, shifter)

   use neb
   use interactions
   use charges
   implicit none

! Argument Declaration and Description
! ===========================================================================
   integer, intent(in)                              :: natoms
   character (len=2), dimension (natoms),intent(in) :: symbol
   real, dimension(3), intent(in)                   :: shifter

! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer iimage
   integer in1

   real, dimension(3)  :: pos

! Procedure
! ===========================================================================
! save the actual atomic configuration into the *.xyz/*.bas file

   open (unit= 90, file='answer_image.xyz', status='unknown')
   open (unit= 91, file='answer_image.bas', status='unknown')

! write number of images into *.bas file
  write (91,*)  nimg_neb

! loop over images
   do iimage = 1,nimg_neb

! write number of atoms
    write (90,*) natoms
    write (91,*) natoms

! comment line in *.xyz
    write (90,100) iimage, iter_neb, etot_neb(iimage)

! loop over atoms
    do iatom = 1,natoms
     in1 = imass(iatom)
     pos(:) = ratom_neb(:,iatom,iimage) - shifter(:)
! wirte atomic position
     write (90,200) symbol(iatom), pos(:)
     write (91,300) nzx(in1), pos(:)
    end do ! do iatoms
   end do ! do iimage

! close files
   close(unit=90)
   close(unit=91)

! Format Statements
! ===========================================================================
100 format ('## image no. ',i2,' NEB iteration no.',i4,' Energy ',f20.6)
200 format (2x, a2, 3(2x,f14.8))
300 format (2x, i2, 3(2x,f14.8))

  return
 end subroutine writeout_image
