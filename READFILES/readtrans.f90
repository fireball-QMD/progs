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
!
! readscript.f90
! Program Description
! ===========================================================================
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine readtrans (natoms)

   use interactions
   use transport
   implicit none

! Argument Declaration and Description
! ===========================================================================
   integer, intent (in)             :: natoms

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================

   integer i
   integer j
   integer iatom
   integer jatom
   integer nint
   integer imu
   integer inu

   real delta

! Procedure
! ===========================================================================

   write (*,*) '  '
   write (*,*) '  '
   write (*,*) ' Welcome to readtrans subtoutine! '
   write (*,*) ' Read in the data from the scriptfile - trans.optional '

   open (unit = 12, file = 'trans.optional', status = 'old')


   write (*,*) ' Select if different imaginary part eta will be used. '
   read (12,*) ieta
   write (*,600) ieta
   write (*,*) ' Select if detail informations will be printed. '
   read (12,*) iwrt_trans
   write (*,700) iwrt_trans
   write (*,*) ' Select if analysis of channels will be performed. '
   read (12,*) ichannel
   write (*,800) ichannel
   write (*,*) ' Perform fitting of hoppings. '
   read (12,*) ifithop
   write (*,900) ifithop

   write (*,*) 'Warnning !!'
   write (*,*) 'Energy scale will be shifted with respect to '
   write (*,*) 'Fermi level equals to zero.  '
! energy range
   read (12,*) Elow
   read (12,*) Eup
   if(Elow .gt. Eup) then
    write (*,*)  ' Lower energy range is higher than Upper energy range!'
    write (*,400) Elow, Eup
    write (*,*)  ' Please, fix it and run again!'
    stop
   endif
   write (*,400) Elow, Eup
! number of the energy steps
   read (12,*) nE
   if (nE .eq. 0) then
    write (*,*) ' Number of the energy steps cannot be zero!!!'
    stop
   endif
   write (*,*) 'nE =',nE
! imaginary part of energy
   read (12,*) eta
! evaluate the energy step
   delta = (Eup - Elow)/real(nE)
! convert energy step into complex variable
   dE = delta*(1.0d0,0.0d0)
! resume
   write (*,401) real(dE), eta
   write (*,402) nE

   close (unit = 12)

   write (*,*) '  '
   write (*,*) '  '
   write (*,*) ' Welcome to readeta.f! '
   write (*,*) ' Read information about different imaginary part eta  '

! Allocate array for list of atoms
   allocate (ideta (natoms))
   ideta = 0
   eta0 = eta
   write (*,*) 'natoms =',natoms
   if (.not. ieta) then

      write (*,*) '  '
      write (*,*) '  No optional eta will be applied.'

   else

      open (unit = 17, file = 'eta.optional', status = 'old')
! Read the basis file
      write (*,100)
! Read value of alternative eta
      read (17,*) eta0
      write (*,*) '  eta0 = ',eta0
! read number of intervals
      read (17,*) nint
      write (*,*) ' nint = ',nint
      do inu = 1,nint
         read (17,*) iatom,jatom
! loop over the interval
         do imu = iatom,jatom
! control check
            if (imu .gt. natoms) then
               write (*,*) 'ERROR, atom ',imu,'exceed max. number of atoms.'
               stop
            endif
            ideta(imu) = 1
         enddo ! do imu
      enddo ! do inu

      write (*,*) 'Information about optional eta0 applied on the selected atoms'
      do inu = 1,natoms
         write (*,*) 'Atom: ',inu,'   eta: ',ideta(inu)
      enddo
! close the input file
      close (unit = 17)
   endif


   return

! Format Statements
! ===========================================================================
100     format (2x, 70('='))
200     format (a30)
600     format (2x,' ieta         = ',i2)
700     format (2x,' iwrtout      = ',i2)
800     format (2x,' ichannel     = ',i2)
900     format (2x,' ifithop      = ',i2)
400     format (2x,'Energy range = ',f16.8,' - ',f16.8,'  [eV]')
401     format (2x,'Energy step = ', f16.8,' eta = ', f16.8)
402     format (2x,'Number of energy step = ',i6)

 end subroutine readtrans
