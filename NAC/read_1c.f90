! copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Institute of Physics, Czech Republic - Pavel Jelinek

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

! read_1c.f90
! Program Description
! ===========================================================================
!       This routine reads in the 1-center (exchange-correlation)
! interactions. These 1-center interactions are contributions as described
! in the Horsefield approach.  This routine also reads in the variables which
! are needed to compute changes of the exchange correlation for the of charge
! transfer
!
! JOM : It also reads the one -center terms for non-adiabatic coupling
! ===========================================================================
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
! JOM-warning
! now the options itheory, itheory_xv, ispin are in MODULES/options.f90
!        subroutine read_1c (nspecies, itheory, itheory_xc, ispin,       &
!     &                      ioff1c)
        subroutine read_1c (nspecies, ioff1c)
        use charges
        use dimensions
        use integrals
        use interactions
        use options
        use nonadiabatic
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: ioff1c
!       integer, intent (in) :: ispin
!       integer, intent (in) :: itheory
!       integer, intent (in) :: itheory_xc
        integer, intent (in) :: nspecies
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iline
        integer in1, in2
        integer ins
        integer issh
        integer isorp
        integer itype
        integer jssh
        integer kssh
        integer kkssh
        integer numsh
 
        character (len=70) message
 
        logical skip_it

        real, dimension (nspecies) :: idshell

        integer, dimension (nsh_max) :: imask
        integer ideriv
        integer iissh, jjssh
 
! Allocate Arrays
! ===========================================================================

! Original statement
!        if( itheory_xc .eq. 0 ) then 

! NOTE for SNXC and OLSXC in Extended Hubbard theory. (OCT. 15. 04)
! For Extended Hubbard theory, one-center XC matrix elements(xcnu1c) were
! created only in the case of Horsfield(itheory_xc=0).
! To avoid segmentation fault, we tentatively use the same matrix elements
! for SNXC and OLSXC.

        if( itheory_xc .eq. 0 .or. itheory .eq.2) then
        allocate(exc1c_0 (nspecies,0:2))
        allocate(exc1c (nspecies,1:nsh_max,1:nsh_max,0:2))
        allocate(xcnu1c (nsh_max, nsh_max, nspecies))
        allocate(xcnu1cs (nsh_max, nsh_max, nspecies))
 
! Procedure
! ===========================================================================
! ***************************************************************************
!
!         H O R S E F I E L D    E X C H A N G E - C O R R E L A T I O N
!
! *************************************************************************** 

! Initialize to zero.
        exc1c_0 = 0.0d0
        exc1c =0.0d0
        xcnu1c = 0.0d0
        xcnu1cs = 0.0d0
! Read data from the 1-center exchange-correlation file.
         open (unit = 36, file = trim(fdataLocation)//'/xc_onecenter.dat', status = 'unknown')
         do iline = 1, 4
          read (36,100) message
         end do
 
         do in1 = 1, nspecies + nsup
          read (36,100) message
         end do
         read (36,100) message

         do isorp = 0, 2
          read (36,100) message
          in2 = 1
          do in1 = 1, nspecies + nsup
           skip_it = .false.
           do ins = 1, nsup
            if (nsu(ins) .eq. in1) skip_it = .true.
           end do
           if (skip_it) then
            read (36,*) itype, numsh
            read (36,100) message
            do jssh = 1, numsh
             read (36,100) message
            end do
           else
            read (36,*) itype, numsh
            if (numsh .ne. nssh(in2)) then
             write (*,*) ' numsh .ne. nssh in read_1c.f90 '
             write (*,*) itype, numsh, in1, nssh(in2)
             stop
            end if
            read (36,*) exc1c_0(in2,isorp)
            do jssh = 1, numsh
             read (36,*) (exc1c(in2,jssh,issh,isorp), issh = 1, numsh)
            end do
            exc1c_0(in2,isorp) = exc1c_0(in2,isorp)*ioff1c
            do jssh = 1, numsh
             do issh = 1, numsh
              exc1c(in2,jssh,issh,isorp) = exc1c(in2,jssh,issh,isorp)*ioff1c
             end do
            end do
            in2 = in2 + 1
           end if
          end do
         end do
         
         in2 = 1
         do in1 = 1, nspecies + nsup
          skip_it = .false.
          do ins = 1, nsup
           if (nsu(ins) .eq. in1) skip_it = .true.
          end do
         if (skip_it) then
          read (36,100) message
         else
          read (36,*) idshell(in2), dq(in2)
          in2 = in2 + 1
         end if
        end do
        close (unit = 36)
 
! Read data from the 1-center exchange-correlation extended hubbard file.
        if (itheory .eq. 2) then
         open (unit = 37, file = trim(fdataLocation)//'/nuxc_onecenter.dat',         &
     &         status = 'unknown')
         do iline = 1, 4
          read (37,100) message
         end do
         
         do in1 = 1, nspecies + nsup
          read (37,100) message
         end do
         read (37,100) message
         
         in2 = 1
         do in1 = 1, nspecies + nsup
          skip_it = .false.
          do ins = 1, nsup
           if (nsu(ins) .eq. in1) skip_it = .true.
          end do
          if (skip_it) then
           read (37,*) itype, numsh
           do issh = 1, numsh
            read (37,100) message
           end do
          else 
           read (37,*) itype, numsh
           if (itype .ne. in1 .or. numsh .ne. nssh(in2)) stop           &
     &      ' Problem reading nuxc_onecenter.dat file '
           do issh = 1, numsh
            read (37,*) (xcnu1c(issh,jssh,in2), jssh = 1, numsh)
           end do
           in2 = in2 + 1
          end if
         end do
         close (unit = 37)

! Read data from the 1-center spin exchange-correlation extended hubbard
! file.
          if (ispin .eq. 1) then
           open (unit = 38, file = trim(fdataLocation)//'/nuxcs_onecenter.dat',              &
     &           status = 'unknown')
           do iline = 1, 4
            read (38,100) message
           end do
 
           do in1 = 1, nspecies + nsup
            read (38,100) message
           end do
           read (38,100) message

           in2 = 1
           do in1 = 1, nspecies + nsup
            skip_it = .false.
            do ins = 1, nsup
             if (nsu(ins) .eq. in1) skip_it = .true.
            end do
            if (skip_it) then
             read (38,*) itype, numsh
             do issh = 1, numsh
              read(38,100) message
             end do
            else
             read (38,*) itype, numsh
             if (itype .ne. in1 .or. numsh .ne. nssh(in2)) stop              &
     &        ' Problem reading nuxcs_onecenter.dat  file '
             do issh = 1, numsh
              read (38,*) (xcnu1cs(issh,jssh,in2),jssh = 1, numsh)
             end do
             in2 = in2 + 1
            end if
           end do
           close (unit = 38)
          end if
         end if
        end if   ! end if (itheory_xc .eq. 0)
! ***************************************************************************
!
!            M c W E D A   E X C H A N G E - C O R R E L A T I O N
!
! *************************************************************************** 
        if (itheory_xc .eq. 2) then 
        
         allocate(exc1c0 (nspecies,nsh_max,nsh_max))
         allocate(nuxc1c (nspecies,nsh_max,nsh_max))
         allocate(dexc1c (nspecies,nsh_max,nsh_max,nsh_max))
         allocate(d2exc1c (nspecies,nsh_max,nsh_max))
         allocate(dnuxc1c (nspecies,nsh_max,nsh_max,nsh_max))
         allocate(d2nuxc1c (nspecies,nsh_max,nsh_max,nsh_max,nsh_max))
        
! jel-nac
!         if (imdet .eq. 1) then 
         if (imdet .ne. 0) then 
           allocate(f1nac1c (nspecies,nsh_max,nsh_max))
           allocate(f2nac1c (nspecies,nsh_max,nsh_max))
         end if

         open (unit = 36, file = trim(fdataLocation)//'/xc1c_dqi.dat', status = 'unknown')
        
         do iline = 1, 4
          read (36,100) message
         end do
        
         do in1 = 1, nspecies + nsup
          read (36,100) message
         end do
         read (36,100) message

         in2 = 1
         do in1 = 1, nspecies + nsup
          skip_it = .false.
          do ins = 1, nsup
           if (nsu(ins) .eq. in1) skip_it = .true.
          end do

          if (skip_it) then
           read (36,*) itype, numsh

           do issh = 1, numsh
            read (36,*)
           end do
           read (36,*)
           do issh = 1, numsh
            read (36,*)
           end do

          else

           read (36,*) itype, numsh
           if (numsh .ne. nssh(in2)) then
            write (*,*) ' numsh .ne. nssh in read_1c.f90 '
            write (*,*) itype, numsh, in1, nssh(in2)
            stop
           end if

!           read (36,*) exc1c0(in2,iss)
           do issh = 1, numsh
            read (36,*) (exc1c0(in2,issh,jssh),jssh=1,numsh)
           end do
           read (36,*)
           do issh = 1, numsh
            read (36,*) (nuxc1c(in2,issh,jssh),jssh=1,numsh)
           end do
!           exc1c0(in2) = exc1c0(in2)*ioff1c
           do issh = 1, numsh
            do jssh = 1, numsh
             nuxc1c(in2,issh,jssh) = nuxc1c(in2,issh,jssh)*ioff1c
             exc1c0(in2,issh,jssh) = exc1c0(in2,issh,jssh)*ioff1c
            end do
           end do

! End loop over 2. derivative (only non-diagonal elements d^2 / dqi dqj)  

! increment 'shadow' ispec counter
           in2 = in2 + 1
          end if ! if (skip_it)
         end do ! in1

! ==================================================================
!       READ FILE             nuxc1crho.dat
! ==================================================================
         open (unit = 36, file = trim(fdataLocation)//'/nuxc1crho.dat', status = 'unknown')

! Read header
         do iline = 1, 4
          read (36,100) message
         end do

         do in1 = 1, nspecies + nsup
          read (36,100) message
         end do
         read (36,100) message

         in2 = 1
! skip unsed species
         do in1 = 1, nspecies + nsup
          skip_it = .false.
          do ins = 1, nsup
           if (nsu(ins) .eq. in1) skip_it = .true.
          end do

          if (skip_it) then

           do 
            read (36,*) itype, numsh, kkssh
            do issh = 1, numsh
             read (36,*)
            end do
            if (numsh .eq. kkssh) exit
           end do ! do kssh 

          else

           do kssh = 1, nssh(in2)
            read (36,*) itype, numsh, kkssh
            if (numsh .ne. nssh(in2)) then
             write (*,*) ' numsh .ne. nssh in read_1c.f90 '
             write (*,*) itype, numsh, in1, nssh(in2)
             stop
            end if

            do issh = 1, numsh
             read (36,*) (dnuxc1c(in2,issh,jssh,kssh),jssh=1,numsh)
            end do

            do issh = 1, numsh
             do jssh = 1, numsh
              dnuxc1c(in2,issh,jssh,kssh) = dnuxc1c(in2,issh,jssh,kssh)*ioff1c
!              dexc1c(in2,issh,jssh,kssh) = exc1c0(in2,issh,jssh.kssh)*ioff1c
             end do
            end do
           end do ! do kssh
! increment 'shadow' ispec counter
           in2 = in2 + 1
          end if ! if (skip_it)
         end do ! in1

! ==================================================================
!       READ FILE             exc1crho.dat
! ==================================================================
         open (unit = 36, file = trim(fdataLocation)//'/exc1crho.dat', status = 'unknown')

! Read header
         do iline = 1, 4
          read (36,100) message
         end do

         do in1 = 1, nspecies + nsup
          read (36,100) message
         end do
         read (36,100) message

         in2 = 1
! skip unsed species
         do in1 = 1, nspecies + nsup
          skip_it = .false.
          do ins = 1, nsup
           if (nsu(ins) .eq. in1) skip_it = .true.
          end do

          if (skip_it) then
           do
            read (36,*) itype, numsh, kkssh
            do issh = 1, numsh
             read (36,*)
            end do
            if (numsh .eq. kkssh) exit
           end do ! do kssh

          else

           do kssh = 1, nssh(in2)
            read (36,*) itype, numsh, kkssh
            if (numsh .ne. nssh(in2)) then
             write (*,*) ' numsh .ne. nssh in read_1c.f90 '
             write (*,*) itype, numsh, in1, nssh(in2)
             stop
            end if

            do issh = 1, numsh
             read (36,*) (dexc1c(in2,issh,jssh,kssh),jssh=1,numsh)
            end do

            do issh = 1, numsh
             do jssh = 1, numsh
              dexc1c(in2,issh,jssh,kssh) = dexc1c(in2,issh,jssh,kssh)*ioff1c
             end do
            end do
           end do ! do kssh
! increment 'shadow' ispec counter
           in2 = in2 + 1
          end if ! if (skip_it)
         end do ! in1

         if (imdet .eq. 1) then
! ==================================================================
!       READ FILES             goverlapf1.dat and goverlapf2.dat
! ==================================================================
         open (unit = 36, file = trim(fdataLocation)//'/goverlapf1.dat', status = 'unknown')

! Read header
         do iline = 1, 4
          read (36,100) message
         end do

         do in1 = 1, nspecies + nsup
          read (36,100) message
         end do
         read (36,100) message

         in2 = 1
! skip unsed species
! JOM-warning : I am not sure that "skip" works well here
         do in1 = 1, nspecies + nsup
          skip_it = .false.
          do ins = 1, nsup
           if (nsu(ins) .eq. in1) skip_it = .true.
          end do

          if (skip_it) then

!          do 
!           read (36,*) itype, numsh, kkssh
            read (36,*) itype, numsh
            do issh = 1, numsh
             read (36,*)
            end do
!           if (numsh .eq. kkssh) exit
!          end do ! do kssh 

          else

!          do kssh = 1, nssh(in2)
!           read (36,*) itype, numsh, kkssh
            read (36,*) itype, numsh
            if (numsh .ne. nssh(in2)) then
             write (*,*) ' numsh .ne. nssh in read_1c.f90 '
             write (*,*) itype, numsh, in1, nssh(in2)
             stop
            end if

            do issh = 1, numsh
             read (36,*) (f1nac1c(in2,issh,jssh),jssh=1,numsh)
            end do

!            do issh = 1, numsh
!             do jssh = 1, numsh
!              dnuxc1c(in2,issh,jssh,kssh) = dnuxc1c(in2,issh,jssh,kssh)*ioff1c
!!              dexc1c(in2,issh,jssh,kssh) = exc1c0(in2,issh,jssh.kssh)*ioff1c
!             end do
!            end do
!          end do ! do kssh
! increment 'shadow' ispec counter
           in2 = in2 + 1
          end if ! if (skip_it)
         end do ! in1

         open (unit = 37, file = trim(fdataLocation)//'/goverlapf2.dat', status = 'unknown')

! Read header
         do iline = 1, 4
          read (37,100) message
         end do

         do in1 = 1, nspecies + nsup
          read (37,100) message
         end do
         read (37,100) message

         in2 = 1
! skip unsed species
! JOM-warning : I am not sure that "skip" works well here
         do in1 = 1, nspecies + nsup
          skip_it = .false.
          do ins = 1, nsup
           if (nsu(ins) .eq. in1) skip_it = .true.
          end do

          if (skip_it) then

!          do 
!           read (36,*) itype, numsh, kkssh
            read (37,*) itype, numsh
            do issh = 1, numsh
             read (37,*)
            end do
!           if (numsh .eq. kkssh) exit
!          end do ! do kssh 

          else

!          do kssh = 1, nssh(in2)
!           read (36,*) itype, numsh, kkssh
            read (37,*) itype, numsh
            if (numsh .ne. nssh(in2)) then
             write (*,*) ' numsh .ne. nssh in read_1c.f90 '
             write (*,*) itype, numsh, in1, nssh(in2)
             stop
            end if

            do issh = 1, numsh
             read (37,*) (f2nac1c(in2,issh,jssh),jssh=1,numsh)
            end do

!            do issh = 1, numsh
!             do jssh = 1, numsh
!              dnuxc1c(in2,issh,jssh,kssh) = dnuxc1c(in2,issh,jssh,kssh)*ioff1c
!!              dexc1c(in2,issh,jssh,kssh) = exc1c0(in2,issh,jssh.kssh)*ioff1c
!             end do
!            end do
!          end do ! do kssh
! increment 'shadow' ispec counter
           in2 = in2 + 1
          end if ! if (skip_it)
         end do ! in1
         end if
! ==================================================================
 end if ! if(itheory_xc.eq.2) 
        
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (a70)
200     format (2x, i3, 4x, i3, 2x, i3, 2x, i3, 2x, i3, 2x, i3, 2x, i3)
!300     format (<numsh>f16.6)       
        return
        end subroutine read_1c
      
