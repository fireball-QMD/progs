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

! psi2es.f90
! Program Description
! ===========================================================================
!       This routine projects the actual time-dependent wave-function in
! MO basis set onto eigenfunction for given Rn (ionic position);
! diagonalization of current Hamiltonian should be done before
!
! ===========================================================================
! Code rewritten by:
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
        subroutine get_psi2es (itime_step)

        use tdse
        use density
        use dimensions
        use interactions
        use neighbor_map
        use configuration
        use kpoints
        use options

        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
       integer, intent (in) :: itime_step

! Local Parameters and Data Declaration
! ===========================================================================
       real, parameter :: trash = 0.05d0

       integer  ielec
       integer  iatom
       integer  jatom
       integer  ineigh
       integer  iband
       integer  imu
       integer  jmu
       integer  inu
       integer  jnu
       integer  ikpoint
       integer  in1
       integer  in2
       complex  a0
       complex  phi
       complex  proj
       complex, dimension (norbitals) :: psitmp
       real, dimension (norbitals) :: rpsi
       real tmp

       character(4)    :: ext1
       character(3)    :: ext2
       character(80)   :: fname


! Local Variable Declaration and Description
! ===========================================================================


! Procedure
! ===========================================================================
! Initialize variables
        a0 = cmplx(0.0d0,0.0d0)

! Loop over k-points
        do ikpoint = 1, nkpoints
         write (*,200) ikpoint
! Loop over electrons
         do ielec = 1, nelec

! make a copy of electron wavefunction
          psitmp(:) = psi(:,ielec,ikpoint)

! project ielec onto n-th eigenstate
          if (icluster .eq. 1) then
           do iband = 1,norbitals_new
            proj = a0
            do imu = 1, norbitals
             phi = cmplx(blowre(imu,iband,ikpoint),0.0d0)
             proj = proj + phi*psitmp(imu)
            enddo ! imu
            psi2es(iband,ielec,ikpoint) = proj
! writeout control output
            write (*,210) ielec,iband,psi2es(iband,ielec,ikpoint)
           enddo ! iband
          else
           do iband = 1,norbitals_new
            proj = a0
            do imu = 1, norbitals
! FIX later: check the sign of imaginary part maybe is wrong???
             phi = cmplx(blowre(imu,iband,ikpoint),-1.0d0*blowim(imu,iband,ikpoint))
             proj = proj + phi*psitmp(imu)
            enddo ! imu
            psi2es(iband,ielec,ikpoint) = proj
! writeout control output
            write (*,210) ielec,iband,psi2es(iband,ielec,ikpoint)
           enddo ! iband
          endif ! if (icluster)
          write (*,*)
         enddo ! ielec
         write (*,100)
        enddo ! ikpoints

! write out projection into a file

! open output files: first time?
        if(isp2es) then
! files do exist, we append
! psi
         open (unit = 18, file = 'psi2es.dat', status = 'unknown',        &
     &             position = 'append')

        else
! file does not exist, we start from the scratch
         open (unit = 18, file = 'psi2es.dat', status = 'unknown')
        endif

! loop over k-points
        do ikpoint = 1, nkpoints

         write (18,200) ikpoint

! Loop over electrons
         do ielec = 1, nelec
          write (18,201) ielec
! Loop over bands
          do iband = 1, norbitals_new
! filter out small contributions
           if (abs(real(psi2es(iband,ielec,ikpoint))) .gt. trash) then
            tmp = real(psi2es(iband,ielec,ikpoint)*conjg(psi2es(iband,ielec,ikpoint)))
            write (18,220) iband,tmp
           endif
          enddo ! do iband

! ==============================
! write out psi-k
! ==============================
          write (ext1,'(i4.4)') ielec
          write (ext2,'(i3.3)') ikpoint
          fname = 'psi'//ext1//'_k'//ext2//'.dat'
! open output files: first time?
          if(isp2es) then
           open (unit = 21, file = fname, status = 'unknown',        &
     &             position = 'append')
          else
           open (unit = 21, file = fname, status = 'unknown')
          endif ! if(isp2es)
          do iband = 1, norbitals
!           rpsi(imu) = (real(psi2es(imu,ielec,ikpoint)))**2 + (imag(psi2es(imu,ielec,ikpoint)))**2
           rpsi(iband) =  real(psi2es(iband,ielec,ikpoint)*conjg(psi2es(iband,ielec,ikpoint)))
          enddo
          write (21,401, advance="no") itime_step 
          do iband=1,norbitals_new
             write (21,400, advance="no") rpsi(iband)
          enddo
             write (21,*) ' '
          close (21)

         enddo ! do ielec

! open eigenstate-k file
! create an unique filename for each k-point
! e-k
         write (ext2,'(i3.3)') ikpoint
         fname = 'ES-t_k'//ext2//'.dat'
! open output files: first time?
         if(isp2es) then
          open (unit = 20, file = fname, status = 'unknown',        &
     &             position = 'append')
         else
          open (unit = 20, file = fname, status = 'unknown')
         endif ! if(isp2es)

! write out eigenvalues
         write (20,400) itime_step,(eigen_k(iband,ikpoint),iband=1,norbitals_new)
         close (20)

        enddo ! do ikpoint
! close file
        close (18)

! indicate that output files do exist for next call
        isp2es = .true.

! Deallocate Arrays
! ===========================================================================


! Format Statements
! ===========================================================================
100     format (2x, 70('='))
200     format (1x,'k-point no.',i4)
201     format (3x,'electron no.',i4)
210     format (3x,'projection of ',i4,'-e onto',i4,'-band =',2f10.4)
220     format (3x,i4,'-band =',f10.4)     
400     format (f12.6)
401     format (i6)


        return
        end subroutine get_psi2es

