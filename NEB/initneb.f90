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

! initneb.f90
! Program Description
! ===========================================================================

!
! Program Declaration
! ===========================================================================
        subroutine initneb (natoms, nspecies, imass, nzx, ratom) 
        use dimensions
        use neb
        use barrier
        use options
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: natoms
        integer, intent (in) :: nspecies

        integer, intent (in), dimension (natoms) :: imass
        integer, intent (in), dimension (nspec_max) :: nzx
	real, intent (inout), dimension (3,natoms) :: ratom

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer ispec
        integer i
        integer image_count
        integer natoms_in
        integer nimg_neb_in
        integer nucz
        integer irestart

        integer, dimension (natoms) :: imass_in
	real, dimension (3,natoms) :: dr

        character (len=30) file_neb

! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Initialize all of the positions of the final structure to zero.

        if (ineb .eq. 1) then 
         write (*,*) '  '
         write (*,100) 
         write (*,*) ' You have chosen to calculate a minimal energy '
         write (*,*) ' path based on Nugged Elastic Band (NEB) method '
         write (*,*) ' between an initial configuration and a final '
         write (*,*) ' configuration. '

         open (unit = 79, file = 'neb.optional', status = 'old')

         write (*,*) '  '
         write (*,*) ' An amount of spring constant applied between '
         write (*,*) ' images ' 
         write (*,*) ' Insert the the spring constant of the "chain": '
         read (79,*) k_neb
         write (*,*) ' k_neb = ', k_neb
        
         write (*,*) '  ' 
         write (*,*) ' The convergence criteria of neb procedure ' 
         write (*,*) ' is define by a set of tolerance parameters '
         write (*,*) ' Please, insert the desired tolerance of    '
         write (*,*) ' the maximal displacement :' 
         read (79,*) tol_displ_neb
         write (*,*) ' tol_displ_neb = ', tol_displ_neb

         write (*,*) '  ' 
         write (*,*) ' The convergence criteria of NEB procedure ' 
         write (*,*) ' is define by a set of tolerance parameters '
         write (*,*) ' Please, insert the desired tolerance of    '
         write (*,*) ' the NEB force :' 
         read (79,*) tol_ftot_neb
         write (*,*) ' tol_ftot_neb = ', tol_ftot_neb

         write (*,*) '  ' 
         write (*,*) ' The convergence criteria of NEB procedure ' 
         write (*,*) ' is define by a set of tolerance parameters '
         write (*,*) ' Please, insert the desired tolerance of    '
         write (*,*) ' the total energy :' 
         read (79,*) tol_etot_neb
         write (*,*) ' tol_etot_neb = ', tol_etot_neb

         write (*,*) '  '
         write (*,*) ' Number of images you will use to map the minimum '
         write (*,*) ' energy path (including the initial and final fixed '
         write (*,*) ' configurations).'
         read (79,*) nimg_neb
         write (*,*) ' nimg_neb = ', nimg_neb

         write (*,*) '  '
         write (*,*) ' Max. number of NEB iterations '
         read (79,*) niter_neb_max
         write (*,*) ' niter_neb_max = ', niter_neb_max

         write (*,*) '  '
         write (*,*) ' Time step used to optimize images atomic '
         write (*,*) ' configuration. As minimizatio technique we use'
         write (*,*) ' the velocity Verlet algorithm configurations.'
         read (79,*) dt_neb
         write (*,*) ' dt_neb = ', dt_neb

         write (*,*) '  '
         write (*,*) ' Do you want restart calculations? '
         write (*,*) ' no = 0; yes = 1 '
         read (79,*) irestart 
         if (irestart .eq. 1) write (*,*) ' NEB restart : YES'
         if (irestart .ne. 1) write (*,*) ' NEB restart : NO'

! RESTART
         if (irestart .eq. 1) then

          write (*,*) '  '
          write (*,*) ' Insert the filename with the images configuration: '
          read (79,200) file_neb
          write (*,201) 'filename =', file_neb
! open file 
          open (unit = 80, file = file_neb, status = 'old')
          read (80, *) nimg_neb_in
          if (nimg_neb_in .ne. nimg_neb) then
           write (*,*) ' Sorry nimg_neb_in = ', nimg_neb_in, ' and nimg_neb = ', nimg_neb
           write (*,*) ' Please check it out! '
           stop
          end if
         
! allocate arrays
          allocate (ratom_neb (3,natoms,nimg_neb)) 
          allocate (ftot_neb (3,natoms,nimg_neb))
          allocate (vatom_neb (3,natoms,nimg_neb))
          allocate (tang (3,natoms))
          allocate (Fs (3,natoms))
          allocate (Frec (3,natoms))
          allocate (Fneb (3,natoms,nimg_neb))
          allocate (etot_neb (nimg_neb))
          allocate (detot_neb (nimg_neb))
          etot_neb(:) = 0.0d0
          vatom_neb(:,:,:) =  0.0d0

! Lopp over images
          do image_count = 1, nimg_neb
! Loop over the number of atoms
           read (80, *) natoms_in
           if (natoms_in .ne. natoms) then
            write (*,*) ' Sorry natoms_in = ', natoms_in, ' and natoms = ', natoms
            write (*,*) ' The initial and final configuration files are not the '
            write (*,*) ' same structures - please check! '
            stop
           endif ! if (natoms_in) 
           do iatom = 1, natoms
            read (80,*) nucz, ratom_neb(:,iatom,image_count)
            do ispec = 1, nspecies
             if (nucz .eq. nzx(ispec)) imass_in(iatom) = ispec
            end do
            if (imass_in(iatom) .ne. imass(iatom)) then
             write (*,*) ' iatom = ', iatom
             write (*,*) ' Sorry imass_in(iatom) = ', imass_in(iatom)
             write (*,*) ' and imass(iatom) = ', imass(iatom)
             write (*,*) ' The initial and final configuration files are not the '
             write (*,*) ' same structures - please check! '
             write (*,*) ' The same ordering must be used in the two files. '
             stop
            end if
           end do ! do iatom
          end do ! do image_count

! end RESTART

         else 
! FROM SCRATCH

          write (*,*) '  '
          write (*,*) ' Insert the filename with the final configuration: '
          read (79,200) file_neb
          write (*,201) 'filename =', file_neb
! open file 
          open (unit = 80, file = file_neb, status = 'old')
          read (80, *) natoms_in
          if (natoms_in .ne. natoms) then
           write (*,*) ' Sorry natoms_in = ', natoms_in, ' and natoms = ', natoms
           write (*,*) ' The initial and final configuration files are not the '
           write (*,*) ' same structures - please check! '
           stop
          end if

! Loop over the number of atoms
          do iatom = 1, natoms
           read (80,*) nucz, ratom_final(:,iatom)
           do ispec = 1, nspecies
            if (nucz .eq. nzx(ispec)) imass_in(iatom) = ispec
           end do
           if (imass_in(iatom) .ne. imass(iatom)) then
            write (*,*) ' iatom = ', iatom
            write (*,*) ' Sorry imass_in(iatom) = ', imass_in(iatom)
            write (*,*) ' and imass(iatom) = ', imass(iatom)
            write (*,*) ' The initial and final configuration files are not the '
            write (*,*) ' same structures - please check! '
            write (*,*) ' The same ordering must be used in the two files. '
            stop
           end if
          end do ! do iatom

! allocate arrays
          allocate (ratom_neb (3,natoms,nimg_neb)) 
          allocate (ftot_neb (3,natoms,nimg_neb))
          allocate (vatom_neb (3,natoms,nimg_neb))
          allocate (tang (3,natoms))
          allocate (Fs (3,natoms))
          allocate (Frec (3,natoms))
          allocate (Fneb (3,natoms,nimg_neb))
          allocate (etot_neb (nimg_neb))
          allocate (detot_neb (nimg_neb))
          etot_neb(:) = 0.0d0
          vatom_neb(:,:,:) =  0.0d0

! copy the initial state
          do iatom = 1, natoms
           ratom_neb(:,iatom,1) = ratom(:,iatom)
          end do

! copy the final state
          do iatom = 1, natoms
           ratom_neb(:,iatom,nimg_neb) = ratom_final(:,iatom)
          end do

! construct an initial guess of the minimal energy path
! doing linear interpolation between the initial and final state
          do iatom = 1,natoms
	   dr(:,iatom) =                                               &
            (ratom_final(:,iatom)-ratom(:,iatom))/real(nimg_neb-1)    
          end do
        
          do i = 2, nimg_neb-1
           do iatom = 1,natoms
            ratom_neb(:,iatom,i) =  ratom_neb(:,iatom,i-1) + dr(:,iatom)
           end do ! do iatom 
          end do ! do i

         end if ! if (irestart .eq. 1) 



! set NEB iteration
         iter_neb = 1

! we do first image, which is allowed to relax (e.g. 2nd image)
         do iatom = 1, natoms
          ratom(:,iatom) = ratom_neb(:,iatom,1)
         enddo

! control output of the intial image guess 
         do image_count = 1, nimg_neb
          write (200,*) natoms
!          write (200,300) image_count
          do iatom = 1,natoms
           ispec = imass(iatom)
           write (200,301) nzx(ispec),ratom_neb(:,iatom,image_count)
          end do ! do iatom 
         end do ! do image_count

        end if ! if (ineb .eq. 1)

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
200     format (a20)
201     format (2x, a30)
300     format ('## image = ',i2)
301     format (i2,3f12.6)
 
        return
        end
