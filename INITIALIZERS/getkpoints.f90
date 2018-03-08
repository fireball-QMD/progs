! copyright info:
!
!                             @Copyright 2006
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
! fireball-qmd is a free (GPLv3) open project.

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.


! getkpoints.f90
! Program Description
! ===========================================================================
!       This routine either generates kpoints automatically or from a user-
!       defined mesh or reads them in from a file.
! ===========================================================================
! Code written by:
! j keith employed by
! James P. Lewis
! Department of Physics and Astronomy
! Brigham Young University
! N233 ESC P.O. Box 24658
! Provo, UT 84602-4658
! FAX (801) 422-2265
! Office Telephone (801) 422-7444
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine getkpoints (icluster, Vouc, a1vec, a2vec,  &
   &                          a3vec, lvsfile, basisfile, iquench, ireducekpts, rescal)
        use constants_fireball
        use dimensions
        use kpoints
        use options, only : verbosity
        implicit none

! Argument Declaration and Description
! ===========================================================================
! Input
        character(len = 40), intent(in) :: lvsfile
        character(len = 40), intent(in) :: basisfile
        integer, intent(in) :: icluster
        integer, intent(in) :: ireducekpts
        integer, intent(in) :: iquench
        real, intent(in), dimension (3) :: a1vec
        real, intent(in), dimension (3) :: a2vec
        real, intent(in), dimension (3) :: a3vec
        real, intent(in) :: Vouc
        real, intent(in) :: rescal
 
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        character (len = 40) kptsfile
        integer ki
        integer ikpoint
        integer, dimension(3) :: mpmesh
        real sum_weight
        real coef
        real, dimension(3,3) :: primitive_lattice, qvec
        real, dimension(3) :: qmag, amag 

! Procedure
! ===========================================================================

        if (icluster .ne. 0) then
          write (*,*) '  '
          write (*,*) ' You are doing a cluster calculation '
          write (*,*) ' k-points are not read in '
          write (*,*) '  '
          nkpoints = 1
          allocate (special_k(3, nkpoints))
          allocate (special_k_orig(3, nkpoints))
          allocate (scale_k(3, nkpoints))
          allocate (weight_k(nkpoints))
          allocate (weight_k_orig(nkpoints))
          special_k_orig(:,1) = 0
          weight_k_orig(1) = 1
          weight_k(1) = 1
          call initkpoints ()
          return
        end if


! read kpoints from file--useful for band structure calcs
        if(index(kptpreference,'.kpts')/=0) then
           kptsfile=kptpreference

! automatic kpoint generation
        else if(index(kptpreference,'a') .eq. 1) then 
print*, 'automatic kpoint generation--> the strategy is to evenly tile' 
print*, 'reciprocal space.  Thus k points are assigned based on reciprocal vector'
print*, 'length (and for illustration cubic real space length below) as so:'
print*, '#kpts Gm->  2-->  4-->  6-->  8-->  10->  12->  14->  16----->'
print*, 'qmag  0     0.31  0.63  1.38  2.13  2.88  3.63  4.38  5.13  '
print*, 'real  inf   20.0  10.0  4.55  2.95  2.18  1.73  1.43  1.22   '
print*, 'feel free to improve this algorithm or adjust kpoint assignments.'  
print*, 'It is only for semiconductors/insulators--metals probably need a'
print*, 'higher density of kpoints. Example: if uc is 15 x 15 x 15 Ang then'
print*, 'mpmesh is 2,2,2'
print*, 'Use this for quick calcs. Do kpoint convergence tests if serious.'
           ! first calculate reciprocal lattice vectors
           primitive_lattice(1,:)=a1vec
           primitive_lattice(2,:)=a2vec
           primitive_lattice(3,:)=a3vec
           coef=2.0*pi/Vouc
           qvec(1,1)=coef*(primitive_lattice(2,2)*primitive_lattice(3,3) &
     &     -primitive_lattice(2,3)*primitive_lattice(3,2))
           qvec(1,2)=coef*(primitive_lattice(2,3)*primitive_lattice(3,1) &
     &     -primitive_lattice(2,1)*primitive_lattice(3,3))
           qvec(1,3)=coef*(primitive_lattice(2,1)*primitive_lattice(3,2) &
     &     -primitive_lattice(2,2)*primitive_lattice(3,1))
           qvec(2,1)=coef*(primitive_lattice(3,2)*primitive_lattice(1,3) &
     &     -primitive_lattice(3,3)*primitive_lattice(1,2))
           qvec(2,2)=coef*(primitive_lattice(3,3)*primitive_lattice(1,1) &
     &     -primitive_lattice(3,1)*primitive_lattice(1,3))
           qvec(2,3)=coef*(primitive_lattice(3,1)*primitive_lattice(1,2) &
     &     -primitive_lattice(3,2)*primitive_lattice(1,1)) 
           qvec(3,1)=coef*(primitive_lattice(1,2)*primitive_lattice(2,3) &
     &     -primitive_lattice(1,3)*primitive_lattice(2,2))
           qvec(3,2)=coef*(primitive_lattice(1,3)*primitive_lattice(2,1) &
     &     -primitive_lattice(1,1)*primitive_lattice(2,3))
           qvec(3,3)=coef*(primitive_lattice(1,1)*primitive_lattice(2,2) &
     &     -primitive_lattice(1,2)*primitive_lattice(2,1))
           ! get their magnitudes
           qmag=(/sqrt(dot_product(qvec(1,:),qvec(1,:))),sqrt(dot_product( &
     &     qvec(2,:),qvec(2,:))),sqrt(dot_product(qvec(3,:),qvec(3,:)))/)
           do ki=1,3
              ! set the first two by 'good judgement'
              if(qmag(ki)<0.31)then
                 mpmesh(ki)=1
              elseif(qmag(ki).ge.0.31 .and. qmag(ki)<0.63)then
                 mpmesh(ki)=2
              ! then just add 2 points every 0.75 Ang^-1   
              elseif(qmag(ki).ge.0.63 .and. qmag(ki)<1.38)then
                 mpmesh(ki)=4
              elseif(qmag(ki).ge.1.38 .and. qmag(ki)<2.13)then
                 mpmesh(ki)=6
              elseif(qmag(ki).ge.2.13 .and. qmag(ki)<2.88)then
                 mpmesh(ki)=8
              elseif(qmag(ki).ge.2.88 .and. qmag(ki)<3.63)then
                 mpmesh(ki)=10
              elseif(qmag(ki).ge.3.63 .and. qmag(ki)<4.38)then
                 mpmesh(ki)=12
              elseif(qmag(ki).ge.4.38 .and. qmag(ki)<5.13)then
                 mpmesh(ki)=14
              ! max out at 16
              else
                 mpmesh(ki)=16
              end if
           end do
           call greatKAuto(lvsfile,basisfile,mpmesh,iquench,ireducekpts,nkpoints)
           kptsfile='greatK.kpts'

! user-specified Monkhorst Pack mesh
        else
            write (*,*) 'This option is not supported at the moment'
            write (*,*) ' Please, change the input and run it again.'
            stop
!           read(kptpreference,'(3i)')mpmesh
!           call greatKAuto(lvsfile,basisfile,mpmesh,iquench,ireducekpts,nkpoints)
!           kptsfile='greatK.kpts'
        end if

! read kpoints from *.kpts
        open (unit = 54, file = kptsfile, status = 'old')
        read (54,*) nkpoints 
        write (*,*) '  '
        write (*,*) ' Reading k-points '
        if (verbosity .ge. 3) write (*,100)
        if (verbosity .ge. 3) write (*,*) ' nkpoints = ', nkpoints
        if (nkpoints .le. 0) then
         write (*,*) ' nkpoints .le. 0 Huh! Fix XXX.kpts file! '
         stop
        end if
! Allocate the kpoints module memory
        allocate (special_k(3, nkpoints))
        allocate (special_k_orig(3, nkpoints))
        allocate (scale_k(3, nkpoints))
        allocate (weight_k(nkpoints))
        allocate (weight_k_orig(nkpoints))
        sum_weight = 0.0d0
        do ikpoint = 1, nkpoints
         read (54,*) special_k_orig(:,ikpoint), weight_k_orig(ikpoint)
         do ki=1,3
           special_k_orig(ki,ikpoint)=special_k_orig(ki,ikpoint)/rescal
         end do
         if (verbosity .ge. 3) write (*,300) ikpoint, special_k_orig(:,ikpoint),                   &
     &              weight_k_orig(ikpoint)
         sum_weight = sum_weight + weight_k_orig(ikpoint)
        end do
        write (*,100)
        if (abs(sum_weight - 1.0d0) .gt. 1.0d-3) then
         write (*,*) ' Sum of k-point weights = ', sum_weight
         write (*,*) ' They do not add to one (within 1.0d-3) '
         write (*,*) ' Sorry -- Stop --- fix XXX.kpts'
         write (*,*) ' The k-points sum of weights does not add to 1! '
         stop
        end if
        close (unit = 54)
! Initialize the kpoints from the special k-points read in from the
! *.kpt file. If this is a restart, then read k-points from a file which
! dumped the k-points from a previous run.
        call initkpoints ()
 
! Format Statements
! ===========================================================================
100     format (2x, 70('='))
300     format (2x, i4, ' - ikpoint = ', 3f11.6, 2x, ' weight = ', f9.4) 


        end subroutine getkpoints
