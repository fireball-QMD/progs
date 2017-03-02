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

! readsa.f90
! Program Description
! ===========================================================================
!       This routine reads the simulated annealing files.
!
! ===========================================================================
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
        subroutine readsa(natoms,isannealing,kseed,iall,tsa,dxmax,  &
     &                    Etotold,etotnew,ianneal,bsa,b)
        use dimensions
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent(in) :: natoms
        integer, intent(in) :: isannealing
        real, intent(in) :: b(3,natoms)
 
! Output
        integer, intent(out) :: kseed
        integer, intent(out) :: iall
        integer, intent(out) :: ianneal(natoms)
 
        real, intent(out) :: tsa
        real, intent(out) :: dxmax
        real, intent(out) :: etotnew
        real, intent(out) :: etotold
        real, intent(out) :: bsa (3,natoms)
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer i,j,ix
 
! Procedure
! ===========================================================================
        write(*,*) '  '
        etotnew=10000000.0e+20
        etotold=10000000.0e+20
        do i=1,natoms
         do ix=1,3
          bSA(ix,i)=b(ix,i)
         end do
        end do
! Skip this if no simulated annealing.
        if(isannealing.eq.1)then
! OK, we are simulated annealing.
         write(*,*) '  '
         write(*,*) ' We fill in details for the simulated annealing.'
         write(*,*) ' We have',natoms,' atoms.'
         write(*,*) '  '
         write(*,*) ' Which atoms are annealed is determined by the'
         write(*,*) ' input file SA.dat. '
         open(24,file='SA.optional',status='old')
         write(*,*) '  '
         write(*,*) ' Insert tsa, the simulated annealing temperature.'
         read(24,*)tsa
         write(*,*) ' tsa (Kelvins)=',tsa
         write(*,*) ' '
         write(*,*) ' Insert dxmax: the maximum displacement per step'
         write(*,*) ' in the x(or y or z) directions. The actual'
         write(*,*) ' displacement will be randomly generated between'
         write(*,*) ' [+dxmax, -dxmax].'
         read(24,*)dxmax
         write(*,*) ' dxmax=',dxmax
         write(*,*) '  '
         write(*,*) ' The next line is: iall. If '
         write(*,*) ' iall=1, then all are annealed. If iall.ne.1,'
         write(*,*) ' Then we read (1/line) i,ianneal, where i is atom number'
         write(*,*) ' and ianneal=1,0 for Y/N anneal.'
         read(24,*)iall
         if(iall.eq.1)then
          do i=1,natoms
           ianneal(i)=1
          end do
         end if
         if(iall.ne.1)then
          do i=1,natoms
           read(24,*)j,ianneal(i)
           if(j.ne.i)then
            write(*,*) ' bad i,j in SA.dat',i,j
            write(*,*) ' sorry... fix up SA.dat'
            stop
           end if
          end do
         end if
         close(unit=24)
!
         write(*,*) '  '
         do i=1,natoms
          write(*,*) ' Atom #',i,' ianneal=',ianneal(i)
         end do
         write(*,*) '  '
         write(*,*) ' Finally, we need a seed for the random number'
         write(*,*) ' generator. We call it kseed. Insert kseed.....'
         read(24,*)kseed
         write(*,*) ' kseed=',kseed
         close(unit=24)
        end if
 
! Format Statements
! ===========================================================================
 
        return
        end
