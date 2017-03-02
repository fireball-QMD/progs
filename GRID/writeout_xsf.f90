! copyright info:
!
!                             @Copyright 2005
!                     Fireball Enterprise Center, BYU
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Institute of Physics, Czech Republic - Pavel Jelinek
!
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
 
! writeout_xsf.f90
! Program Description
! ===========================================================================
!       The subroutine writes out data in xsf-format file used by xcrysden
!
! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162  
! FAX +420-2-33343184
! Office telephone  +420-2-20318530
! email: jelinekp@fzu.cz
!
! ===========================================================================
!
! Program Declaration
! ===========================================================================
 subroutine writeout_xsf (xsfname, message, aa)
 
   use configuration
   use grid
   use charges
   use interactions
   use options
   implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input

   real, dimension (:), pointer, intent (in) :: aa
   character (len=40) xsfname
   character (len=30) message
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer i
   integer j
   integer k
   integer i0
   integer j0
   integer k0
   integer index
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! write out dden into *.xsf file (format of xcrysden visual code)
!  for details see www.xcrysden.org

   write (*,*) '  Write out xsf file ',xsfname
   write (*,100)
! ope file
   open ( unit = 302, file = xsfname, status = 'unknown' )
! print the list of atoms
!   if (icluster .eq. 1) then

!   write (302,*) 'ATOMS'
!   do iatom = 1,natoms
!    write (302,'(i2,3f14.8)') nzx(imass(iatom)),(ratom2g(i,iatom),i=1,3)
!   enddo

!  else

   write (302,*) 'CRYSTAL'
   write (302,*) 'PRIMVEC'
   write (302,*) (a1vec(i),i=1,3)
   write (302,*) (a2vec(i),i=1,3)
   write (302,*) (a3vec(i),i=1,3)

   write (302,*) 'PRIMCOORD'
   write (302,*) natoms,1
   do iatom = 1,natoms
    write (302,'(i2,3f14.8)') nzx(imass(iatom)),(ratom2g(i,iatom),i=1,3)
   enddo

!  endif

  write (302,*)
  write (302,*) 'BEGIN_BLOCK_DATAGRID_3D'
  write (302,*) message
  write (302,*) 'DATAGRID_3D_DENSITY'
  write (302,*) rm1+1, rm2+1, rm3+1
! print origin of the grid
  write (302,*) 0.0d0, 0.0d0, 0.0d0
! print lattice vector
  write (302,*) (a1vec(i),i=1,3)
  write (302,*) (a2vec(i),i=1,3)
  write (302,*) (a3vec(i),i=1,3)

! print values of the grid point
  do k = 0, rm3
    if (k .eq. rm3) then
     k0 = 0
    else
     k0 = k
    endif
    do j = 0, rm2
      if (j .eq. rm2) then
       j0 = 0
      else
       j0 = j
      endif
      do i = 0, rm1
        if (i .eq. rm1) then
         i0 = 0
        else
         i0 = i
        endif
! mapping index within the regular mesh
        index = i0 + rm1*j0 + rm1*rm2*k0
        write (302,200) aa(index)
      enddo ! do i
    enddo ! do j
  enddo ! do k
  write (302,*) 'END_DATAGRID_3D'
  write (302,*) 'END_BLOCK_DATAGRID_3D'
! close file
  close (302)
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (70('='))
200     format (2e16.8) 

        return
      end subroutine writeout_xsf

      
 subroutine writeout_cxsf (xsfname, message, aa)
 
   use configuration
   use grid
   use charges
   use interactions
   use options
   implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input

   complex, dimension (:), pointer, intent (in) :: aa
   character (len=40) xsfname
   character (len=30) message
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
   integer iatom
   integer i
   integer j
   integer k
   integer i0
   integer j0
   integer k0
   integer index
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! write out dden into *.xsf file (format of xcrysden visual code)
!  for details see www.xcrysden.org

   write (*,*) '  Write out xsf file ',xsfname
   write (*,100)
! ope file
   open ( unit = 302, file = xsfname, status = 'unknown' )
! print the list of atoms
   
   write (302,*) 'CRYSTAL'
   write (302,*) 'PRIMVEC'
   write (302,*) (a1vec(i),i=1,3)
   write (302,*) (a2vec(i),i=1,3)
   write (302,*) (a3vec(i),i=1,3)

   write (302,*) 'PRIMCOORD'
   write (302,*) natoms,1
   do iatom = 1,natoms
    write (302,'(i2,3f14.8)') nzx(imass(iatom)),(ratom2g(i,iatom),i=1,3)
   enddo


   write (302,*)
   write (302,*) 'BEGIN_BLOCK_DATAGRID_3D'
   write (302,*) message
   write (302,*) 'DATAGRID_3D_DENSITY'
   write (302,*) rm1+1, rm2+1, rm3+1
! print origin of the grid
   write (302,*) 0.0d0, 0.0d0, 0.0d0
! print lattice vector
   write (302,*) (a1vec(i),i=1,3)
   write (302,*) (a2vec(i),i=1,3)
   write (302,*) (a3vec(i),i=1,3)

! print values of the grid point
   do k = 0, rm3
    if (k .eq. rm3) then
     k0 = 0
    else
     k0 = k
    endif
    do j = 0, rm2
      if (j .eq. rm2) then
       j0 = 0
      else
       j0 = j
      endif
      do i = 0, rm1
        if (i .eq. rm1) then
         i0 = 0
        else
         i0 = i
        endif
! mapping index within the regular mesh
        index = i0 + rm1*j0 + rm1*rm2*k0
        write (302,200) aa(index)
      enddo ! do i
    enddo ! do j
   enddo ! do k
   write (302,*) 'END_DATAGRID_3D'
   write (302,*) 'END_BLOCK_DATAGRID_3D'
! close file
   close (302)
! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (70('='))
200     format (2e16.8) 

        return
      end subroutine writeout_cxsf
      
