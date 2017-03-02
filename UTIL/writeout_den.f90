! iopyright info:
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

! writeout_den.f90
! Program Description
! ===========================================================================
!       This routine writes out the density projected on the grid.
!
! ===========================================================================
! Code written by:
! ===========================================================================
!
! Program Declaration
! ===========================================================================
        subroutine writeout_den (natoms, icluster, a1vec, a2vec, a3vec)
        use dimensions
        use configuration
        use grid
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: icluster
        integer, intent (in) :: natoms
   
        real, intent (inout), dimension (3)         :: a1vec
        real, intent (inout), dimension (3)         :: a2vec
        real, intent (inout), dimension (3)         :: a3vec
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
        integer iatom
        integer i
        integer j
        integer k
        integer index
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
        write (*,*) '  '
        write (*,*) ' Writing the density into a file: '

! write out density into density.xsf file (format of xcrysden visual code)
!  for details see www.xcrysden.org
   open ( unit = 301, file = 'density.xsf', status = 'unknown' )

! print the list of atoms
   if (icluster .eq. 1) then
 
    write (301,*) 'ATOMS'
    do iatom = 1,natoms
     write (301,'(i2,3f14.8)') nzx(imass(iatom)),(ratom2g(i,iatom),i=1,3)
    enddo

   else

    write (301,*) 'CRYSTAL'
    write (301,*) 'PRIMVEC'
    write (301,*) (a1vec(i),i=1,3)
    write (301,*) (a2vec(i),i=1,3)
    write (301,*) (a3vec(i),i=1,3)

    write (301,*) 'PRIMCOORD'
    write (301,*) natoms,1
    do iatom = 1,natoms
     write (301,'(i2,3f14.8)') nzx(imass(iatom)),(ratom2g(i,iatom),i=1,3)
    enddo

   endif
           
   write (301,*)
   write (301,*) 'BEGIN_BLOCK_DATAGRID_3D'
   write (301,*) 'density_3D'
   write (301,*) 'DATAGRID_3D_DENSITY'
   write (301,*) rm1, rm2, rm3
! print origin of the grid
   write (301,*) 0.0d0, 0.0d0, 0.0d0
! print lattice vector
   write (301,*) (a1vec(i),i=1,3)
   write (301,*) (a2vec(i),i=1,3)
   write (301,*) (a3vec(i),i=1,3)
   index = 0
! print values of the grid point
   do k = 0, rm3-1
     do j = 0, rm2-1
       do i = 0, rm1-1
          write (301,*) drhoG(index)
          index = index + 1
       enddo ! do i
     enddo ! do j
   enddo ! do k
   write (301,*) 'END_DATAGRID_3D'

   write (301,*) 'atomic_density_3D'
   write (301,*) 'DATAGRID_3D_DENSITY'
   write (301,*) rm1, rm2, rm3
! print origin of the grid
   write (301,*) 0.0d0, 0.0d0, 0.0d0
! print lattice vector
   write (301,*) (a1vec(i),i=1,3)
   write (301,*) (a2vec(i),i=1,3)
   write (301,*) (a3vec(i),i=1,3)
   index = 0
! print values of the grid point
   do k = 0, rm3-1
     do j = 0, rm2-1
       do i = 0, rm1-1
          write (301,*) rhoG0(index)
          index = index + 1
       enddo ! do i
     enddo ! do j
   enddo ! do k
   write (301,*) 'END_DATAGRID_3D'

   write (301,*) 'END_BLOCK_DATAGRID_3D'
! close file
   close (301)



! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
100     format (2x, i4)
101     format (2x, i6, 2x, f9.3)
102     format (2x, i2, 3(2x,f16.9), 3(2x,f16.9))
 
        close (unit = 63)
        return
      end subroutine writeout_den
