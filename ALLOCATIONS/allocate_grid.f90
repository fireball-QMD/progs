! copyright info:
!
!                             @Copyright 2001
!                           Fireball Committee
! University of Utah - James P. Lewis, Chair
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
! Ohio State University - Dave Drabold
 
!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
! allocate_rho.f90
! Program Description
! ===========================================================================
!       This routine allocates the arrays which store information about 
! wavefunctions, neutral atomic potentials (vna) of each specie and finally
! the arrays of numerical grid 
!
! ===========================================================================
! Code written by:

! ===========================================================================
!
! Program Declaration
! ===========================================================================
  subroutine allocate_grid (natoms, nspecies)

    use wavefunction
    use vnneutral
    use interactions
    use grid
    implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input

    integer, intent (in) :: natoms
    integer, intent (in) :: nspecies
 
! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
 
! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================

! wave function part
    allocate (mesh_wf (nsh_max, nspecies))
    allocate (drr_wf (nsh_max, nspecies))
    allocate (rmax_wf (nsh_max, nspecies))
    allocate (rr_wf (wfmax_points, nsh_max, nspecies))
    allocate (wf_spline (wfmax_points, nsh_max, nspecies))
    allocate (wf (wfmax_points, nsh_max, nspecies))

! vneutral part
    allocate (mesh_na (nspecies))
    allocate (rmax_na (nspecies))
    allocate (drr_na (nspecies))
    allocate (rr_na (max_vna_points, nspecies))
    allocate (vnna_spline (max_vna_points, nspecies))
    allocate (vnna (max_vna_points, nspecies))

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
 
    return
  end subroutine allocate_grid
