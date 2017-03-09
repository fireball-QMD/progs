! copyright info:
!
!                             @Copyright 1999
!                          Fireball2000 Committee
!
! ASU - Otto F. Sankey
!       Kevin Schmidt
!       Jian Jun Dong
!       John Tomfohr
!       Gary B. Adams
!
! Motorola - Alex A. Demkov
!
! University of Regensburg - Juergen Fritsch
!
! University of Utah - James P. Lewis
!                      Kurt R. Glaesemann
!
! Universidad de Madrid - Jose Ortega
!
!                  made possible by support from Motorola
!                         fireball@fermi.la.asu.edu

!
! fireball-qmd is a free (GPLv3) open project.

!
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


! sa.f
! Program Description
! ===========================================================================
! "true" simulated annealing. see S. Kirkpatrick, C.D. Gelatt, Jr., and
! M.P. Vecchi, Science 220, 671 (1983); N Metropolis, A.W. Rosenbluth, M.N.
! Rosenbluth, A.H. Teller, and E. Teller, J. Chem. Phys. 21, 1087 (1953).
!
! the procedure is a monte carlo procedure and not molecular dynamics.
! it finds a miminimu energy state by using forces and random displcements
! to head to low energy state. we put in a simulated annealing temperature, tsa,
! which allow a probabilty exp(-de/k*tsa) as defined below.
!
! the metropolis paper is the famous Monte Carlo paper, and the Kirkpatrick
! paper is THE paper on simulated annealing. I think Kirkpatrick coined the
! term, but i'm not sure.
!
! Program Declaration
! ===========================================================================
!
        subroutine sa(tsa,dxmax,kseed,natoms,Etotold,etotnew,ianneal,bsa,b)
        use dimensions
        use fragments
        implicit none

! Argument Declaration and Description
! ===========================================================================
        real tsa      ! simulated annealing temperature
        real dxmax    ! maximum allowed move per step in a direction
        integer kseed ! seed for the random number generato
        integer natoms! number of atoms
        real Etotold
        real etotnew
        integer ianneal(natoms) ! which atoms are annealed
        real bSA(3,natoms)   ! old positions on exit
        real b(3,natoms)     ! old positions, then new

! care must be taken to that dx is not too small (we never get anywhere),
! but also not too large (so that all "uphill" moves are basically
! forbidden).
!
! we allow a monte carlo move if -F dot dx is less than zero,
! and we allow with probability exp( - (-F dot dx)/k*tsa ) if -F dot dx
! is greater than zero.
!
! b(ix,n) the basis positions. these are also input.
! we replace the "old" value of b by b + (dx,dy,dz).

! Local Parameters and Data Declaration
! ===========================================================================

! Local Variable Declaration and Description
! ===========================================================================
        integer i
        integer in
        integer ix
        integer jx
        real dummy
        real dx,dy,dz
        real etherm

! Procedure
! ===========================================================================
        print*,'I know nothing about fixed fragments'
        stop ' i suck'

        write(*,*) '  '
        write(*,*) ' welcome .... simulated annealing (SA)'
        write(*,*) '  '
        write(*,*) ' etotold=',etotold
        write(*,*) ' etotnew=',etotnew
! ==============================================================
! Loop over all atoms in our basis.
! Did the energy go up or down from the last move.
        if(etotnew.lt.etotold)then
          write(*,*) ' Keep our position and reshuffle'
          goto 40
        end if
! Ooops. The energy went up. So lets go back to
! where we were.
        write(*,*) ' Go back to where we were and shuffle again.'
        do i=1,natoms
          do ix=1,3
            b(ix,i)=bSA(ix,i)
          end do
        end do
        etotnew=etotold
40      continue
! Now we prepare for the next SA step.
! First save where we are.
        do i=1,natoms
          do ix=1,3
            bSA(ix,i)=b(ix,i)
          end do
        end do
!
! ================================================
! Next move the atoms.
! etherm is the thermal energy in eV.
        etherm=(1./40.)*(tsa/300.)
        do 100 in=1,natoms
! We only anneal those if ianneal(i)=1
          if(ianneal(in).eq.0)go to 100
          write(*,*) ' '
! now move atoms dx,dy,dz by random # generator.
          call randomize(kseed,dummy)
          dx=(2*dummy-1)*dxmax
          call randomize(kseed,dummy)
          dy=(2*dummy-1)*dxmax
          call randomize(kseed,dummy)
          dz=(2*dummy-1)*dxmax
          write(*,*) ' i=',in,' dx,dy,dz=',dx,dy,dz
          write(*,*) ' Atom ',in
          write(*,89)(b(jx,in),jx=1,3)
 
! now guess the energy change if we allow this step.
          de=-(ftot(1,in)*dx+ftot(2,in)*dy+ftot(3,in)*dz)
! make the move if de<0
          if(de.ge.0.0d0)then
! make the move with probablity exp(-e/kT) is e>0.
            prob=dexp(-de/etherm)
            call randomize(kseed,proran)
            if(proran.gt.prob)go to 100
          end if
 
! if we end up here, we make the move. b ---> b+(dx,dy,dz)
          b(1,in)=b(1,in)+dx
          b(2,in)=b(2,in)+dy
          b(3,in)=b(3,in)+dz
          write(*,88)(b(ix,in),ix=1,3),(bSA(jx,in),jx=1,3)
100     continue
!
        write(*,*) ' completed SA simulated annealing.............'

! Format Statements
! ===========================================================================
89        format(' Initial b=',3f10.5)
88        format(' Final   b=',3f10.5,' bSA=',3f10.5)

        return
        end
