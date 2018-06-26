! copyright info:
!
!                             @Copyright 2005
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad de Madrid - Jose Ortega
! Universidad de Madrid - Pavel Jelinek
!
! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Motorola - Jun Wang
! Ohio State University - Dave Drabold

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


! cgo.f90
! Program Description
! ===========================================================================
!
!     FUNCTION CGO
!     Conjugate Gradient Optimization      
!     Description:
!     The optimization of geometry by Conj. Gradient Method.   
!     Used quadratic interpolation from three values of tot. energy
!     to find minima in given search direction
!     
!     input file: cgtol.input (placed in the actual directory)
!     input:
!            val ..  the total energy at the position
!            force .. the vector of forces
!            flag .. flag of the cycle
!            natoms .. # of atoms in system
!
!     output:
!            cgo .. the flag indicating the next action
!              return value: 10 .. convergence is reached
!                            1 ..  next scf step with forces
!                            2 ..  next scf step without forces
!
! ===========================================================================
! Code written by:
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
 subroutine cgo (val, force, iqout, xvfile, iforce, itheory, iwrtxyz)

   use optimization
   use interactions
   use configuration
   use charges

   implicit none
   
! Argument Declaration and Description
! ===========================================================================
   integer, intent (in)                             :: iqout
   integer, intent (in)                             :: itheory
   character (len=30),intent(in)                    :: xvfile
   integer, intent(in)                              :: iwrtxyz

   integer, intent (inout)                          :: iforce
   real, intent(inout)                              :: val
   real, dimension(3,natoms), intent(inout)         :: force

! Local Parameters and Data Declaration
! ===========================================================================
 
! Local Variable Declaration and Description
! ===========================================================================
   integer      :: i,ix
   real         :: Acoeff,Bcoeff,Ccoeff,beta
   real         :: Fi_max, res_etot
   real         :: dx,fx,fy,fz
   real         :: gg0,gg1,sqrh,orth,frs
   real         :: Fproj

   real, parameter   :: cg_eps2 = 1.0d-8
   real, parameter   :: cg_eps3 = 1.0d-8

! Allocate Arrays
! ===========================================================================
 
! Procedure
! ===========================================================================
! Auxilliary parameters

   select  case (istatus)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     ISTATUS = 0
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
   case (0)
!     we're going to generate the search direction
!     and to save actual the position, the energy and the forces
!     remeber that force_i = -dE/dx_i
!     we return the new position in the search direction
!     for localization of the minima in the given direction

      minflag = 0
      cg_iter = 0

      if(wrtout) then 
         write(*,201) cg_iter, val 
         write(*,*) ' ++ initial position :'
         do i = 1,natoms
            if (ishiftO .eq. 1) then
             write(*,36) i,ratom(:,i) - shifter(:),force(:,i)*mask(:,i)
            else
             write(*,36) i,ratom(:,i) ,force(:,i)*mask(:,i)
            endif
         enddo
      endif
! write temporal configuration into xv-file 
!      call writeout_xv (natoms, nspecies, xvfile, cg_iter, 0.0,    &
!     &                      imass, nzx, force)

! We save energy of actual configuration, we need it to evaluate quadratic 
! interpolation later. 
      etot0 = val
! set search vector
      do i = 1,natoms
         do ix = 1,3
!            g(ix,i) = force(ix,i)*mask(ix,i)
            g(ix,i) = force(ix,i)
            h(ix,i) = g(ix,i)
! save the actual position
            x0(ix,i) = ratom(ix,i)
            f0(ix,i) = force(ix,i)
         enddo
      enddo
! seve scf-charges
      if(itheory .eq. 1) then 
         do i = 1,natoms
            Q0(:,i) = Qin(:,i) 
         enddo
      endif


! Built the max. lenght of the next step in the search direction.
      dx = 0.0d0
      do i = 1,natoms
       do ix = 1,3
         fx = sqrt (h(1,i)**2.0d0 + h(2,i)**2.0d0 + h(3,i)**2.0d0)*mask(ix,i)
         dx = max (dx, fx)
!         do ix = 1,3 
!            if(abs(h(ix,i)) .gt. dx) dx = abs(h(ix,i))
!         enddo
       enddo
      enddo
! Fit maximal lenght of the step
      if(dx .lt. cg_drmax) then 
         alpha_cg = 1.0d0
      else
         alpha_cg = cg_drmax / dx
      endif

      if(wrtout) then 
         write(*,*) ' ++ CG parameters :'
         write(*,*) ' ++========================++'
         write(*,*) '  ++ alpha = ',alpha_cg, '  alpha*dx  = ', alpha_cg*dx
         write(*,*) '  ++ dx    = ',dx, '  drmax = ',cg_drmax
      endif

! Write new trial position
      do i = 1,natoms
       do ix = 1,3
         ratom(ix,i) = x0(ix,i) + alpha_cg*h(ix,i)*mask(ix,i)
       enddo
      enddo
      if(wrtout) then 
         write(*,*) ' ++==============================++'
         write(*,*) ' ++ Proposed trial position : '
         do i = 1,natoms
            if (ishiftO .eq. 1) then
              write( *,37) i, ratom(:,i) - shifter(:)
            else
              write( *,37) i, ratom(:,i)
            endif
         end do
      endif

! Set status flags 
      istatus = 2
      initer  = 1
      iforce  = 1
      
!++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     ISTATUS=1
!++++++++++++++++++++++++++++++++++++++++++++++++++++++
   case (1)
!     we're going to check criteria of convergence
!     if there are not satisfied we evaluate the new search direction
!     and we save the energy and the forces at the point
!     finally, we generate max. length of the next step and we evaluate
!     the new point according to the proposed search direction

      frs = 0.d0
      fx = 0.0d0 
      fy = 0.0d0
      fz = 0.0d0
      Fi_max = 0.0d0 
      do i = 1,natoms
         do ix = 1,3
            frs = frs + (force(ix,i)*mask(ix,i))**2.
            if (Fi_max .lt. abs(force(ix,i)*mask(ix,i))) then 
              Fi_max = abs(force(ix,i))
            endif
         enddo
         fx = fx + force(1,i)
         fy = fy + force(2,i)
         fz = fz + force(3,i)
      enddo

      frs = sqrt(frs)
      fx = fx**2.0d0 + fy**2.0d0 + fz**2.0d0
      fx = sqrt(fx) / natoms
      
!     dump information
      write(*,200) initer, val, Fi_max
      if(wrtout) then
         write(*,*) '++ ========================================='
         write(*,*) '++ Trial position:'
         write(*,*) '++ Etot= ', val,'[eV]'
         do i = 1,natoms
            if (ishiftO .eq. 1) then
              write(*,36) i, ratom(:,i) - shifter(:), force(:,i) * mask(:,i)
            else
              write(*,36) i, ratom(:,i), force(:,i) * mask(:,i)
            endif
         enddo
      endif

!-----------------------------------------------------------------------
!          CHECK UP HILL DISPLACEMENT 
!-----------------------------------------------------------------------     
      if(etot0 .lt. val) then
         
         initer = initer + 1
         write(*,*) '++ WARNNING: UP HILL DISPLACEMENT'
!         write(*,*) '++ I am going decrease step and repeat calc'
!     test to max. inner iter. 
         if(initer .gt. cg_minint) then

!     next step will be free of this resctriction
            if(minflag .eq. 1) then
               write(*,*) '+++++++++++++++++++++++++++++++++++++  '
               write(*,*) '+++    LOCAL MINIMUM REACHED      +++  '
               write(*,*) '+++++++++++++++++++++++++++++++++++++  '
! reconstruct the last fixed position and the charges
               write(*,*)'ETOT= ', etot0*natoms
               val = etot0
               do i = 1,natoms
                  ratom(:,i) = x0(:,i)
                  force(:,i) = f0(:,i)
               enddo
! restore minimal scf-charges
               if(itheory .eq. 1) then 
                  do i = 1,natoms
                     Qin(:,i) = Q0(:,i) 
                  enddo
               endif
               if (icg2md .eq. 1) then 
! set terminate flag (go to MD)
                istatus = 20
                write (*,*) '+++ the minimization switch to quenching MD +++ '
               else
                istatus = 10
               endif
               return

            endif ! if(minflag)

!     the next exceed of inter. cycle will cause stop of program
            minflag = 1
            write(*,*) '++ number of inner iter. exceeded'
            write(*,*) '++ I am going to reset the search dir.'
            do i = 1,natoms
               h(:,i) = f0(:,i)
            enddo
           
           
!     built up the max.lenght of the next step
!     alpha = cg_drmax/h2
!            cg_drmax = cg_drmax / 2.0d0
            dx = 0.0d0
            do i = 1,natoms
             do ix = 1,3
               fx = sqrt(h(1,i)**2.0d0 + h(2,i)**2.0d0 + h(3,i)**2.0d0)*mask(ix,i)
               dx = max (dx, fx)
!               do ix = 1,3 
!                  if(abs(h(ix,i)) .gt. dx) dx = abs(h(ix,i))
!               enddo
             enddo !do ix
            enddo
            if(dx .lt. cg_drmax) then 
               alpha_cg = 1.0d0
            else
               alpha_cg = cg_drmax / dx
            endif

            if(wrtout) then
               write(*,*) '++ alpha= ',alpha_cg,' alpha*dx= ', alpha_cg*dx
               write(*,*) '++ dx=    ',dx,' ++ drmax= ',cg_drmax
            endif
!     write trial position
            do i = 1,natoms
             do ix = 1,3
               ratom(ix,i) = x0(ix,i) + alpha_cg*h(ix,i)*mask(ix,i)
             enddo ! do ix
            enddo
            
!     set gamma to one for future
            gamma = 1.0d0
!     return point
            istatus = 2
            initer = 1
            return

         endif ! if(initer)
! -----------------------------------------------------------------------
! Trial tot. energy is higher than fixed tot. energy of fixed configuration 
! Let's decrease the lenght of the previous step by 'dummy' and generate 
! new position


         gamma = gamma * cg_dummy
         do i = 1,natoms
          do ix = 1,3
            ratom(ix,i) = x0(ix,i) + gamma*alpha_cg*h(ix,i)*mask(ix,i)
          enddo ! do ix 
         enddo ! do i

         if(wrtout) then 
            write(*,*) '++ gamma= ',gamma
            write(*,*) '++ New trial (reduced) position  '
            write(*,*) '++ -------------------------------'
            do i = 1,natoms
             if (ishiftO .eq. 1) then
               write( *,37) i, ratom(:,i) - shifter(:)
             else
               write( *,37) i, ratom(:,i)
             endif
            enddo
         endif

!     go back and solve scf again with new set of coord.
         istatus = 1
         return

!     end of 'if(etot0.lt.val)'
      endif

!     =======================================================
!     TRIAL POSITION IS ACCEPTED
!     =======================================================
!     increment number of cg iteration
      cg_iter = cg_iter + 1
!     dump the header
      if(wrtout) then 
         write(*,*) '++ -------------------------------------'
         write(*,*) '++    TRIAL POSITION ACCEPTED         ++'
         write(*,*) '++ -------------------------------------'
      endif
      write(*,201) cg_iter, val
      
!     now we save on disk the actual minimized configuration
! writeout GEOM
      open (unit= 86, file='answer.bas', status='unknown')
      write (86,*) natoms
      do i=1,natoms
        if (ishiftO .eq. 1) then
          write (86,700) nzx(imass(i)), ratom(:,i) - shifter(:)
        else
          write (86,700) nzx(imass(i)), ratom(:,i)
        endif
      end do
      close(unit=86)

! writeout XYZ
      if (iwrtxyz .eq. 1)  then
       open (unit= 87, file='answer.xyz', status='unknown', position = 'append')
       write (87,*) natoms
       write (87,*) '## CG step no.',cg_iter
       do i=1,natoms
        if (ishiftO .eq. 1) then
          write (87,800) symbol(i), ratom(:,i) - shifter(:)
         else
          write (87,800) symbol(i), ratom(:,i) 
         endif
       end do
       close(unit=87)
      endif

! writeout CHARGES
      open (unit= 87, file='CHARGES', status='unknown')
      write (87,600) natoms,basisfile,iqout
      do i = 1,natoms
         write(87,601) (Qin(ix,i),ix=1,nssh(imass(i)))
      end do
      close(87)
      
!     write fixed postion in output file

      if(wrtout) then 
         write (*,*) '++ ======================================'
         write (*,*) '++ The fixed position:'
         write (*,*) '++ Etot[eV] =', val
         do i = 1,natoms
          if (ishiftO .eq. 1) then
            write (*,37) i,ratom(:,i) - shifter(:)
           else
            write (*,37) i,ratom(:,i) 
           endif
         enddo
      endif
!      call writeout_xv (natoms, nspecies, xvfile, cg_iter, 0.0, imass, & 
!     &                  nzx, force)           
!     =========================================================

!     refresh flag 
      minflag = 0
!     =========================================================
!     check the criteria convergence
      res_etot = abs(etot0 - val)
      write (*,*) ' ++ CONVERGENCE STATUS :  '
      write (*,*) ' ++ ============================================'
! etot
      if (res_etot .le. energy_tol) then
         write (*,100) res_etot, energy_tol
      else
         write (*,101) res_etot, energy_tol 
      endif
! ftot
      if (Fi_max .le. force_tol) then
         write (*,110) Fi_max, force_tol
      else
         write (*,111) Fi_max, force_tol 
      endif

! Test of convergence
      if (Fi_max .le. force_tol .and. res_etot .le. energy_tol) then
         istatus = 10
         return
      endif

! Check number of CG iterations
      if (cg_iter .eq. cg_maxstep) then 
          istatus = 10
          return 
      endif


!     =========================================================

!     evaluate the search path correction
!     ----------------------------------------------------------
      gg0=0.0
      gg1=0.0
      do i = 1,natoms
         do ix = 1,3
            gg0 = gg0 + mask(ix,i)*g(ix,i)**2
!     Fletcher-Reeves version
!            gg1=gg1+(force(ix,i)*mask(ix,i))**2
            gg1 = gg1 + mask(ix,i)*force(ix,i)**2
!     Polak-Ribiere version
!            gg1=gg1+(force(ix,i)*mask(ix,i)-g(ix,i))*force(ix,i)*mask(ix,i)
!            gg1 = gg1 + (force(ix,i) - g(ix,i)) * force(ix,i)
         enddo
      enddo
      beta = gg1 / gg0
!     ----------------------------------------------------------------

!     update vectors for new itertive cycle
      orth = 0.0d0
      sqrh = 0.0d0
      frs  = 0.0d0
      do i = 1,natoms
         do ix = 1,3
!            g(ix,i) = force(ix,i)*mask(ix,i)
            g(ix,i) = force(ix,i)
            h(ix,i) = g(ix,i) + beta*h(ix,i)
            orth = orth + h(ix,i)*g(ix,i)
            sqrh = sqrh + h(ix,i)**2.0d0
            frs = frs + g(ix,i)**2.0d0 
         enddo
      enddo
      
      sqrh = sqrt(sqrh)
!     test of orthogonality
!     -----------------------------------------------------------
      orth = abs(orth / gg1)
      if(wrtout) then 
         write(*,*) '++ gg0= ',gg0
         write(*,*) '++ gg1= ',gg1
         write(*,*) '++ orth = ',orth
         write(*,*) '++ beta = ',beta
         write(*,*) '++ sqrh = ',sqrh
         write(*,*) '++ frs  = ',frs
      endif
      if(orth .lt. cg_eps3) then
         write(*,*) '++ warninng generated direction is nearly' 
         write(*,*) '++ orthogonal to the steep descent'
         write(*,*) '++ algorithm will be restarted'
         do i = 1,natoms
            do ix = 1,3
               h(ix,i) = g(ix,i)
            enddo
         enddo
      endif
!     ------------------------------------------------------------

!     test if the derivation is small than required precision
!     very accurate condition, but cause more outer iteration
      if(sqrh.lt.cg_eps2) then
         write(*,*) '++ Minimum is reached !!'
         write(*,*) '++ sqrh is less than cg_eps2'
         istatus = 10
         return
      endif

!     save the actual configuration
!     save the energy
      etot0 = val
      do i =1,natoms
         x0(:,i) = ratom(:,i)
         f0(:,i) = force(:,i)
      enddo
! save scf charges
      if(itheory .eq. 1) then 
         do i = 1,natoms
            Qin(:,i) = Q0(:,i) 
         enddo
      endif
           
!     built up the max.lenght of the next step
      dx = 0.0d0
      do i = 1,natoms
        do ix = 1,3
         fx = sqrt (h(1,i)**2.0d0 + h(2,i)**2.0d0 + h(3,i)**2.0d0)*mask(ix,i)
         dx = max (dx, fx)
!         do ix = 1,3 
!            if(abs(h(ix,i)).gt.dx) dx = abs(h(ix,i))
!         enddo
       enddo ! do ix
      enddo
      if (dx .lt. cg_drmax) then
         alpha_cg = 1.0d0
      else
         alpha_cg = cg_drmax / dx
      endif
           
      if(wrtout) then 
         write(*,*) '++ alpha= ',alpha_cg, 'alpha*dx= ', alpha_cg*dx
         write(*,*) '++ dx=   ',dx
         write(*,*) '++ drmax= ',cg_drmax
      endif

!     write trial position
      do i = 1,natoms
       do ix = 1,3
         ratom(ix,i) = x0(ix,i) + alpha_cg*h(ix,i)*mask(ix,i)
       enddo ! do ix
      enddo

!     return point
      istatus = 2
      initer = 1

!++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     ISTATUS = 2
!++++++++++++++++++++++++++++++++++++++++++++++++++++++
   case(2) 
!     we're going to search minimum in the direction
!     and find the new position for the next step
!     also we check the max. lenght of the step
!     we return in icgforce 
!----------------------------------------------------
!     quadratic interpolation is used to reach minimum
!     the function in search direction is approximated by
!     projection N-D --> 1-D; xi --> gamma
!
!     Version used two Etot and one force
!     F(gamma)= A*gamma^2 + B*gamma + C :
!     const. C = etot0
!            B = -f1*dr12 
!            A = e2 - B - C  
  
      etot1 = val
      if(wrtout) then
         write(*,*) ' ++ ========================================='
         write(*,*) ' ++ Etot= ', etot1,'[eV]'
         do i = 1,natoms
           if (ishiftO .eq. 1) then
             write(*,37) i,ratom(:,i) - shifter(:)
           else
             write(*,37) i,ratom(:,i) 
            endif
         enddo
      endif

      
!c--------------------------------------------------------
!c ver.2. 2xEtot & 1Force
!c you nedd also avoid to calculate next Etot point (comment 
!c above if-body), correct charge mixing (simple mixing between 
!c two points) and define 'Fproj' variable
!c      
! evaluate distance between two points
      Fproj = 0.
      do i = 1,natoms
         do ix = 1,3
            Fproj = Fproj - alpha_cg*h(ix,i)*g(ix,i)*mask(ix,i)
         enddo
      enddo

! evaluate C coeff.
      Ccoeff = etot0
! evaluate B coeff.
      Bcoeff = Fproj
! evaluate A coeff.; e2 ~ val
      Acoeff = val - Bcoeff - Ccoeff

      if(wrtout) then
         write (*,*) ' ++ CG parameters of quadratic interpolation :'
         write (*,*) ' ++ ==========================================='
         write (*,*) ' ++ etot0  = ',etot0
         write (*,*) ' ++ etot1  = ',etot1
         write (*,*) ' ++ Fproj  = ',Fproj
         write (*,*) ' ++ alpha  = ',alpha_cg
         write (*,*) ' ++ A      = ',Acoeff
         write (*,*) ' ++ B      = ',Bcoeff
         write(*,*) ' ++ C      = ',Ccoeff
      endif

! Let's check the concave shape of the interpolation function
      if (Acoeff .le. 0.0d0) then 

         write(*,*) '++ Warninng: concave shape of the fitting parabola'
         gamma = 1.0d0
!     if we go uphill then we refine the step
         if (etot1 .gt. etot0) then
            if (val .gt. etot0) gamma = gamma*0.5*cg_dummy
         endif

      else 
         
! find minimum: dF(x2)=0
         gamma = -1.0d0 * Bcoeff / (2.0d0 * Acoeff)
! we bound the maximal displacement to const*cg_drmax 
         if (gamma .gt. 2.0d0) gamma = 2.0d0
         if (wrtout) write(*,*) '++ gamma_min= ',gamma
      endif

! Get new postion based on quadratic guess  
      do i = 1,natoms
       do ix = 1,3
         ratom(ix,i) = x0(ix,i) + gamma*alpha_cg*h(ix,i)*mask(ix,i)
       enddo ! do ix
      enddo ! do iatom


! Set auxilliary variables for next cycle
      istatus = 1
      iforce  = 1
   
   end select

   return

! Deallocate Arrays
! ===========================================================================
 
! Format Statements
! ===========================================================================
36 format(' ++ atom= ',i4,6f12.5) 
37 format(' ++ atom= ',i4,3f12.5)
38 format(i3,3f18.9)
100 format (2x,'  +++ Etot  RES =', f16.8,'  TOL = ',f16.8,'   CONVERGED ')
101 format (2x,'  +++ Etot  RES =', f16.8,'  TOL = ',f16.8,'   NOT CONVERGED ')
110 format (2x,'  +++ Fmax  RES =', f16.8,'  TOL = ',f16.8,'   CONVERGED ')
111 format (2x,'  +++ Fmax  RES =', f16.8,'  TOL = ',f16.8,'   NOT CONVERGED ')
200 format (2x,' ++++ initer= ',i3,' Etot= ',f16.8,' Fi_max= ',f14.6)
201 format (2x,' ++++ iter. no.',i3,' etot = ',f16.8 )
600 format (2x, i5, 2x, a40, 2x, i2)
601 format (2x, 10f14.8)
700 format (2x, i2, 3(2x,f18.8))
800 format (2x, a2, 3(2x,f12.6))
968 format(3f16.8,i2,2f10.6,2f4.1,i2,f12.1)
   

   return
        
 end subroutine cgo


