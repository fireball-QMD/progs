        module constants_fireball
!$ volatile cfac, delk, xlevi

! Bohr constant
         real, parameter :: abohr = 0.529177249
         real, parameter :: eq2 = 14.39975d0
         real, parameter :: Hartree = 14.39975d0/abohr

! Rydberg/eV conversion
         real, parameter :: ryd = 13.6d0

! Boltzmann's constant
real, parameter :: kb = 8.617343693082104d-5! 1.380662d-23 J/K * 6.24150974*10^18eV/J

! fovermp = f/mp = (1eV/1angstrom)/(mass proton) = angstrom/(fs**2).
! JPL 2000 I believe that this is incorrect after some lengthy discussion
! with Srini and Kurt.  The true conversion should be atomic mass units.
!        real, parameter :: fovermp = 0.0095790d0
         real, parameter :: fovermp = 0.009648957597d0

         real, parameter :: kconvert = 11604.49558d0
         real, parameter :: pi = 3.141592653589793238462643
         real, parameter :: spin = 2.0d0

! Gear algorithm constants_fireball
!         integer, parameter :: gear_order = 5  !  (only use 2-7)
!         real, dimension (0:gear_order) :: cfac
! JOM-info I change this so that I can run verlet (gear_order=2)
        integer, parameter :: gear_order = 2  !  (only use 2-7)
        real, dimension (0:5) :: cfac

! Kronecker delta
         real, dimension (3, 3) :: delk

! Levi-Civita parameter
         real, dimension (3, 3, 3) :: xlevi

! This is debatable, but they better be set high enough that the answers do not
! vary enough to concern you!  The forces are most susceptible to the norders.
! Higher order does not necessarily imply higher accuracy.  In fact too high
! of order can cause "ringing" (rapid occilations of the function).
! Negative numbers tell it to use a cubic spline of order abs(norder1)
! Splines provide better energy convservation and ovoid the ringing problem.
! Splines also avoid the problem with double numerical basis sets having
! unphysically low eigenvalues. 
! Note: what we used to call 6th order was really 5th!
! We generally use -5 which puts 3 points on each side of where you are
! interpolating.  Total number of points=abs(norder1)+1.
         integer, parameter :: norder1 = -5
! If we do one big 100+ point spline the we just build it once
! This gives continous f, f', f'', but not f'''.  All higher derivatives are 
! continous since they are zero.
!         logical, parameter :: superspline = .false.
         logical, parameter :: superspline = .true.

! What method do we use to interpolate 2-D grids
         integer, parameter :: D2intMeth = 1 ! 1 -> polynomials x then y
                                             ! All other methods sucked and were deleted
! twister data
        real amat (3, 3, 5) 
        logical haveDorbitals

        end module
