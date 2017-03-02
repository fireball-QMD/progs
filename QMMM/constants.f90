! <compile=optimized>

!************************************************************************
!                              AMBER                                   **
!                                                                      **
!                        Copyright (c) 2009                            **
!                    Ross Walker and David A. Case                     **
!                       All Rights Reserved.                           **
!                                                                      **
!************************************************************************

!+ Specification and control of Amber's working precision


! Description:
! Preprocessor directives that characterize the floating-point
! working precision as single or double precision.
! The current scheme guarantees internal consistency at the expense
! of flexibility.  A need for flexibility has yet to appear.
! The preprocessor guard  prevents multiple, and thus
! inconsistent, definitions.
! The default working precision is double precision.
! User control of the working precision at build time should be
! exercised via the preprocessor name _REAL_.
! To build a single precision Amber use
!     make -e AMBERBUILDFLAGS=' -D_REAL_ '
! The preprocessor names that characterize the precision are

!   _REAL_     precision type specifier.
!              Use  _REAL_ foo  as the precision independent
!              notation for  double precision foo  and  real foo.

!   AMBER_MPI_REAL
!              MPI precision type specifier.
!              Use AMBER_MPI_REAL as the precision independent
!              notation for MPI_DOUBLE_PRECISION and MPI_REAL.

!   D_OR_S()   precision prefix for the BLAS and LAPACK Library routines.
!              Use, e.g.,  D_OR_S()axpy(...)  as the precision independent
!              notation for daxpy(...) and saxpy(...).

!   DPREC      defined when the working precision is double;
!              undefined when the working precision is single.

!   VD_OR_VS() precision prefix for the Intel Vector Math Library routines.
!              Use, e.g.,  VD_OR_VS()exp(...)  as the precision independent
!              notation for vdexp(...) and vsexp(...).






!+++++++++++++++++++++++++++++++++++++++
!This module contains various parameters
!and constants used by the different 
!routines that make up sander.
!
!If you want to use one of the constants
!in your routine you should include the
!line:
!
!use constants, only : xxx, yyy, zzz
!
!where xxx,yyy,zzz are the constants you plan
!to use in your routine.
!This line needs to go before the
!implicit none declaration.
!
! Written by: Ross Walker (TSRI, 2005)
!++++++++++++++++++++++++++++++++++++++++

module constants

!------------------------------------------------------------
! Generic Floating Point Constants
double precision, parameter :: TEN_TO_MINUS2  = 1.0d-2
double precision, parameter :: TEN_TO_MINUS3  = 1.0d-3
double precision, parameter :: TEN_TO_MINUS4  = 1.0d-4
double precision, parameter :: TEN_TO_MINUS5  = 1.0d-5
double precision, parameter :: TEN_TO_MINUS6  = 1.0d-6
double precision, parameter :: TEN_TO_MINUS8  = 1.0d-8
double precision, parameter :: TEN_TO_MINUS10 = 1.0d-10
double precision, parameter :: TEN_TO_PLUS3   = 1.0d+3
double precision, parameter :: TEN_TO_PLUS10  = 1.0d+10

double precision, parameter :: zero      = 0.0d0
double precision, parameter :: one       = 1.0d0
double precision, parameter :: two       = 2.0d0
double precision, parameter :: three     = 3.0d0
double precision, parameter :: four      = 4.0d0
double precision, parameter :: five      = 5.0d0
double precision, parameter :: six       = 6.0d0
double precision, parameter :: seven     = 7.0d0
double precision, parameter :: eight     = 8.0d0
double precision, parameter :: nine      = 9.0d0
double precision, parameter :: ten       = 10.0d0
double precision, parameter :: eleven    = 11.0d0
double precision, parameter :: twelve    = 12.0d0
double precision, parameter :: sixteen   = 16.0d0
double precision, parameter :: twenty    = 20.0d0
double precision, parameter :: thirtytwo = 32.0d0
double precision, parameter :: sixtyfour = 64.0d0

double precision, parameter :: half         = one/two
double precision, parameter :: third        = one/three
double precision, parameter :: fourth       = one/four
double precision, parameter :: fifth        = one/five
double precision, parameter :: sixth        = one/six
double precision, parameter :: seventh      = one/seven
double precision, parameter :: eighth       = one/eight
double precision, parameter :: ninth        = one/nine
double precision, parameter :: tenth        = one/ten
double precision, parameter :: eleventh     = one/eleven
double precision, parameter :: twelfth      = one/twelve
double precision, parameter :: sixteenth    = one/sixteen
double precision, parameter :: thirtysecond = one/thirtytwo
double precision, parameter :: sixtyfourth  = one/sixtyfour

double precision, parameter :: thirtieth    = one/30.0d0

!------------------------------------------------------------


!------------------------------------------------------------
! Physical Constants
double precision, parameter :: HBAR = 627.509d0 * 0.0241888d-3 * 20.455d0 !Planck's constant in internal units
double precision, parameter :: J_PER_CAL = 4.184d0
double precision, parameter :: JPKC = J_PER_CAL * 1000.0d0 !kilocalories per joule
double precision, parameter :: BOLTZMANN = 1.380658d-23 !Boltzmann's constant in J/K
double precision, parameter :: AVOGADRO = 6.0221367d+23 !Avogadro's number
double precision, parameter :: KB = (BOLTZMANN * AVOGADRO) / JPKC !Boltzmann's constant in internal units
double precision, parameter :: AMBER_ELECTROSTATIC = 18.2223d0
double precision, parameter :: AMBER_ELECTROSTATIC2 = AMBER_ELECTROSTATIC * AMBER_ELECTROSTATIC
!Ratio by which to scale amber charges to get electron charges - amberchg * oneqscale = electron charges
! = 1.0 / 18.2223d0
double precision, parameter :: INV_AMBER_ELECTROSTATIC = 1.0d0/AMBER_ELECTROSTATIC
double precision, parameter :: INV_AMBER_ELECTROSTATIC2 = 1.0d0/AMBER_ELECTROSTATIC2

double precision, parameter :: CHARGE_ON_ELEC = 1.60217733d-19 !Charge on an electron in Coulombs
double precision, parameter :: BOHRS_TO_A = 0.529177249D0   ! Bohrs * this = angstroms - Same constants as used in dynamo v2.
!double precision, parameter :: BOHRS_TO_A = 0.52917706D0   ! Bohrs * this = angstroms - Same constants as used in Gaussian 98
!double precision, parameter :: BOHRS_TO_A = 0.529177D0     ! Bohrs * this = angstroms - Same constants as used in Mopac6 hcore.f
!double precision, parameter :: BOHRS_TO_A = 0.529167D0     !                            as used in Mopac6 repp.f
double precision, parameter :: A_TO_BOHRS = 1.0d0 / BOHRS_TO_A
!double precision, parameter :: A_TO_BOHRS = 1.88976D0      !Same constants as used in Mopac6 gover.f
double precision, parameter :: A2_TO_BOHRS2 = A_TO_BOHRS * A_TO_BOHRS !Mopac6 uses 3.5711928576D0 in gover.f for this.
double precision, parameter :: A3_TO_BOHRS3 = A2_TO_BOHRS2 * A_TO_BOHRS
double precision, parameter :: A4_TO_BOHRS4 = A2_TO_BOHRS2 * A2_TO_BOHRS2

!double precision, parameter :: ONE_AU = 27.2113962d0 !One atomic unit of energy in eV.
!double precision, parameter :: AU_TO_EV = ONE_AU  !Conversion from AU to EV - not used because we match dynamo v2 below.
double precision, parameter :: AU_TO_EV = 27.21d0 !Conversion from AU to EV - Same as dynamo v2 uses and Gaussian 98
                                        !Note (RCW+MC): more precise would be: 1 a.u. 27.211396 eV
                                        !Mopac6 uses 27.21D0 in calpar.f, delri.f and repp.f but in
                                        !ffhpol.f it uses 27.2107 and in the manual it quotes 27.211
double precision, parameter :: HALF_AU_TO_EV = AU_TO_EV * half
double precision, parameter :: FOURTH_AU_TO_EV = AU_TO_EV * fourth
double precision, parameter :: EIGHTH_AU_TO_EV = AU_TO_EV * eighth
double precision, parameter :: SXNTH_AU_TO_EV = EIGHTH_AU_TO_EV*half
double precision, parameter :: A2_TO_BOHRS2xAU_TO_EV = A2_TO_BOHRS2*AU_TO_EV

!double precision, parameter :: EV_TO_KCAL = 23.060362D0      !Conversion from EV to KCAL/MOL
!Dynamo parameter
double precision, parameter :: EV_TO_KCAL = 23.061d0  !Dynamo's conversion
                                            !Mopac6 uses 23.061 in ffhpol.f analyt.f compfg.f datin.f dcart.f
                                            !                      delri1.f delri2.f deritr.f interp.f iter.f
                                            !                      moldat.f mopac.f
double precision, parameter :: KCAL_TO_EV = one / EV_TO_KCAL

double precision, parameter :: AU_TO_KCAL = AU_TO_EV*EV_TO_KCAL !1 hartree.

!------------------------------------------------------------
!Numeric Constants
double precision, parameter :: PI      = 3.1415926535897932384626433832795d0

!The BOOK says :
!
!2Chronicles 4:2 reads thus, 'Also he made a molten sea of ten cubits 
!from brim to brim, round in compass, and five cubits the height thereof; 
!and a line of thirty cubits did compass it round about.'
!
!Hence, Pi is exactly equal to three and there is nothing more to discuss!
!
!If you want to use the value of PI defined by 'the BOOK' then uncomment
!the following line and comment out the definition above...
!double precision, parameter :: PI = 3.0d0

double precision, parameter :: PI2     = PI*PI
double precision, parameter :: HALFPI = PI * 0.5d0
double precision, parameter :: TWOPI  = 2.0d0 * PI
double precision, parameter :: FOURPI = 4.0d0 * PI
double precision, parameter :: INVPI  = 1.0d0 / PI
double precision, parameter :: SQRTPI = 1.77245385090551602729816748334d0 !sqrt(PI)
double precision, parameter :: INVSQRTPI = 1.0d0 / SQRTPI
double precision, parameter :: DEG_TO_RAD = PI / 180.0d0
double precision, parameter :: RAD_TO_DEG = 180.0d0 / PI 
double precision, parameter :: LN_TO_LOG = 2.30258509299404568402d0  ! log(1.0d1)

double precision, parameter :: SQRT2     = 1.4142135623730950488016887242097d0
double precision, parameter :: INVSQRT2  = 1.0d0 / SQRT2

!------------------------------------------------------------
!Generalised Born Constants
double precision, parameter :: alpb_alpha = 0.571412d0 !Alpha prefactor for alpb_alpha

!------------------------------------------------------------
! Unusual Constants
integer, parameter :: RETIRED_INPUT_OPTION = -10301 ! first 5 digit palindromic prime
integer, parameter :: NO_INPUT_VALUE = 12344321  ! from Bob Duke


end module constants

