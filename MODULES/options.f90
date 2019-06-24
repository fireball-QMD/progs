 module options
! This module provide information about basic options of the Fireball task
! ===========================================================================

! ---------------------------------------------------------------------------
! Toggles read from the options.input file.
! ---------------------------------------------------------------------------
! Level of theory (0=> harris, 1=> idogs, 2=> extended hubbard)
        integer itheory
! Level of exchange-correlation theory
! 0 => Horsfield,
! 1 => Generalized Sankey-Niklewski (GSN)
! 2 => McWEDA
        integer itheory_xc

        integer iharris
        integer idogs
        integer ihubbard
        integer imcweda
        integer igsn
        integer ihorsfield
        integer iks
        integer itdse
! JOM-add
        integer imdet

        integer icluster               ! Calculate gas-phase
        integer ifixcharge             ! Fix the charges
        integer ifixneigh              ! Fix the neighbors
        integer igauss                 ! Use gaussian approximation to 3-center
        integer iimage                 ! How often we image the box.
        integer ipathintegral          ! Add quantum effects - path integral
        integer iqout                  ! Charges to use (Lowdin or Mulliken)
        integer iquench                ! Quenching option
        integer ispin                  ! Spin interaction
        integer iensemble              ! Which ensemble?
        integer iordern                ! Perform linear-scaling algorithm
        integer iumbrella              ! Do umbrella sampling (on distances)
        integer ivdw                   ! Include van der Waals interactions
        integer idynmat                ! Dynamical matrix simulation
        integer iharmonic              ! whether to attach harmonic oscillators
        integer ithermoint             ! do thermodynamic integration
        integer ireducekpts            ! whether to reduce kpts by atm symmetry
        integer iendtemp               ! toggles T_final for MD calculations
        integer ineb                   ! do Nudged Elastic Band Method
        integer itrans                 ! do transport calculations
        integer igrid                  ! the grid projection
        integer ibias                  ! to apply the bias voltage
        integer igap                   ! to introduce LDA gap corrections GAP ENRIQUE-FF
        integer icDFT                  ! do the constrain DFT
        integer iProjWF                ! do projection within MDET simulations 
        integer iqmmm                  ! QM/MM Electrostatic embedding 
        integer mix_embedding          ! mix electrostatic and mechanical embedding
	integer iephc                  ! do e-ph coupling
        integer idftd3                 ! DFTD3 corrections
        integer imbd                   ! MBD corrections
        integer idipole                ! Long range term with XYZ dipole
        integer isocket                ! socket for i-pi
        integer ixczw                  ! Second order XC theory
! ---------------------------------------------------------------------------
! Other controls
! ---------------------------------------------------------------------------
        integer iforce
        integer imix
!CHROM
	integer iclassicMD
!END CHROM
        integer itestrange
        integer iconstraints (4)             ! big (4) constraints
        integer ioff2c (1:24)              ! for diagnostic purposes
        integer ioff3c (1:4)
        real testrange
        real  cut_embedding          ! cutoff for the mix embedding
        real  dftd3_s6, dftd3_rs6, dftd3_s18, dftd3_rs18, dftd3_alp
        character (len = 40) dftd3_func
        integer dftd3_version
        logical dftd3_tz 
        real dftd3_params (5)
        integer verbosity ! 0 minimal, 10 max         
        integer ntpr ! Every ntpr steps, information will be printed , if = 0 nothing will be write until the end        
        integer restartxyz ! if = 1, start the simulation from restart.xyz
        integer inputxyz ! if = 1, the coordinate input file is a .xyz

 end module options
