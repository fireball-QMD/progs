 module md
! This module defines variables controlling MD  
! ===========================================================================

! ---------------------------------------------------------------------------
! Information for MD - equations of motion. 
! ---------------------------------------------------------------------------
        real dt
        real T_average
        real T_initial
        real T_final
        real T_instantaneous
        real T_previous
        real T_want
        real T_wantPrev
        real time
        real tkinetic
        real T_fermi
        real endtempdifference
        real T_increment
        
        real taurelax
        
        integer nstepf
        integer nstepi

        integer itime_step_g
! input/output files
        character (len = 30) acfile
        character (len = 30) xvfile

! constraints --------------------------------------------------------------

        logical :: fixCenOfMass = .false.
        real, dimension (3) :: rcmOld, rcmNew, rcmDiff
        real xmassTot, rcmDiffMag

 
 end module md
