        module barrier

!$ volatile ibarrier,bar_how_often,barrier_push,barrier_tol,bar_stop_push,bar_too_much,bar_sav
!$ volatile ratom_final,barrier_achieved

! ---------------------------------------------------------------------------
! Information for energy barrier calculations. 
! ---------------------------------------------------------------------------
! Are even doing a barrier calculation
        integer :: ibarrier
! how_often defines how often to quench velocites.  0 means never.
        integer :: bar_how_often

! barrier_push is how much to push
        real :: barrier_push
! barrier_tol determines if we have reached the end
        real :: barrier_tol
! If stop_push=0 then push if dotp<0.  If stop_push>0 then for small values
! of dotp give a push also (0...1).  If you set to a very large number, say
! about 10000, and you set sav to 0.0, then you have only the push, and 
! the real forces are ignored.  Values like 1.6 (with sav=0.3 are useful to
! force it against its will, but still allow some flexability.
        real :: bar_stop_push    ! 0 to 1
! too_much defines when to slow down a good force, if it's too big and 
! will overshoot the goal.  This is a damping effect.
        real :: bar_too_much     ! A big number
! If sav=1 then just give a push.  If sav<1 then partially quench "bad" forces
        real :: bar_sav          ! 0 to 1

! The final geometry
        real, dimension (:, :), allocatable :: ratom_final

! Are we there yet?
        logical :: barrier_achieved

        end module
