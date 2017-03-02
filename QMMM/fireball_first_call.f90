subroutine fireball_first_call


use qmmm_module

implicit none



allocate ( qmmm_struct%Qresp(qmmm_struct%nquant_nlink))
allocate ( qmmm_struct%Qneutral_TOT(qmmm_struct%nquant_nlink))


! Read amber information

             call amber_fireball ()

! Read basic informations about the task

             call initbasics ()
        
! Read data
             call readdata ()

end subroutine fireball_first_call
