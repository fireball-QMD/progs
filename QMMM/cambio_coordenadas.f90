subroutine cambio_coordenadas (a) 


use qmmm_module, only : mm_struct


implicit none

integer :: iatom, h, k
integer, intent (in) :: a

allocate ( mm_struct%totcoords_array(3,mm_struct%totatoms))


do h = 1, mm_struct%totatoms
   do k = 1, 3
      mm_struct%totcoords_array(k,h) = mm_struct%totcoords (h*3+k-3)  !cambio a matriz para sacarlo por pantalla
   end do
end do


        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' Atom Coordinates from mm_struct%totcoords before change: '
        write (*,200)
        write (*,201)
        write (*,200)
        do iatom = 1, mm_struct%totatoms
         write (*,202)  mm_struct%totcoords_array(:,iatom)
        end do
        write (*,200)
        write (*,*) '  '



mm_struct%totcoords=mm_struct%totcoords*(a+1) ! cambio de coordenadas a todos los atomos de la dinamica
 

do h = 1, mm_struct%totatoms
   do k = 1, 3
      mm_struct%totcoords_array(k,h) = mm_struct%totcoords (h*3+k-3)  !cambio a matriz para sacarlo por pantalla
   end do
end do



! Now write out the basis file information.
        write (*,*) '  '
        write (*,*) '  '
        write (*,*) ' Atom Coordinates from mm_struct%totcoords after change: '
        write (*,200)
        write (*,201)
        write (*,200)
        do iatom = 1, mm_struct%totatoms
         write (*,202)  mm_struct%totcoords_array(:,iatom)
        end do
        write (*,200)
        write (*,*) '  '
 
! Format Statements
! ===========================================================================
200     format (2x, 70('='))
201     format (8x,' x ', 8x, ' y ', 8x, ' z ')
202     format (2x,3(2x,f9.3))



end subroutine cambio_coordenadas


