! Fortran 95/90 Module including procedure
module cproc
  interface
     subroutine cclient1()
       !DEC$ ATTRIBUTES C :: cclient1
     end subroutine cclient1     
     subroutine cclient2()
       !DEC$ ATTRIBUTES C :: cclient2
     end subroutine cclient2
     subroutine cclient()!marker)
       !DEC$ ATTRIBUTES C :: cclient
       ! !DEC$ ATTRIBUTES REFERENCE :: marker
       ! character marker
     end subroutine cclient
  end interface
end module cproc
