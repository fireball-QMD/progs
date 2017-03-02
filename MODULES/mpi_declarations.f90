        module mpi_declarations
! FIXME replace these volatile attributes with thread private when supported
!$ volatile mpi_whatever_real, mpi_whatever_double

! MPI real type indicator; depends on the default real type
         integer mpi_whatever_real
         integer mpi_whatever_double

         integer mpi_on
        end module
