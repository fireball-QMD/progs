! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.

! charges_L2.f90
! Program Description
! ===========================================================================
!   This routine computes the occupation numbers by minimizing the average quadratic error
! ===========================================================================
! ===========================================================================
! Program Declaration
! ===========================================================================

                subroutine charges_L2(Kscf,igauss)


        use charges
        use configuration
        use density
        use forces
        use interactions
        use neighbor_map
        implicit none


! Argument Declaration and Description
! ===========================================================================
! Input
        integer, intent (in) :: Kscf
        integer, intent (in) :: igauss

! Local Parameters and Data Declaration
! ===========================================================================
! Define in the module interactions

! Local Variable Declaration and Description
! ===========================================================================
        integer ialp
        integer iatom
        integer ialpha
        integer ibeta
        integer imu
        integer in1
        integer in2
        integer in3
        integer indna
        integer isorp1
        integer isorp2
        integer ineigh
        integer interaction
        integer interaction0
        integer inu
        integer isorp
        integer issh
        integer jatom
        integer iforce
        integer matom
        integer mneigh

        integer info
        integer, dimension(nsp_tot+1) :: pivots         

        real cost
        real x
        real y
        real I

        real, dimension (3, 3, 3) :: deps
        real, dimension (3, 3) :: eps
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: rhat


        ! molecular spehrical matrices
        real, dimension (numorb_max, numorb_max) :: rhomx
        real, dimension (3, numorb_max, numorb_max) :: rhompx

        real, dimension (nsh_max, nsh_max) :: rhomm
        real, dimension (3, nsh_max, nsh_max) :: rhompm
        real, dimension (3) :: rnabc
        real, dimension (3) :: rna
        real, dimension (3) :: sighat

        real, dimension (nsp_tot+1,1) :: b
       !real, dimension (nsp_tot+1,1) :: bb
         
        real, dimension(nsp_tot,nsp_tot) :: J    
        real, dimension(nsp_tot+1,nsp_tot+1) :: Atemp   
        
    ! Procedure
! ===========================================================================

        ! Initialize matrices
    Qout=0.0d0

    !#######WRITE FOR DEBUGGING
          do iatom = 1,natoms
              in1 = imass(iatom)
              do issh = 1,nsp(in1)   !NEW: nsp here!!!
                 ialpha = issh + degelec_sp(iatom) !NEW: degelec_sp here!!!
                 if ( in1 .eq. 1 ) then
                   write(*,*)'1: The OUTPUT CHARGES for this H are:',Qout(1,iatom),Qout(2,iatom)
                 end if
                ! write(*,*) 'In charges_L2: issh,iatom,Qout', issh, iatom,
                ! Qout(issh,iatom)
              enddo
          enddo
     !#######WRITE FOR DEBUGGING



    if (Kscf .eq. 1) then 
     J=0.0d0 
     Aglob=0.0d0
    endif
   
     b=0.0d0

!===================================
!Computations start now:

        write(*,*) 'we are for the ',Kscf,'th time in the subroutine charges_L2'
!====================================
!====================================
!Next build the matrix of spherical two center integrals, J, with coefficients J(ialpha,ibeta)
       if (Kscf .eq. 1) then
        do iatom = 1,natoms
         in1 = imass(iatom)
         r1 = ratom(:,iatom)
          do ineigh = 1,neighn(iatom)
           jatom = neigh_j(ineigh,iatom)
           in2=imass(jatom)
           !======================================================================
    ! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
    ! This gives us the distance dbc (or y value in the 2D grid).
            r2 = ratom(:,jatom)
            r21(:) = r2(:) - r1(:)
            y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))

             ! Find the unit vector in sigma direction.
            if (y .lt. 1.0d-05) then
             sighat(1) = 0.0d0
             sighat(2) = 0.0d0
             sighat(3) = 1.0d0
           else
            sighat(:) = r21(:)/y
           end if
           call epsilon (r2, sighat, eps)
           call deps2cent (r1, r2, eps, deps)
            iforce = 0
            interaction0 = 22
            in3=in1


            do isorp2 = 1,nsp(in2)   !NEW:  nsp here!!!
             ibeta = isorp2+degelec_sp(jatom) !NEW: degelec_sp here !!!
   !===
    !Now we compute the (spherical) two center integral for Jisorp2,in1

 
            call doscentrosS (interaction0, isorp2, iforce, in1, in2, &
     &                                in3, y, eps, rhomm, rhompm)

          
    !===================================================================
    !End of computation of the two center integral
   !===
           do isorp1= 1,nsp(in1)  !NEW: nsp here!!!
            ialpha = isorp1+degelec_sp(iatom) !NEW: degelec_sp here!!!
            

            J(ialpha,ibeta)=rhomm(isorp1,isorp1)  
 
          enddo !end loop over isorp1
         enddo !end loop over isorp2
        enddo !end loop over ineigh
       enddo !end loop over iatom   
      endif !end if Kscf = 1
       

!===================================
!===================================
!Next build the matrix of coefficients, A, for the linear system
!NOTE: theres gotta be a better way to achieve this
       if (Kscf .eq. 1) then
        do ialpha = 1,nsp_tot  !NEW
         do ibeta = 1,nsp_tot  !NEW: nsp here !!!
           Aglob(ialpha,ibeta)=J(ialpha,ibeta)
         enddo
        enddo

      

        Aglob(nsp_tot+1,nsp_tot+1)=0 !NEW: nsp_tot here!!!

       do ialpha = 1,nsp_tot !NEW: nsp_tot here!!!
        Aglob(ialpha,nsp_tot+1) = 1 !NEW: nsp here!!!
       enddo

       do ibeta = 1,nsp_tot !
        Aglob(nsp_tot+1,ibeta) = 1 !NEW: nsp here!!!
       enddo  
     endif !end if Kscf = 1

    ! write(*,*) 'Construction of A succesful'
    ! write(*,*) 'Next we write A:'
    ! do ialpha = 1,nsp_tot+1   !NEW: nsp here
    ! write(*,'(*(F14.7))')( Aglob(ialpha,ibeta) , ibeta=1,nsp_tot+1)
    ! enddo 
!===================================
!===================================
!Next we compute the vector of coefficients, b:

      write(*,*) 'NEXT begins the computation of b'
 !********************FIRST: CASE OF PURELY THREE-CENTERS INTEGRALS ***
       do ialp = 1,natoms
       rna(:) = ratom(:,ialp)
       indna = imass(ialp)

      ! Loop over the neighbors of each ialp.
      ! Now look at all common neighbor pairs which were figured out in main.
      ! The common neighbors were figured out in common_neighbors.f90
      !BEWARE!!  We must CHANGE this when we fully implement the SYMMETRIZATION
      !techniqe !!!!!!!!1
        do ineigh = 1, neigh_comn(ialp)     !Beware!! This index here tags pairs!
          mneigh = neigh_comm(ineigh,ialp)
      ! The second atom (jatom) is the mneigh'th neighbor of iatom.

        if (mneigh .ne. 0) then


          iatom = neigh_comj(1,ineigh,ialp)
           r1(:) = ratom(:,iatom) 
           in1 = imass(iatom)

           jatom = neigh_comj(2,ineigh,ialp)
           r2(:) = ratom(:,jatom)
           in2 = imass(jatom)
!We already have all the atoms involved in this integral. Next we set up stuff:

! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
           r21(:) = r2(:) - r1(:)
           y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))

! Find the unit vector in sigma direction.
           if (y .lt. 1.0d-05) then
            sighat(1) = 0.0d0
            sighat(2) = 0.0d0
            sighat(3) = 1.0d0
            write (*,*) ' There is an error here in assemble_3c.f '
            write (*,*) ' r1 = r2!!!! BAD!!!! '
           else
            sighat(:) = r21(:)/y
           end if

! Find rnabc = vector pointing from center of bondcharge to the neutral atom.
! The center of the bondcharge is at r1 + r21/2.  This gives us the distance
! dnabc (x value in 2D grid).
           rnabc(:) = rna(:) - (r1(:) + r21(:)/2.0d0)
           x = sqrt(rnabc(1)**2 + rnabc(2)**2 + rnabc(3)**2)

! Find the unit vector in rnabc direction.
           if (x .lt. 1.0d-05) then
            rhat(1) = 0.0d0
            rhat(2) = 0.0d0
            rhat(3) = 0.0d0
           else
            rhat(:) = rnabc(:)/x
           end if
           cost = sighat(1)*rhat(1) + sighat(2)*rhat(2) + sighat(3)*rhat(3)
           call epsilon (rhat, sighat, eps)
!All is ready to compute the three-centers integral
          interaction = 3
           do isorp = 1, nsp(indna) !NEW : nsp here !!!
            ialpha = isorp+degelec_sp(ialp) !NEW: degelec_sp here!!!

            call trescentros (interaction, isorp, isorpmax, in1, in2,   &
     &                        indna, x, y, cost, eps, rhomx, nspecies)

             

                 do inu = 1,num_orb(in2)
                  do imu = 1,num_orb(in1)
                   I = rhomx(imu,inu)
                   b(ialpha,1)=b(ialpha,1)+I*rho(imu,inu,mneigh,iatom)

                  enddo !end loop over inu
                 enddo !end loop over imu
           enddo !end loop over isorp

        endif !end if mneigh .ne. 0 
      enddo !end loop over ineigh
     enddo !end loop over ialp

!***************SECOND: TWO-CENTERS INTEGRALS **********

       do iatom = 1,natoms
        r1(:) = ratom(:,iatom)
        in1 = imass(iatom)
       
        ! Loop over the neighbors of each iatom.
         do ineigh = 1, neighn(iatom)          ! <==== loop over i's neighbors
          jatom = neigh_j(ineigh,iatom)
          r2(:) = ratom(:,jatom)
          in2 = imass(jatom)


          if (iatom .eq. jatom ) then    !ONE-CENTER CASE
              
            do isorp = 1,nsp(in1) !NEW: nso here!!!

               ialpha = isorp + degelec_sp(iatom) !NEW: degelec_sp here!!!
               interaction= 17
               y=0
               sighat(1) = 0.0d0
               sighat(2) = 0.0d0
               sighat(3) = 1.0d0

                call epsilon (r2, sighat, eps)
                call deps2cent (r1, r2, eps, deps)
              
                call doscentros (interaction, isorp, iforce, in1, in1, in1, &
     &                       y, eps, deps, rhomx, rhompx)
                
               do inu = 1, num_orb(in1)
                do imu = 1, num_orb(in1)
                    
                   I = rhomx(imu,inu)
                   b(ialpha,1) = b(ialpha,1) + I*rho(imu,inu,ineigh,iatom)

                enddo
               enddo
           
             enddo !end loop over isorp
                    
        
          else   !Here begins the PURELY TWO-CENTERS CASE:


       ! SET-UP STUFF
! ****************************************************************************
! Find r21 = vector pointing from r1 to r2, the two ends of the bondcharge.
! This gives us the distance dbc (or y value in the 2D grid).
           r21(:) = r2(:) - r1(:)
           y = sqrt(r21(1)*r21(1) + r21(2)*r21(2) + r21(3)*r21(3))

! Find the unit vector in sigma direction.
           if (y .lt. 1.0d-05) then
            sighat(1) = 0.0d0
            sighat(2) = 0.0d0
            sighat(3) = 1.0d0
           else
            sighat(:) = r21(:)/y
           end if

           call epsilon (r2, sighat, eps)
           call deps2cent (r1, r2, eps, deps)




! ****************************************************************************
! Left piece: den_ontopl <i|n_i|j> (den1,2) part of <i|n|j> (denmx)
           

           interaction = 15
           in3 = in1
           do isorp = 1, nsp(in3) !NEW: nsp here!!!
            call doscentros (interaction, isorp, iforce, in1, in3, in2, &
     &                       y, eps, deps, rhomx, rhompx)
               

            ialpha = isorp + degelec_sp(iatom) !NEW: degelec_sp here!!!

            do inu = 1, num_orb(in2)
             do imu = 1, num_orb(in1)

                 I = rhomx(inu,imu)
                 b(ialpha,1) = b(ialpha,1) + I*rho(imu,inu,ineigh,iatom)

             enddo
            enddo


           enddo !end of loop over isorp


! Right piece: den_ontopr <i|n_j|j> (den0) part of <i|n|j> (denmx)
          

           interaction = 16
           in3 = in2
           do isorp = 1, nsp(in3) !NEW: nsp here!!!
            call doscentros (interaction, isorp, iforce, in1, in3, in2, y,   &
     &                       eps, deps, rhomx, rhompx)
            

            ialpha = isorp + degelec_sp(jatom)    !NEW: degelec_sp here!!!

            do inu = 1, num_orb(in2)
             do imu = 1, num_orb(in1)
                 
                  I = rhomx(inu,imu)
                  b(ialpha,1) = b(ialpha,1)+I*rho(imu,inu,ineigh,iatom)

             enddo
            enddo


           enddo !end of loop over isorp

! Center piece: den_atomcase <i|n_j|i>
           

             interaction = 17
           
           do isorp = 1,nsp(in2) !NEW: nsp here!!!
             call doscentros (interaction, isorp, iforce, in1, in2, in1, y, &
      &                       eps, deps, rhomx, rhompx)


            ialpha = isorp + degelec_sp(jatom) !NEW: degelec_sp here!!!      

            do inu = 1, num_orb(in1)
             do imu = 1, num_orb(in1)

                  I = rhomx(inu,imu)
                  matom = neigh_self(iatom)
                  b(ialpha,1) = b(ialpha,1)+I*rho(imu,inu,matom,iatom)

             enddo
            enddo


           enddo !end of loop over isorp
        endif !end if on line 267 ( iatom .eq. jatom and corresponding else )
       enddo !end loop over ineigh
      enddo !end loop over iatom


      b(nsp_tot+1,1)=sum(nelectron)+qstate !Number of electrons here  !NEW: nsp_tot here!!!
write(*,*) 'b is:'
write(*,'(1000F14.7)') , (b(ialpha,1) ,ialpha=1,nsp_tot+1) !NEW: nsp_tot here!!!
!==================================
!==================================
!Next solve the system Ax=b


        Atemp=Aglob 
       ! bb=b

      call dgesv(nsp_tot+1,1,Atemp,nsp_tot+1,pivots,b,nsp_tot+1,info) !NEW:nsp_tot here!!!

       if (info .ne. 0) then
        
        write(*,*) 'Error while trying to solve linear system in charges_L2.f90'

       else

 !         write(*,*) 'The error between the charges before and the charges now
!is ',norm2(b-bb)


!The occupation numbers are now stored in the b vector.

!==================================
!Define the new Qout



          do iatom = 1,natoms
              in1 = imass(iatom)
              do issh = 1,nsp(in1)   !NEW: nsp here!!!
                 ialpha = issh + degelec_sp(iatom) !NEW: degelec_sp here!!!
                 write(*,*) 'Here at the END: in1= ',in1
                 write(*,*) 'And issh = ',issh
                 Qout(issh,iatom) = b(ialpha,1)
                 if ( in1 .eq. 1 ) then
                   write(*,*)'The OUTPUT CHARGES for this H are: ',Qout(1,iatom),Qout(2,iatom)
                 end if
                ! write(*,*) 'In charges_L2: issh,iatom,Qout', issh, iatom, Qout(issh,iatom)
              enddo
          enddo 
                 


       endif



write(*,*) 'we are exiting the subroutine charges_L2 for the ',Kscf,'th time'






return
end
