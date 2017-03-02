! copyright info:
!
!                             @Copyright 2004
!                           Fireball Committee
! Brigham Young University - James P. Lewis, Chair
! Arizona State University - Otto F. Sankey
! Universidad Autonoma de Madrid - Jose Ortega
! Lawrence Livermore National Laboratory - Kurt Glaesemann
! Universidad de Madrid - Pavel Jelinek
! Brigham Young University - Hao Wang

! Other contributors, past and present:
! Auburn University - Jian Jun Dong
! Arizona State University - Gary B. Adams
! Arizona State University - Kevin Schmidt
! Arizona State University - John Tomfohr
! Motorola, Physical Sciences Research Labs - Alex Demkov
! Motorola, Physical Sciences Research Labs - Jun Wang
! Ohio University - Dave Drabold
! University of Regensburg - Juergen Fritsch
 
!
! RESTRICTED RIGHTS LEGEND
! Use, duplication, or disclosure of this software and its documentation
! by the Government is subject to restrictions as set forth in subdivision
! { (b) (3) (ii) } of the Rights in Technical Data and Computer Software
! clause at 52.227-7013.
 
! Dassemble_olsxc_3c_mdet.f90
! Program Description
! ===========================================================================
!       This routine assembles the OFF-SITE OLSXC interactions and
! three-center forces.  The matrix elements look like <Psi(1)|V(3)|Psi(2)>.
!
!  OFF-SITE means that (1) and (2) are at different atoms.
!  (3) must be different from (1) and (2)
!
! The generic form of the matrix element is
! v = <mu(r-(r1+ratm)) ! v(r-ratm) ! nu(r-(r2+ratm))>.
! We want the derivatives with respect to ratm, r1 and r2
! These are (resp.) "a" , "b" , and "c" cases
!
! The derivative f3 looks like for the neutral atom case:
! f3na = - sum (all neighbors of atom alpha at (li,bi) but bi.ne.balph)
!    * sum (all neighbors m of (li,bi), and not having b value balph)
!    * rho(mu,nu,i,m)* deriv wrt balpha <i! v(balph) !j>.
! Note the minus sign to make it "force-like".
!
! This program gets fa, fb, and fc pieces from Dtrescentros.
!
! See notes
! "general theory of third party terms" for the definition on p. 2 of
! these three terms. Breifly fa=-d/dratm, fb=-d/d(bi), fc=-d/d(bj)
! where matrix elements are < psi(r-bi) ! v(r-ratm) ! psi(r-bj)>
!
! ===========================================================================
! JOM : adapted to calculate also gh_3c ( Grad <  mu | H | nu > ) for
! the non-adiabatic calculation,
! gh_3c (ix,katom,imu,inu,ineigh,iatom)
! gh_3c is the gradient wrt to katom of the 3-C-contribution to
! matrix element 
! < mu | H | nu >, mu in iatom, and nu in ineigh of iatom
!
! ===========================================================================

! ===========================================================================
! Code written by:
! P. Jelinek
! Department of Thin Films
! Institute of Physics AS CR
! Cukrovarnicka 10
! Prague 6, CZ-162  
! FAX +420-2-33343184
! Office telephone  +420-2-20318528
! email: jelinekp@fzu.cz
!
! J.Ortega & J.P.Lewis
! ===========================================================================
 
! Program Declaration
! ===========================================================================
! Old statement
!        subroutine Dassemble_olsxc_3c (natoms, nspecies, nprocs, iordern)
! MHL to include gausssians.
        subroutine Dassemble_olsxc_3c_mdet (nprocs, iordern, igauss)
        use configuration
        use constants_fireball
        use density
        use charges
        use dimensions
        use forces
        use interactions
        use neighbor_map
        use nonadiabatic
        implicit none
 
! Argument Declaration and Description
! ===========================================================================
! Input

        integer, intent (in) :: iordern
        integer, intent (in) :: nprocs
! MHL
        integer, intent (in) :: igauss
 
! Local Parameters and Data Declaration
! ===========================================================================
! Define in the module interactions
!        real, parameter ::  xc_overtol = 1.0d-6
 
! Local Variable Declaration and Description
! ===========================================================================
        integer ialp
        integer iatom
        integer iatomstart
        integer ibeta
        integer ierror
        integer imu
        integer inu
        integer in1
        integer in2
        integer indna
        integer index1
        integer index2
        integer ineigh
        integer interaction
        integer isorp
        integer issh
        integer ix
        integer jssh
        integer jatom
        integer jbeta
        integer l1
        integer l2
        integer mneigh
        integer my_proc
        integer n1
        integer n2
        integer natomsp

        integer itestg


        real cost
        real x
        real y
        real muxc
        real dmuxc
        real d2muxc
        real exc
        real dexc
        real d2exc
        real sx
        real sm
        real rho_av

        real, dimension (3, 3, 3) :: depsA
        real, dimension (3, 3, 3) :: depsB
        real, dimension (3, 3) :: eps
        real, dimension (nsh_max, nsh_max) :: rho_3c
        real, dimension (3, nsh_max, nsh_max) :: rhop_3ca
        real, dimension (3, nsh_max, nsh_max) :: rhop_3cb
        real, dimension (3, nsh_max, nsh_max) :: rhop_3cc
        real, dimension (3, numorb_max, numorb_max) :: rhoxpa
        real, dimension (3, numorb_max, numorb_max) :: rhoxpb
        real, dimension (3, numorb_max, numorb_max) :: rhoxpc
        real, dimension (3, nsh_max, nsh_max) :: rhompa
        real, dimension (3, nsh_max, nsh_max) :: rhompb
        real, dimension (3, nsh_max, nsh_max) :: rhompc
        real, dimension (3, numorb_max, numorb_max) :: rhoinpa
        real, dimension (3, numorb_max, numorb_max) :: rhoinpb
        real, dimension (3, numorb_max, numorb_max) :: rhoinpc
        real, dimension (3, numorb_max, numorb_max) :: avrhop_a
        real, dimension (3, numorb_max, numorb_max) :: avrhop_b
        real, dimension (3, numorb_max, numorb_max) :: avrhop_c
        real, dimension (3, numorb_max, numorb_max) :: mxca
        real, dimension (3, numorb_max, numorb_max) :: mxcb
        real, dimension (3, numorb_max, numorb_max) :: mxcc
        real, dimension (numorb_max, numorb_max) :: rhoin
        real, dimension (nsh_max, nsh_max) :: rhomm
        real, dimension (3) :: spm
        real, dimension (3) :: rhop_avb
        real, dimension (3) :: rhop_avc
        real, dimension (3) :: r1
        real, dimension (3) :: r2
        real, dimension (3) :: r21
        real, dimension (3) :: rhat
        real, dimension (3) :: rna
        real, dimension (3) :: rnabc
        real, dimension (3) :: sighat

! BTN communication domain
        integer MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
        common  /btnmpi/ MPI_BTN_WORLD, MPI_OPT_WORLD, MPI_BTN_WORLD_SAVE
 
! Procedure
! ========================================================================
! Initialize the force contributions to zero.
        f3xca = 0.0d0
        f3xcb = 0.0d0
        f3xcc = 0.0d0
! JOM
! JOM-info : if we merge gh_3c and gh_xc_3c into gh_3c, then we do not
! have to initialize here! (only in Dassemble_3c.f90)
! merge large arrays VLADA
     !   gh_xc_3c = 0.0d0

! Determine which atoms are assigned to this processor.
        if (iordern .eq. 1) then
         call MPI_COMM_RANK (MPI_BTN_WORLD, my_proc, ierror)
         natomsp = natoms/nprocs
         if (my_proc .lt. mod(natoms,nprocs)) then
          natomsp = natomsp + 1
          iatomstart = natomsp*my_proc + 1
         else
          iatomstart = (natomsp + 1)*mod(natoms,nprocs)                 &
                      + natomsp*(my_proc - mod(natoms,nprocs)) + 1
         end if
        else
         iatomstart = 1
         natomsp = natoms
        end if
 
!
!****************************************************************************
!                      T H R E E - C E N T E R   P A R T
!                                 S N X C
!****************************************************************************  
! Choose atom ialp in the central cell. This is the atom whose position
! we take the derivative, and is the atom who has the the neutral atom
! potential.

! Loop over the atoms in the central cell.
!$omp parallel do private (rna, indna, ineigh, mneigh, iatom, ibeta, r1, in1) &
!$omp&  private (jatom, jbeta, r2, in2, r21, y, sighat, rnabc, x, rhat, cost) &
!$omp&  private (eps, depsA, depsB, rho_3c, rhop_3ca, rhop_3cb, rhop_3cc )    &
!$omp&  private (rhoinpa, rhoinpb, rhoinpc, avrhop_a, avrhop_b, avrhop_c)     &
!$omp&  private (isorp, interaction, rhoxpa, rhoxpb, rhoxpc, rhomm, rhompa)   &
!$omp&  private (rhompb, rhompc, imu, inu, issh, jssh, sm, spm, n1, l1)       &
!$omp&  private (rhop_avb, rhop_avc, index1, n2, l2, index2, sx, mxcb, mxcc)  &
!$omp&  private (mxca, exc, rho_av, muxc, dexc, d2exc, dmuxc, d2muxc, rhoin)
        do ialp = iatomstart, iatomstart - 1 + natomsp
         rna(:) = ratom(:,ialp)
         indna = imass(ialp)

! Loop over the neighbors of each ialp.
! Now look at all common neighbor pairs which were figured out in main.
! The common neighbors were figured out in common_neighbors.f90
         do ineigh = 1, neigh_comn(ialp)
          mneigh = neigh_comm(ineigh,ialp)
! The second atom (jatom) is the mneigh'th neighbor of iatom.

          if (mneigh .ne. 0) then
           iatom = neigh_comj(1,ineigh,ialp)
           ibeta = neigh_comb(1,ineigh,ialp)
           r1(:) = ratom(:,iatom) + xl(:,ibeta)
           in1 = imass(iatom)

           jatom = neigh_comj(2,ineigh,ialp)
           jbeta = neigh_comb(2,ineigh,ialp)
           r2(:) = ratom(:,jatom) + xl(:,jbeta)
           in2 = imass(jatom)
!JOM-test
!          write(*,*)'----------------------------------------'
!          write(*,*)'ialp=',ialp,'iatom=',iatom,'jatom=',jatom
!          write(*,*)'----------------------------------------'
!JOM-test

! ****************************************************************************
!
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
            write (*,*) ' There is an error here in assemble_olsxc_off.f '
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
           call deps3center (r1, r2, r21, y, rna, rnabc, eps, depsA, depsB)
 
! ****************************************************************************
!
! PERFORM ACTUAL CALCULATIONS
! CALL TRESCENTROS FOR RHO_OFF and ITS DERIVATIVES, THREE-CENTER PART
! ****************************************************************************
           interaction = 3
           rho_3c = 0.0d0
           rhop_3ca = 0.0d0
           rhop_3cb = 0.0d0
           rhop_3cc = 0.0d0
           rhoinpa = 0.0d0
           rhoinpb = 0.0d0
           rhoinpc = 0.0d0
           avrhop_a = 0.0d0 
           avrhop_b = 0.0d0 
           avrhop_c = 0.0d0 
           do isorp = 1, nssh(indna)

! MHL : include gaussians
           if (igauss.eq.0) then
            call Dtrescentros (interaction, isorp, isorpmax, in1,       &
     &                         in2, indna, x, y, cost, eps, depsA,      &
     &                         depsB, rhat, sighat, rhoin, rhoxpa,      &
     &                         rhoxpb, rhoxpc, nspecies)

!
! The terms rhompa, rhompb, and rhompc are already force-like ( - ) !!
            call trescentrosS (isorp, isorpmax, in1, in2, indna, x, y,  &
     &                         cost, eps, rhomm, nspecies)


            call DtrescentrosS (isorp, isorpmax, in1, in2, indna, x, y, &
     &                          cost, rhat, sighat, rhomm, rhompa,      &
     &                          rhompb, rhompc, nspecies)

            else if (igauss.eq.1) then
            call DtrescentrosG_VXC (isorp, in1, in2, indna, x, y, cost, &
     &                              eps, depsA, depsB, rhat, sighat,    &
     &                              rhoin, rhoxpa, rhoxpb, rhoxpc, rcutoff)

! The terms rhompa, rhompb, and rhompc are already force-like ( - ) !!
            call trescentrosGS_VXC (isorp, in1, in2, indna, x, y, cost, &
     &                              eps, rhomm, rcutoff)

            call DtrescentrosGS_VXC (isorp, in1, in2, indna, x, y, cost,&
     &                               rhat, sighat, rhomm, rhompa,       &
     &                               rhompb, rhompc, rcutoff)
        end if

! rhoin is input density in crystal coordinates
! The terms rhoxpa, rhoxpb, and rhoxpc ARE force-like ( - ) !!
! We now calculate rhoinpa, rhoinpb, rhoinpc as the derivatives
! of rhoin (so, NOT FORCE-LIKE: another ( - )  )
! The same for rhop_3ca, rhop_3cb, rhop_3cc
            rhoin(:,:) = rho_off(:,:,mneigh,iatom)

            do inu = 1, num_orb(in2)
             do imu = 1, num_orb(in1)
              rhoinpa(:,imu,inu) =                                           &
     &         rhoinpa(:,imu,inu) - rhoxpa(:,imu,inu)*Qneutral(isorp,indna)
              rhoinpb(:,imu,inu) =                                           &
     &         rhoinpb(:,imu,inu) - rhoxpb(:,imu,inu)*Qneutral(isorp,indna)
              rhoinpc(:,imu,inu) =                                           &
     &         rhoinpc(:,imu,inu) - rhoxpc(:,imu,inu)*Qneutral(isorp,indna)
             end do
            end do

            do inu = 1, nssh(in2)
             do imu = 1, nssh(in1)
              rho_3c(imu,inu) = rho_3c(imu,inu)                              &
                                + rhomm(imu,inu)*Qneutral(isorp,indna)
              rhop_3ca(:,imu,inu) =                                          &
     &         rhop_3ca(:,imu,inu) - rhompa(:,imu,inu)*Qneutral(isorp,indna)
              rhop_3cb(:,imu,inu) =                                          &
     &         rhop_3cb(:,imu,inu) - rhompb(:,imu,inu)*Qneutral(isorp,indna)
              rhop_3cc(:,imu,inu) =                                          &
     &         rhop_3cc(:,imu,inu) - rhompc(:,imu,inu)*Qneutral(isorp,indna)
! JOM-test
!              write(*,*)'issh,jssh=',imu,inu
!           write(*,*)'rhompb=',rhompb(1,imu,inu)
!           write(*,*)'rhompc=',rhompc(1,imu,inu)
             end do
            end do

           end do

!----------------------------------------------------------------------------
! Now, assemble the derivative of the three center part of the
! average density
! Loop over shells.
           do issh = 1, nssh(in1)
            do jssh = 1, nssh(in2)
             sm = sm_mat(issh,jssh,mneigh,iatom)
             if (abs(sm) .gt. xc_overtol) then 
              spm(:) = spm_mat(:,issh,jssh,mneigh,iatom)
              avrhop_b (:,issh,jssh) = avrhop_b (:,issh,jssh) +              &
     &          (sm*rhop_3cb(:,issh,jssh) - rho_3c(issh,jssh)*spm(:))/(sm*sm)
              avrhop_c (:,issh,jssh) = avrhop_c (:,issh,jssh) +              &
     &         (sm*rhop_3cc(:,issh,jssh) + rho_3c(issh,jssh)*spm(:))/(sm*sm)
! JOM-test
!               write(*,*)'xxxxxxxxxxxxxxxxxxxxxx'
!               write(*,*)'issh,jssh=',issh,jssh
!               write(*,*)'sm,spm=',sm,spm(1)
!               write(*,*)'rho_3c=',rho_3c(issh,jssh)
!               write(*,*)'rhoinpb=',rhoinpb(1,issh,jssh)
!               write(*,*)'rhoinpc=',rhoinpc(1,issh,jssh)
!               write(*,*)'rhop_3cb=',rhop_3cb(1,issh,jssh)
!               write(*,*)'rhop_3cc=',rhop_3cc(1,issh,jssh)
!               write(*,*)'xxxxxxxxxxxxxxxxxxxxxx'
               
! JOM-test
             endif
            end do
           end do                                              

! ****************************************************************************
! SNXC - PART
! CALCULATE DERIVATIVE OF VXC MATRIX ELEMENT
! ****************************************************************************
! Loop over shells.
           n1 = 0
           do issh = 1, nssh(in1)
            l1 = lssh(issh,in1)
            n1 = n1 + l1 + 1
            n2 = 0
            do jssh = 1, nssh(in2)
             l2 = lssh(jssh,in2)
             n2 = n2 + l2 + 1
             rho_av =  arho_off(issh,jssh,mneigh,iatom)

             rhop_avb(:) =  avrhop_b(:,issh,jssh)
             rhop_avc(:) =  avrhop_c(:,issh,jssh)
             call cepal (rho_av, exc, muxc, dexc, d2exc, dmuxc, d2muxc)
             do index1 = -l1, l1
              imu = n1 + index1
              do index2 = -l2, l2
               inu = n2 + index2      
               sx = s_mat(imu,inu,mneigh,iatom)

! derivatives of XC matrix elements: mxca, mxcb and mxcc
               mxcb(:,imu,inu) =                                             &
     &           rhop_avb(:)*d2muxc*(rhoin(imu,inu) - rho_av*sx)             &
     &          + rhoinpb(:,imu,inu)*dmuxc    
               mxcc(:,imu,inu) =                                             &
     &          + rhop_avc(:)*d2muxc*(rhoin(imu,inu) - rho_av*sx)            &
     &          + rhoinpc(:,imu,inu)*dmuxc          
               mxca(:,imu,inu) = - mxcb(:,imu,inu) - mxcc(:,imu,inu)
! JOM-test
!              write(*,*)'imu,inu=',imu,inu
!              write(*,*)'d2muxc, dmuxc=',d2muxc, dmuxc
!              write(*,*)'rhoin=',rhoin(imu,inu)
!              write(*,*)'rhoinpb=',rhoinpb(1,imu,inu)
!              write(*,*)'rhoinpc=',rhoinpc(1,imu,inu)
!              write(*,*)'rho_av=',rho_av
!              write(*,*)'sx=',sx
!              write(*,*)'rhop_avb=',rhop_avb(1)
!              write(*,*)'rhop_avc=',rhop_avc(1)
!              write(*,*)'MXCB=',mxcb(1,imu,inu)
!              write(*,*)'MXCC=',mxcc(1,imu,inu)
!              write(*,*)'MXCA=',mxca(1,imu,inu)
! JOM-test
 
! End loop over shells.
              end do
             end do
             n2 = n2 + l2
            end do
            n1 = n1 + l1
           end do  
!$omp critical (Dxc3)
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
             do ix = 1, 3
              f3xca(ix,ialp) = f3xca(ix,ialp)                                &
     &         - rho(imu,inu,mneigh,iatom)*mxca(ix,imu,inu)
              f3xcb(ix,iatom) = f3xcb(ix,iatom)                              &
     &        - rho(imu,inu,mneigh,iatom)*mxcb(ix,imu,inu)
              f3xcc(ix,jatom) = f3xcc(ix,jatom)                              &
     &         - rho(imu,inu,mneigh,iatom)*mxcc(ix,imu,inu)
             end do
            end do
           end do   
! JOM-add : add to gh_3c or gh_xc_3c
           do inu = 1, num_orb(in2)
            do imu = 1, num_orb(in1)
! merge large arrays VLADA
!              gh_xc_3c(:,ialp,imu,inu,mneigh,iatom) =                           &
!    &        gh_xc_3c(:,ialp,imu,inu,mneigh,iatom) + mxca(:,imu,inu)
!              gh_xc_3c(:,iatom,imu,inu,mneigh,iatom) =                          &
!     &        gh_xc_3c(:,iatom,imu,inu,mneigh,iatom) + mxcb(:,imu,inu)
!              gh_xc_3c(:,jatom,imu,inu,mneigh,iatom) =                          &
!     &        gh_xc_3c(:,jatom,imu,inu,mneigh,iatom) + mxcc(:,imu,inu)
             gh_3c(:,ialp,imu,inu,mneigh,iatom) =                           &
    &        gh_3c(:,ialp,imu,inu,mneigh,iatom) + mxca(:,imu,inu)
             gh_3c(:,iatom,imu,inu,mneigh,iatom) =                          &
    &        gh_3c(:,iatom,imu,inu,mneigh,iatom) + mxcb(:,imu,inu)
             gh_3c(:,jatom,imu,inu,mneigh,iatom) =                          &
    &        gh_3c(:,jatom,imu,inu,mneigh,iatom) + mxcc(:,imu,inu)

            end do
           end do   
! JOM-end
!$omp end critical (Dxc3)
! ****************************************************************************
!           write(*,*)'**********GH_XC_3C**************'
!           write(*,*)'ialp, ineigh', ialp, ineigh, iatom, jatom
!           do inu = 1, num_orb(in2)
!            do imu = 1, num_orb(in1)
!       write(*,*)'gh_xc_3c=', gh_xc_3c(1,ialp,imu,inu,mneigh,iatom), imu,inu, ialp
!       write(*,*)'gh_xc_3c=', gh_xc_3c(1,iatom,imu,inu,mneigh,iatom), imu,inu, iatom
!       write(*,*)'gh_xc_3c=', gh_xc_3c(1,jatom,imu,inu,mneigh,iatom), imu,inu, jatom
!              end do
!              end do

! End loop over ialp and its common neighbors.
          end if
         end do
        end do

! Format Statements
! ===========================================================================
 
        return
      end subroutine Dassemble_olsxc_3c_mdet
