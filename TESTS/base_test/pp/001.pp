!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#info ../tools/psgen.com ncpp 15:20:17 Jul 01 1999 gayathri   
                                                                      
fhi pseudopotential tool gncpp - version rev 051697                   
                                                                      
               chemical symbol  H                                     
                nuclear charge   1.00                                 
                  total charge   0.00                                 
         number of core states   0                                    
      number of valence states   1                                    
    exchange-correlation model   9  GGA X Becke C Lee/Yang/Parr 
      scalar-relativistic mode                                        
        parameters radial mesh   387    1.024700  0.625000E-02        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    9   0.00                            ! iexc , exmix
    1                                   ! nshells
    0                                   ! L values
    1.000000                            ! Zval
    50.00000                            ! alpha   (Antes era 800.00000) 
    0.200000                            ! Rcut_PP   
            5
   0.0000000E+00   0.000000    
   4.9999999E-03   0.000000    
   9.9999998E-03   0.000000    
   1.5000000E-02   0.000000    
   2.0000000E-02   0.000000    
 L=0 Npoint=    5
   0.0000000E+00  0.0000000E+00
   4.9999999E-03  0.0000000E+00
   9.9999998E-03  0.0000000E+00
   1.5000000E-02  0.0000000E+00
   2.0000000E-02  0.0000000E+00
 L=0 Npoint=    5 cl=     0.0000000
   0.0000000E+00  0.0000000E+00
   4.9999999E-03  0.0000000E+00
   9.9999998E-03  0.0000000E+00
   1.5000000E-02  0.0000000E+00
   2.0000000E-02  0.0000000E+00
#info ../tools/psgen.com ncpp 15:20:17 Jul 01 1999 gayathri

fhi pseudopotential tool gncpp - version rev 051697

               chemical symbol  H 
                nuclear charge   1.00
                  total charge   0.00
         number of core states   0
      number of valence states   1
    exchange-correlation model   9  GGA X Becke C Lee/Yang/Parr     
      scalar-relativistic mode
        parameters radial mesh   387    1.024700  0.625000E-02

       === all-electron atom ===

<        n     l      occupation  eigenvalue(eV)

<  1     1     0        1.0000           -6.3588

                                  (Hartree a.u.)
                  total energy          -0.44591
                kinetic energy           0.42486
                orbital energy          -0.23368
                coulomb energy          -0.92083
                hartree energy           0.28279
   exchange-correlation energy          -0.23273
           xc potential energy          -0.30329
          number of iterations                28   convergence  0.0E+00
            integrated density           1.00000
  ... 1st derivative test 1 =            1.00000
  ... 2nd derivative test 1 =            1.00000

 gncpp - all-electron atom done

    === HAMANN mode ===    h    

  l  n     radius:     node      peak       default core
x 0  1                0.006      1.050      0.420
x 1  2                0.006      0.006      0.420

          === pseudo atom ===

  l  type  rcore       rmatch          eigenvalue(eV)      norm test   slope test
                                 all-electron     pseudo
  0  hamann                   0.4154499   2.1832160  -6.3587890  -6.3587890   1.0000000   1.0000035
  1  hamann                   0.4154499   1.0247027  -6.3587890  -6.3587890   1.0000000   1.0000000


                                  (Hartree a.u.)
                  total energy          -0.44587
                kinetic energy           0.42050
              potential energy          -0.91647
                hartree energy           0.28272
                     xc energy          -0.23262
    integrated valence density           1.00000
   ... 1st derivative test 1 =           1.00000
   ... 2nd derivative test 1 =           1.00000

 gncpp - done for input @

@  1.00  0  1  30.00E+00 : z nc nv iexc rnlc
@     1  0   1.00        : n l f
@ 1  h                   : ltmx spptype
