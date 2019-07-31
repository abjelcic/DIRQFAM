c======================================================================c

      subroutine fam_drhodkappa( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /fam/ omega_start, omega_end, delta_omega, omega_print,
     &             omega, gamma_smear,
     &             i_calculation_type, i_coulomb, i_pairing,
     &             J_multipole, K_multipole, ISO;

      CHARACTER fg_spx;
      common /simplex/ N_total         , N_blocks        ,
     &                 ia_spx(NBX)     , id_spx(NBX)     ,
     &                 nf_size(NBX)    , ng_size(NBX)    ,
     &                 nz_spx(NBSX,NBX), nr_spx(NBSX,NBX),
     &                 ml_spx(NBSX,NBX), fg_spx(NBSX,NBX);

      common /fam_energies/ E_fam( NTX , 2 );

      COMPLEX*16 v;
      common /v_matrix/ v( NBSX , NBSX , NBX , 2 );

      COMPLEX*16 u;
      common /u_matrix/ u( NBSX , NBSX , NBX , 2 );

      LOGICAL dh_nnz, dDelta_nnz, dkappa_nnz, f_nnz;
      common /nnz_blocks/ dh_nnz    ( NBX , NBX ),
     &                    dDelta_nnz( NBX , NBX ),
     &                    dkappa_nnz( NBX , NBX ),
     &                    f_nnz     ( NBX , NBX );

      COMPLEX*16 f20, f02;
      common /f02_f20_matrix/ f20( NTX , NTX , 2 ),
     &                        f02( NTX , NTX , 2 );

      COMPLEX*16 dh_1, dh_2;
      common /delta_h/ dh_1( NTX , NTX , 2 ),
     &                 dh_2( NTX , NTX , 2 );

      COMPLEX*16 dDelta_pl, dDelta_mi;
      common /dDelta/ dDelta_pl( NTX , NTX , 2 ),
     &                dDelta_mi( NTX , NTX , 2 );

      COMPLEX*16 drho_1, drho_2;
      common /delta_rho/ drho_1( NTX , NTX , 2 ),
     &                   drho_2( NTX , NTX , 2 );

      COMPLEX*16 dkappa_pl, dkappa_mi;
      common /delta_kappa/ dkappa_pl( NTX , NTX , 2 ),
     &                     dkappa_mi( NTX , NTX , 2 );

      COMPLEX*16 r20, p20,
     &           RcmPcm_commutator,
     &           lamR, lamP;
      common /spurious/ r20( NTX , NTX , 2 ),
     &                  p20( NTX , NTX , 2 ),
     &                  RcmPcm_commutator(2),
     &                  lamR, lamP;



      COMPLEX*16 z_tmp;
      COMPLEX*16   h20( NTX , NTX , 2 ),
     &             h02( NTX , NTX , 2 ),
     &           x_fam( NTX , NTX , 2 ),
     &           y_fam( NTX , NTX , 2 ),
     &            Temp( NTX , NTX     );



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_drhodkappa() **********************';
      write(6,*) '';
      endif


c----------------------------------------------------------------------c
c     NOTE:                                                            c
c     In code some matrices are altered by addition of a sign, or by   c
c     transposition to make some formulae simpler, and to avoid many   c
c     unnecessary transpositions of matrices in the code               c
c                                                                      c
c                                                                      c
c     Matrix in paper                   Matix in code                  c
c                                                                      c
c     f_1                               + f1_JK                        c
c     f_2                               + f2_JK                        c
c     f^{20}                            - f20                          c
c     f^{02}                            - transp(f02)                  c
c     \delta h^{20}(\omega)             - h20                          c
c     \delta h^{02}(\omega)             - transp(h02)                  c
c     x(\omega)                         - x_fam                        c
c     y(\omega)                         - transp(y_fam)                c
c     r^{20}                            - r20                          c
c     p^{20}                            - p20                          c
c                                                                      c
c                                                                      c
c     Other matrices have consitent definition with ones in the paper  c
c----------------------------------------------------------------------c






c-----Calculation of h20 = + hermconj(u) * dh_1                * v
c-----                     + hermconj(v) * transp(dh_2)        * u
c-----                     - hermconj(u) * dDelta_pl           * u
c-----                     + hermconj(v) * hermconj(dDelta_mi) * v
      h20 = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);

                  if( dh_nnz(ib1,ib2) ) then

                  call zgemm( 'c'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        u(1,1,ib1,it)             , NBSX        ,
     &                        dh_1(i0,j0,it)            , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        v(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        h20(i0,j0,it)             , NTX         );



                  call zgemm( 'c'         , 't'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        v(1,1,ib1,it)             , NBSX        ,
     &                        dh_2(j0,i0,it)            , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        u(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        h20(i0,j0,it)             , NTX         );

                  endif

                  if( dDelta_nnz(ib1,ib2) ) then

                  call zgemm( 'c'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( -1.D0 , 0.D0 )   ,
     &                        u(1,1,ib1,it)             , NBSX        ,
     &                        dDelta_pl(i0,j0,it)       , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        u(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        h20(i0,j0,it)             , NTX         );



                  call zgemm( 'c'         , 'c'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        v(1,1,ib1,it)             , NBSX        ,
     &                        dDelta_mi(j0,i0,it)       , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        v(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        h20(i0,j0,it)             , NTX         );

                  endif

              enddo
          enddo
      enddo






c-----Calculation of h02 = + hermconj(v) * dh_1                * u
c-----                     + hermconj(u) * transp(dh_2)        * v
c-----                     + hermconj(v) * dDelta_pl           * v
c-----                     - hermconj(u) * hermconj(dDelta_mi) * u
      h02 = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);

                  if( dh_nnz(ib1,ib2) ) then

                  call zgemm( 'c'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        v(1,1,ib1,it)             , NBSX        ,
     &                        dh_1(i0,j0,it)            , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        u(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        h02(i0,j0,it)             , NTX         );



                  call zgemm( 'c'         , 't'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        u(1,1,ib1,it)             , NBSX        ,
     &                        dh_2(j0,i0,it)            , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        v(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        h02(i0,j0,it)             , NTX         );

                  endif

                  if( dDelta_nnz(ib1,ib2) ) then

                  call zgemm( 'c'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        v(1,1,ib1,it)             , NBSX        ,
     &                        dDelta_pl(i0,j0,it)       , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        v(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        h02(i0,j0,it)             , NTX         );



                  call zgemm( 'c'         , 'c'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( -1.D0 , 0.D0 )   ,
     &                        u(1,1,ib1,it)             , NBSX        ,
     &                        dDelta_mi(j0,i0,it)       , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        u(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        h02(i0,j0,it)             , NTX         );

                  endif

              enddo
          enddo
      enddo






c-----Calculation of xmn = -( f20mn + h20mn )/( Em + En - hw - igamma )
c-----Calculation of ymn = -( f02mn + h02mn )/( Em + En + hw + igamma )
      do it = 1 , 2
          do j = 1 , N_total
              E_mu = E_fam(j,it);
              do i = 1 , N_total
                  E_nu = E_fam(i,it);

                  x_fam(i,j,it) = - ( f20(i,j,it) + h20(i,j,it) ) /
     &                   COMPLEX( E_mu + E_nu - omega , - gamma_smear );

                  y_fam(i,j,it) = - ( f02(i,j,it) + h02(i,j,it) ) /
     &                   COMPLEX( E_mu + E_nu + omega , + gamma_smear );

              enddo
          enddo
      enddo






c-----Elimination of spurious Nambu-Goldstone mode
      if( J_multipole .eq. 1 ) then

          lamR = COMPLEX( 0.D0 , 0.D0 );
          lamP = COMPLEX( 0.D0 , 0.D0 );
          do it = 1 , 2
              do j = 1 , N_total
                  do i = 1 , N_total
                      lamR = lamR + ( - x_fam(i,j,it) + y_fam(i,j,it) )
     &                            * DCONJG( p20(i,j,it) );
                      lamP = lamP + ( + x_fam(i,j,it) + y_fam(i,j,it) )
     &                            * DCONJG( r20(i,j,it) );
                  enddo
              enddo
          enddo
          lamR = lamR / ( RcmPcm_commutator(1) + RcmPcm_commutator(2) );
          lamP = lamP / ( RcmPcm_commutator(1) + RcmPcm_commutator(2) );


          do it = 1 , 2
              do j = 1 , N_total
                  do i = 1 , N_total
                      x_fam(i,j,it) = x_fam(i,j,it) - lamR*r20(i,j,it)
     &                                              - lamP*p20(i,j,it);
                      y_fam(i,j,it) = y_fam(i,j,it) + lamR*r20(i,j,it)
     &                                              - lamP*p20(i,j,it);
                  enddo
              enddo
          enddo

      endif






c-----Calculation of drho_1 = + u * x_fam * hermconj(v)
c-----                        + v * y_fam * hermconj(u)
      drho_1 = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);

                  if( dh_nnz(ib1,ib2) ) then

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        u(1,1,ib1,it)             , NBSX        ,
     &                        x_fam(i0,j0,it)           , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'c'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        v(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        drho_1(i0,j0,it)          , NTX         );



                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        v(1,1,ib1,it)             , NBSX        ,
     &                        y_fam(i0,j0,it)           , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'c'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        u(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        drho_1(i0,j0,it)          , NTX         );

                  endif

              enddo
          enddo
      enddo






c-----Calculation of drho_2 = transp( + v * x_fam * hermconj(u)
c-----                                + u * y_fam * hermconj(v) )
      drho_2 = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);

                  if( dh_nnz(ib1,ib2) ) then

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        v(1,1,ib1,it)             , NBSX        ,
     &                        x_fam(i0,j0,it)           , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'c'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        u(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        drho_2(i0,j0,it)          , NTX         );



                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        u(1,1,ib1,it)             , NBSX        ,
     &                        y_fam(i0,j0,it)           , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'c'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        v(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        drho_2(i0,j0,it)          , NTX         );

                  endif

              enddo
          enddo

          do j = 1 , N_total
              do i = 1 , j-1
                          z_tmp  = drho_2(i,j,it);
                  drho_2(i,j,it) = drho_2(j,i,it);
                  drho_2(j,i,it) = z_tmp;
              enddo
          enddo

      enddo






c-----Calculation of dkappa_pl = - u * x_fam * hermconj(u)
c-----                           + v * y_fam * hermconj(v)
      dkappa_pl = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);

                 !if( dh_nnz    (ib1,ib2) ) then
                  if( dkappa_nnz(ib1,ib2) ) then

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( -1.D0 , 0.D0 )   ,
     &                        u(1,1,ib1,it)             , NBSX        ,
     &                        x_fam(i0,j0,it)           , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'c'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        u(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        dkappa_pl(i0,j0,it)       , NTX         );



                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        v(1,1,ib1,it)             , NBSX        ,
     &                        y_fam(i0,j0,it)           , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'c'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        v(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        dkappa_pl(i0,j0,it)       , NTX         );

                  endif

              enddo
          enddo
      enddo






c-----Calculation of dkappa_mi = + v * hermconj(x_fam) * hermconj(v)
c-----                           - u * hermconj(y_fam) * hermconj(u)
      dkappa_mi = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);

                 !if( dh_nnz    (ib1,ib2) ) then
                  if( dkappa_nnz(ib1,ib2) ) then

                  call zgemm( 'n'         , 'c'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        v(1,1,ib1,it)             , NBSX        ,
     &                        x_fam(j0,i0,it)           , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'c'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        v(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )   ,
     &                        dkappa_mi(i0,j0,it)       , NTX         );



                  call zgemm( 'n'         , 'c'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( -1.D0 , 0.D0 )   ,
     &                        u(1,1,ib1,it)             , NBSX        ,
     &                        y_fam(j0,i0,it)           , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'c'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        u(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        dkappa_mi(i0,j0,it)       , NTX         );

                  endif

              enddo
          enddo
      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_drhodkappa() ************************';
      write(6,*) '';
      endif

      return;
      end;
