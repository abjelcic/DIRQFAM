c======================================================================c

      subroutine init_multipole( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      CHARACTER*2 nucnam;
      common /nucnuc/ amas, nneu, npro, nmas, nucnam;
      common /mathco/ zero, one, two, half, third, pi;
      common /gaussh/ xh(0:NGH), wh(0:NGH), zb(0:NGH);
      common /gaussl/ xl(0:NGL), wl(0:NGL), sxl(0:NGL), rb(0:NGL);

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

      COMPLEX*16 v;
      common /v_matrix/ v( NBSX , NBSX , NBX , 2 );

      COMPLEX*16 u;
      common /u_matrix/ u( NBSX , NBSX , NBX , 2 );

      LOGICAL dh_nnz, dDelta_nnz, dkappa_nnz, f_nnz;
      common /nnz_blocks/ dh_nnz    ( NBX , NBX ),
     &                    dDelta_nnz( NBX , NBX ),
     &                    dkappa_nnz( NBX , NBX ),
     &                    f_nnz     ( NBX , NBX );

      common /basis/ qhql ( -NGH:NGH , 1:NGL , NTX ),
     &               qh1ql( -NGH:NGH , 1:NGL , NTX ),
     &               qhql1( -NGH:NGH , 1:NGL , NTX ),
     &               wqhql(    1:NGH , 1:NGL , NTX );

      COMPLEX*16 f1_JK, f2_JK;
      common /f_matrix/ f1_JK( NTX , NTX , 2 ),
     &                  f2_JK( NTX , NTX , 2 );

      COMPLEX*16 f20, f02;
      common /f02_f20_matrix/ f20( NTX , NTX , 2 ),
     &                        f02( NTX , NTX , 2 );



      CHARACTER fg1, fg2;
      COMPLEX*16 Temp( NTX , NTX );
      REAL*8 f( -NGH:NGH , 1:NGL , 2 , 0:J_max , 0:J_MAX );
      REAL*8 fac_iso(2);



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN init_multipole() ***********************';
      write(6,*) '';
      endif


c----------------------------------------------------------------------c
c     Excitation operator has the following form:                      c
c                                                                      c
c     f = fac_iso(t_z) * 1/sqrt(2+2*delta_{K,0})                       c
c                      * r^J * ( Y_{J,K} + (-1)^K * Y_{J,-K} ),        c
c                                                                      c
c     with only J=0,K=0 exception: f = fac_iso(t_z) * r^2.             c
c                                                                      c
c                                                                      c
c     For isoscalar excitation: fac_iso(t_z) = +1, and for isovector:  c
c     fac_iso(+1/2) = +1 (protons), fac_iso(-1/2) = -1 (neutrons),     c
c     with only exception for J=1 excitations, where we used:          c
c     fac_iso(+1/2) = + N/(N+Z) (protons),                             c
c     fac_iso(-1/2) = - Z/(N+Z) (neutrons).                            c
c                                                                      c
c     We coded it=1 for neut.(t_z=-1/2), and it=2 for prot.(t_z=+1/2)  c
c----------------------------------------------------------------------c


      if( J_multipole .gt. J_MAX ) then
          stop 'Error: Wrong value of J!';
      endif
      if( K_multipole .gt. J_multipole ) then
          stop 'Error: Wrong value of K!';
      endif

      Temp  = COMPLEX( 0.D0 , 0.D0 );
      f1_JK = COMPLEX( 0.D0 , 0.D0 );
      f2_JK = COMPLEX( 0.D0 , 0.D0 );
      f20   = COMPLEX( 0.D0 , 0.D0 );
      f02   = COMPLEX( 0.D0 , 0.D0 );






c-----Initializing multipole operators in coordinate space
      fac_iso(1) = + 1.D0;
      fac_iso(2) = + 1.D0;
      if( ISO .eq. 1 ) then
          if( J_multipole .eq. 1 ) then
              fac_iso(1) = - DBLE(npro) / DBLE(npro+nneu);
              fac_iso(2) = + DBLE(nneu) / DBLE(npro+nneu);
          else
              fac_iso(1) = - 1.D0;
              fac_iso(2) = + 1.D0;
          endif
      endif

      do it = 1 , 2
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;

                  r = rb(il);
                  z = DBLE(ISIGN(1,ih)) * zb(ABS(ih));

                  fac00 = + 1.D0 ! sometimes 1/sqrt(4*pi) is used
                  fac10 = + 0.5D0   * DSQRT( 3.D0   / pi );
                  fac11 = - 0.5D0   * DSQRT( 3.D0   / pi );
                  fac20 = + 0.25D0  * DSQRT( 5.D0   / pi );
                  fac21 = - 0.5D0   * DSQRT( 15.D0  / pi );
                  fac22 = + 0.25D0  * DSQRT( 15.D0  / pi );
                  fac30 = + 0.25D0  * DSQRT( 7.D0   / pi );
                  fac31 = - 0.125D0 * DSQRT( 42.D0  / pi );
                  fac32 = + 0.25D0  * DSQRT( 105.D0 / pi );
                  fac33 = - 0.125D0 * DSQRT( 70.D0  / pi );

                  fac00 = fac00 * fac_iso(it);
                  fac10 = fac10 * fac_iso(it);
                  fac11 = fac11 * fac_iso(it);
                  fac20 = fac20 * fac_iso(it);
                  fac21 = fac21 * fac_iso(it);
                  fac22 = fac22 * fac_iso(it);
                  fac30 = fac30 * fac_iso(it);
                  fac31 = fac31 * fac_iso(it);
                  fac32 = fac32 * fac_iso(it);
                  fac33 = fac33 * fac_iso(it);

                  f(ih,il,it,0,0) = fac00 * (r*r+z*z);      ! cos(0*phi)
                  f(ih,il,it,1,0) = fac10 * z;              ! cos(0*phi)
                  f(ih,il,it,1,1) = fac11 * r;              ! cos(1*phi)
                  f(ih,il,it,2,0) = fac20 * (2*z*z-r*r);    ! cos(0*phi)
                  f(ih,il,it,2,1) = fac21 * r*z;            ! cos(1*phi)
                  f(ih,il,it,2,2) = fac22 * r*r;            ! cos(2*phi)
                  f(ih,il,it,3,0) = fac30 * z*(2*z*z-3*r*r);! cos(0*phi)
                  f(ih,il,it,3,1) = fac31 * r*(4*z*z-r*r);  ! cos(1*phi)
                  f(ih,il,it,3,2) = fac32 * r*r*z;          ! cos(2*phi)
                  f(ih,il,it,3,3) = fac33 * r*r*r;          ! cos(3*phi)

              enddo
          enddo
      enddo






c-----Calculation of f1_JK matrix
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  j0 = ia_spx(ib2);
                  j1 = ia_spx(ib2)+id_spx(ib2)-1;
                  do j = j0 , j1
                      fg2 = fg_spx(j-j0+1,ib2);
                      nz2 = nz_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);

                      i0 = ia_spx(ib1);
                      i1 = ia_spx(ib1)+id_spx(ib1)-1;
                      do i = i0 , i1
                          fg1 = fg_spx(i-i0+1,ib1);
                          nz1 = nz_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1 .ne. fg2 ) CYCLE;
                          if( abs(ml1-ml2) .ne. K_multipole ) CYCLE;


                          JJ = J_multipole;
                          KK = K_multipole;
                          acc = 0.D0;
                          if( MOD(nz1+nz2,2) .eq. 0 ) then
                              do il = 1 , NGL
                                  do ih = 1 , NGH
                                      x = +f( +ih , il , it , JJ , KK )
     &                                    +f( -ih , il , it , JJ , KK );

                                      acc = acc + x * wqhql(ih,il,i)
     &                                              * wqhql(ih,il,j);
                                  enddo
                              enddo
                          else
                              do il = 1 , NGL
                                  do ih = 1 , NGH
                                      x = +f( +ih , il , it , JJ , KK )
     &                                    -f( -ih , il , it , JJ , KK );

                                      acc = acc + x * wqhql(ih,il,i)
     &                                              * wqhql(ih,il,j);
                                  enddo
                              enddo
                          endif


                          if( K_multipole .ne. 0 ) then
                              acc = acc * 0.5D0;
                          endif

                          f1_JK(i,j,it) = COMPLEX( acc , 0.D0 );

                      enddo

                  enddo

              enddo
          enddo
      enddo






c-----Calculation of f2_JK matrix
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  j0 = ia_spx(ib2);
                  j1 = ia_spx(ib2)+id_spx(ib2)-1;
                  do j = j0 , j1
                      fg2 = fg_spx(j-j0+1,ib2);
                      nz2 = nz_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);

                      i0 = ia_spx(ib1);
                      i1 = ia_spx(ib1)+id_spx(ib1)-1;
                      do i = i0 , i1
                          fg1 = fg_spx(i-i0+1,ib1);
                          nz1 = nz_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1 .ne. fg2 ) CYCLE;
                          if( abs(ml1-ml2) .ne. K_multipole ) CYCLE;


                          JJ = J_multipole;
                          KK = K_multipole;
                          acc = 0.D0;
                          if( MOD(nz1+nz2,2) .eq. 0 ) then
                              do il = 1 , NGL
                                  do ih = 1 , NGH
                                      x = +f( +ih , il , it , JJ , KK )
     &                                    +f( -ih , il , it , JJ , KK );

                                      acc = acc + x * wqhql(ih,il,i)
     &                                              * wqhql(ih,il,j);
                                  enddo
                              enddo
                          else
                              do il = 1 , NGL
                                  do ih = 1 , NGH
                                      x = +f( +ih , il , it , JJ , KK )
     &                                    -f( -ih , il , it , JJ , KK );

                                      acc = acc + x * wqhql(ih,il,i)
     &                                              * wqhql(ih,il,j);
                                  enddo
                              enddo
                          endif


                          if( K_multipole .ne. 0 ) then
                              acc = acc * 0.5D0;
                          endif

                          f2_JK(i,j,it) = COMPLEX( acc , 0.D0 );

                      enddo

                  enddo

              enddo
          enddo
      enddo






c-----Calculation of f20 = + hermconjg(u) * f1_JK         * v
c-----                     + hermconjg(v) * transp(f2_JK) * u
      do it = 1 , 2

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  if( f_nnz(ib1,ib2) .eqv. .false. ) CYCLE;

                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);

                  call zgemm( 'c'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        u(1,1,ib1,it)             , NBSX        ,
     &                        f1_JK(i0,j0,it)           , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        v(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        f20(i0,j0,it)             , NTX         );



                  call zgemm( 'c'         , 't'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        v(1,1,ib1,it)             , NBSX        ,
     &                        f2_JK(j0,i0,it)           , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        u(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        f20(i0,j0,it)             , NTX         );


              enddo
          enddo



      enddo






c-----Calculation of f02 = + hermconjg(v) * f1_JK         * u
c-----                     + hermconjg(u) * transp(f2_JK) * v
      do it = 1 , 2

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  if( f_nnz(ib1,ib2) .eqv. .false. ) CYCLE;

                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);

                  call zgemm( 'c'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        v(1,1,ib1,it)             , NBSX        ,
     &                        f1_JK(i0,j0,it)           , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        u(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        f02(i0,j0,it)             , NTX         );



                  call zgemm( 'c'         , 't'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        u(1,1,ib1,it)             , NBSX        ,
     &                        f2_JK(j0,i0,it)           , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        v(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        f02(i0,j0,it)             , NTX         );


              enddo
          enddo



      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END init_multipole() ************************';
      write(6,*) '';
      endif

      return;
      end;
