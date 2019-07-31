c======================================================================c

      subroutine fam_ddensdcurr( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /mathco/ zero, one, two, half, third, pi;
      common /baspar/ hom, hb0, b0;
      common /defbas/ beta0, q, bp, bz;
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

      common /basis/ qhql ( -NGH:NGH , 1:NGL , NTX ),
     &               qh1ql( -NGH:NGH , 1:NGL , NTX ),
     &               qhql1( -NGH:NGH , 1:NGL , NTX ),
     &               wqhql(    1:NGH , 1:NGL , NTX );

      common /k_PHI/ k_PHI, k_PHIr, k_dzPHI, k_drPHI;
      common /PHI/ PHI_U    (    NTX , KTRUNC ),
     &             PHI_SVt  ( KTRUNC , NCOORD ),
     &
     &             PHIr_U   (    NTX , KTRUNC ),
     &             PHIr_SVt ( KTRUNC , NCOORD ),
     &
     &             dzPHI_U  (    NTX , KTRUNC ),
     &             dzPHI_SVt( KTRUNC , NCOORD ),
     &
     &             drPHI_U  (    NTX , KTRUNC ),
     &             drPHI_SVt( KTRUNC , NCOORD );

      LOGICAL dh_nnz, dDelta_nnz, dkappa_nnz, f_nnz;
      common /nnz_blocks/ dh_nnz    ( NBX , NBX ),
     &                    dDelta_nnz( NBX , NBX ),
     &                    dkappa_nnz( NBX , NBX ),
     &                    f_nnz     ( NBX , NBX );

      COMPLEX*16 drho_1, drho_2;
      common /delta_rho/ drho_1( NTX , NTX , 2 ),
     &                   drho_2( NTX , NTX , 2 );

      COMPLEX*16 drho_v, drho_s, ldrho_v, ldrho_s;
      common /ind_dens/ drho_v ( -NGH:NGH , 1:NGL , 2 ),
     &                  drho_s ( -NGH:NGH , 1:NGL     ),
     &                  ldrho_v( -NGH:NGH , 1:NGL , 2 ),
     &                  ldrho_s( -NGH:NGH , 1:NGL     );

      COMPLEX*16 dj_r, dj_p, dj_z;
      common /ind_curr/ dj_r( -NGH:NGH , 1:NGL , 2 ),
     &                  dj_p( -NGH:NGH , 1:NGL , 2 ),
     &                  dj_z( -NGH:NGH , 1:NGL , 2 );



      REAL*8 Ar( NTX , NTX );
      REAL*8 Ai( NTX , NTX );
      REAL*8 AUX1(    NTX , KTRUNC );
      REAL*8 AUX2( KTRUNC , KTRUNC );
      COMPLEX*16 D( NCOORD );
      CHARACTER fg1, fg2;
      COMPLEX*16 z;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_ddensdcurr() ***********************';
      write(6,*) '';
      endif

c-----Note that drho_1 and drho_2, have block structure given by dh_nnz






c-----Selection rules test
      if( .false. ) then

      K = K_multipole;
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh_nnz(ib1,ib2) ) then

                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);

                  do j = j0 , j0-1+id_spx(ib2)
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);
                      do i = i0 , i0-1+id_spx(ib1)
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1 .eq. fg2 ) then
                              z = drho_1(i,j,it) + drho_2(i,j,it);
                              if( ABS(z) .gt. 1.D-10 ) then
                                  if( abs(ml1-ml2) .ne. K ) then
                                      write(6,*)'Should be zero: ', z;
                                      stop 'Error: |ml1-ml2| =/= K!';
                                  endif
                              endif
                          endif

                          if( fg1 .ne. fg2 ) then
                              z = + drho_1(i,j,it) + drho_2(i,j,it)
     &                            - drho_1(j,i,it) - drho_2(j,i,it);
                              if( ABS(z) .gt. 1.D-10 ) then
                                  if( abs(ml1+ml2+1) .ne. K ) then
                                      write(6,*)'Should be zero: ', z;
                                      stop 'Error: |ml1+ml2+1| =/= K!';
                                  endif
                              endif
                              z = + drho_1(i,j,it) - drho_2(i,j,it)
     &                            + drho_1(j,i,it) - drho_2(j,i,it);
                              if( ABS(z) .gt. 1.D-10 ) then
                                  if( abs(ml1-ml2) .ne. K ) then
                                      write(6,*)'Should be zero: ', z;
                                      stop 'Error: |ml1-ml2| =/= K!';
                                  endif
                              endif
                          endif


                      enddo
                  enddo


                  endif
              enddo
          enddo
      enddo

      endif






c-----Calculation of drho_v
      drho_v = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2

          Ar = 0.D0;
          Ai = 0.D0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh_nnz(ib1,ib2) ) then
                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);
                  do j = j0 , j0-1+id_spx(ib2)
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);
                      do i = i0 , i0-1+id_spx(ib1)
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1.ne.fg2  .or.
     &                        abs(ml1-ml2).ne.K_multipole ) then
                              CYCLE;
                          endif

                          z = drho_1(i,j,it) + drho_2(i,j,it);
                          Ar(i,j) = DREAL(z);
                          Ai(i,j) = DIMAG(z);
                      enddo
                  enddo
                  endif
              enddo
          enddo


          call diagUVtAUV(  N_total , NCOORD ,   k_PHI , N_blocks ,
     &                           Ar ,    NTX ,      Ai ,      NTX ,
     &                       dh_nnz ,    NBX ,
     &                       ia_spx , id_spx ,
     &                        PHI_U ,    NTX , PHI_SVt ,   KTRUNC , D ,
     &                         AUX1 ,    NTX ,    AUX2 ,   KTRUNC     );

          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  drho_v(ih,il,it) = D(ihl) / (2.D0*pi);
              enddo
          enddo

      enddo






c-----Calculation of drho_s (protons + neutrons)
      drho_s = COMPLEX( 0.D0 , 0.D0 );
      do it = 3 , 3

          Ar = 0.D0;
          Ai = 0.D0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh_nnz(ib1,ib2) ) then
                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);
                  do j = j0 , j0-1+id_spx(ib2)
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);
                      do i = i0 , i0-1+id_spx(ib1)
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1.ne.fg2  .or.
     &                        abs(ml1-ml2).ne.K_multipole ) then
                              CYCLE;
                          endif

                          z =  + drho_1(i,j,1) + drho_2(i,j,1)
     &                         + drho_1(i,j,2) + drho_2(i,j,2);
                          if( fg1.eq.'g' .and. fg2.eq.'g' ) then
                              z = -z;
                          endif
                          Ar(i,j) = DREAL(z);
                          Ai(i,j) = DIMAG(z);
                      enddo
                  enddo
                  endif
              enddo
          enddo


          call diagUVtAUV(  N_total , NCOORD ,   k_PHI , N_blocks ,
     &                           Ar ,    NTX ,      Ai ,      NTX ,
     &                       dh_nnz ,    NBX ,
     &                       ia_spx , id_spx ,
     &                        PHI_U ,    NTX , PHI_SVt ,   KTRUNC , D ,
     &                         AUX1 ,    NTX ,    AUX2 ,   KTRUNC     );

          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  drho_s(ih,il) = D(ihl) / (2.D0*pi);
              enddo
          enddo

      enddo






c-----Calculation of dj_r
      dj_r = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2

          Ar = 0.D0;
          Ai = 0.D0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh_nnz(ib1,ib2) ) then
                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);
                  do j = j0 , j0-1+id_spx(ib2)
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);
                      do i = i0 , i0-1+id_spx(ib1)
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1.eq.fg2  .or.
     &                        abs(ml1+ml2+1).ne.K_multipole ) then
                              CYCLE;
                          endif

                          z = drho_1(i,j,it) + drho_2(i,j,it);
                          if( fg1.eq.'g' .and. fg2.eq.'f' ) then
                              z = -z;
                          endif
                          Ar(i,j) = DREAL(z);
                          Ai(i,j) = DIMAG(z);
                      enddo
                  enddo
                  endif
              enddo
          enddo


          call diagUVtAUV(  N_total , NCOORD ,   k_PHI , N_blocks ,
     &                           Ar ,    NTX ,      Ai ,      NTX ,
     &                       dh_nnz ,    NBX ,
     &                       ia_spx , id_spx ,
     &                        PHI_U ,    NTX , PHI_SVt ,   KTRUNC , D ,
     &                         AUX1 ,    NTX ,    AUX2 ,   KTRUNC     );

          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  dj_r(ih,il,it) = D(ihl) / COMPLEX( 0.D0 , 2.D0*pi );
              enddo
          enddo

      enddo






c-----Calculation of dj_p
      dj_p = COMPLEX( 0.D0 , 0.D0 );
      if( K_multipole .ne. 0 ) then
      do it = 1 , 2

          Ar = 0.D0;
          Ai = 0.D0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh_nnz(ib1,ib2) ) then
                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);
                  do j = j0 , j0-1+id_spx(ib2)
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);
                      do i = i0 , i0-1+id_spx(ib1)
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1.eq.fg2  .or.
     &                        abs(ml1+ml2+1).ne.K_multipole ) then
                              CYCLE;
                          endif

                          z = drho_1(i,j,it) + drho_2(i,j,it);
                          if( fg1.eq.'g' .and. fg2.eq.'f' ) then
                              z = -z;
                          endif
                          if( ml1+ml2+1 .lt. 0 ) then
                              z = -z;
                          endif
                          Ar(i,j) = DREAL(z);
                          Ai(i,j) = DIMAG(z);
                      enddo
                  enddo
                  endif
              enddo
          enddo


          call diagUVtAUV(  N_total , NCOORD ,   k_PHI , N_blocks ,
     &                           Ar ,    NTX ,      Ai ,      NTX ,
     &                       dh_nnz ,    NBX ,
     &                       ia_spx , id_spx ,
     &                        PHI_U ,    NTX , PHI_SVt ,   KTRUNC , D ,
     &                         AUX1 ,    NTX ,    AUX2 ,   KTRUNC     );

          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  dj_p(ih,il,it) = D(ihl) / COMPLEX( 0.D0 , -2.D0*pi );
              enddo
          enddo

      enddo
      endif






c-----Calculation of dj_z
      dj_z = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2

          Ar = 0.D0;
          Ai = 0.D0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh_nnz(ib1,ib2) ) then
                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);
                  do j = j0 , j0-1+id_spx(ib2)
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);
                      do i = i0 , i0-1+id_spx(ib1)
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1.eq.fg2  .or.
     &                        abs(ml1-ml2).ne.K_multipole ) then
                              CYCLE;
                          endif

                          z = drho_1(i,j,it) - drho_2(i,j,it);
                          Ar(i,j) = DREAL(z);
                          Ai(i,j) = DIMAG(z);
                      enddo
                  enddo
                  endif
              enddo
          enddo


          call diagUVtAUV(  N_total , NCOORD ,   k_PHI , N_blocks ,
     &                           Ar ,    NTX ,      Ai ,      NTX ,
     &                       dh_nnz ,    NBX ,
     &                       ia_spx , id_spx ,
     &                        PHI_U ,    NTX , PHI_SVt ,   KTRUNC , D ,
     &                         AUX1 ,    NTX ,    AUX2 ,   KTRUNC     );

          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  dj_z(ih,il,it) = D(ihl) / (2.D0*pi);
              enddo
          enddo

      enddo






c-----Calculation of ldrho_v (protons only)
      ldrho_v = COMPLEX( 0.D0 , 0.D0 );
      do it = 2 , 2






          do il = 1 , NGL
              rr = rb(il);
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  zz = zb(abs(ih));

                  fac = - DBLE(K_multipole*K_multipole)/(rr*rr);
                  fac = fac + 2.D0*rr*rr/(b0*b0*b0*b0*bp*bp*bp*bp);
                  fac = fac + 2.D0*zz*zz/(b0*b0*b0*b0*bz*bz*bz*bz);

                  ldrho_v(ih,il,it) = fac*drho_v(ih,il,it);
              enddo
          enddo






          Ar = 0.D0;
          Ai = 0.D0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh_nnz(ib1,ib2) ) then
                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);
                  do j = j0 , j0-1+id_spx(ib2)
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);
                      do i = i0 , i0-1+id_spx(ib1)
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1.ne.fg2  .or.
     &                        abs(ml1-ml2).ne.K_multipole ) then
                              CYCLE;
                          endif

                          z = drho_1(i,j,it) + drho_2(i,j,it);
                          Ar(i,j) = DREAL(z);
                          Ai(i,j) = DIMAG(z);
                      enddo
                  enddo
                  endif
              enddo
          enddo

          call diagUVtAUV( N_total , NCOORD ,  k_dzPHI , N_blocks ,
     &                          Ar ,    NTX ,       Ai ,      NTX ,
     &                      dh_nnz ,    NBX ,
     &                      ia_spx , id_spx ,
     &                     dzPHI_U ,    NTX , dzPHI_SVt,   KTRUNC , D ,
     &                        AUX1 ,    NTX ,     AUX2 ,   KTRUNC     );
          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  ldrho_v(ih,il,it) = ldrho_v(ih,il,it) + D(ihl)/pi;
              enddo
          enddo

          call diagUVtAUV( N_total , NCOORD ,  k_drPHI , N_blocks ,
     &                          Ar ,    NTX ,       Ai ,      NTX ,
     &                      dh_nnz ,    NBX ,
     &                      ia_spx , id_spx ,
     &                     drPHI_U ,    NTX , drPHI_SVt,   KTRUNC , D ,
     &                        AUX1 ,    NTX ,     AUX2 ,   KTRUNC     );
          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  ldrho_v(ih,il,it) = ldrho_v(ih,il,it) + D(ihl)/pi;
              enddo
          enddo






          Ar = 0.D0;
          Ai = 0.D0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh_nnz(ib1,ib2) ) then
                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);
                  do j = j0 , j0-1+id_spx(ib2)
                      fg2 = fg_spx(j-j0+1,ib2);
                      nz2 = nz_spx(j-j0+1,ib2);
                      nr2 = nr_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);
                      do i = i0 , i0-1+id_spx(ib1)
                          fg1 = fg_spx(i-i0+1,ib1);
                          nz1 = nz_spx(i-i0+1,ib1);
                          nr1 = nr_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1.ne.fg2  .or.
     &                        abs(ml1-ml2).ne.K_multipole ) then
                              CYCLE;
                          endif

                          z = drho_1(i,j,it) + drho_2(i,j,it);

                          fac = DBLE(abs(ml1)+abs(ml2))/(b0*b0*bp*bp);
                          fac = fac + DBLE(nz1+nz2+1)/(b0*b0*bz*bz);
                          fac = fac + DBLE(2*(nr1+nr2+1))/(b0*b0*bp*bp);

                          Ar(i,j) = fac * DREAL(z);
                          Ai(i,j) = fac * DIMAG(z);
                      enddo
                  enddo
                  endif
              enddo
          enddo

          call diagUVtAUV( N_total , NCOORD ,    k_PHI , N_blocks ,
     &                          Ar ,    NTX ,       Ai ,      NTX ,
     &                      dh_nnz ,    NBX ,
     &                      ia_spx , id_spx ,
     &                       PHI_U ,    NTX ,   PHI_SVt,   KTRUNC , D ,
     &                        AUX1 ,    NTX ,     AUX2 ,   KTRUNC     );
          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  ldrho_v(ih,il,it) = ldrho_v(ih,il,it) + D(ihl)/(-pi);
              enddo
          enddo






          Ar = 0.D0;
          Ai = 0.D0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh_nnz(ib1,ib2) ) then
                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);
                  do j = j0 , j0-1+id_spx(ib2)
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);
                      do i = i0 , i0-1+id_spx(ib1)
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1.ne.fg2  .or.
     &                        abs(ml1-ml2).ne.K_multipole ) then
                              CYCLE;
                          endif

                          z = drho_1(i,j,it) + drho_2(i,j,it);

                          fac = 0.5D0 * DBLE( ml1*ml1 + ml2*ml2 );

                          Ar(i,j) = fac * DREAL(z);
                          Ai(i,j) = fac * DIMAG(z);
                      enddo
                  enddo
                  endif
              enddo
          enddo

          call diagUVtAUV( N_total , NCOORD ,   k_PHIr , N_blocks ,
     &                          Ar ,    NTX ,       Ai ,      NTX ,
     &                      dh_nnz ,    NBX ,
     &                      ia_spx , id_spx ,
     &                      PHIr_U ,    NTX ,  PHIr_SVt,   KTRUNC , D ,
     &                        AUX1 ,    NTX ,     AUX2 ,   KTRUNC     );
          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  ldrho_v(ih,il,it) = ldrho_v(ih,il,it) + D(ihl)/pi;
              enddo
          enddo






      enddo






c-----Calculation of ldrho_s (protons + neutrons)
      ldrho_s = COMPLEX( 0.D0 , 0.D0 );
      do it = 3 , 3






          do il = 1 , NGL
              rr = rb(il);
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  zz = zb(abs(ih));

                  fac = - DBLE(K_multipole*K_multipole)/(rr*rr);
                  fac = fac + 2.D0*rr*rr/(b0*b0*b0*b0*bp*bp*bp*bp);
                  fac = fac + 2.D0*zz*zz/(b0*b0*b0*b0*bz*bz*bz*bz);

                  ldrho_s(ih,il) = fac*drho_s(ih,il);
              enddo
          enddo






          Ar = 0.D0;
          Ai = 0.D0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh_nnz(ib1,ib2) ) then
                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);
                  do j = j0 , j0-1+id_spx(ib2)
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);
                      do i = i0 , i0-1+id_spx(ib1)
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1.ne.fg2  .or.
     &                        abs(ml1-ml2).ne.K_multipole ) then
                              CYCLE;
                          endif

                          z = + drho_1(i,j,1) + drho_2(i,j,1)
     &                        + drho_1(i,j,2) + drho_2(i,j,2);
                          if( fg1.eq.'g' .and. fg2.eq.'g' ) then
                              z = -z;
                          endif
                          Ar(i,j) = DREAL(z);
                          Ai(i,j) = DIMAG(z);
                      enddo
                  enddo
                  endif
              enddo
          enddo

          call diagUVtAUV( N_total , NCOORD ,  k_dzPHI , N_blocks ,
     &                          Ar ,    NTX ,       Ai ,      NTX ,
     &                      dh_nnz ,    NBX ,
     &                      ia_spx , id_spx ,
     &                     dzPHI_U ,    NTX , dzPHI_SVt,   KTRUNC , D ,
     &                        AUX1 ,    NTX ,     AUX2 ,   KTRUNC     );
          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  ldrho_s(ih,il) = ldrho_s(ih,il) + D(ihl)/pi;
              enddo
          enddo

          call diagUVtAUV( N_total , NCOORD ,  k_drPHI , N_blocks ,
     &                          Ar ,    NTX ,       Ai ,      NTX ,
     &                      dh_nnz ,    NBX ,
     &                      ia_spx , id_spx ,
     &                     drPHI_U ,    NTX , drPHI_SVt,   KTRUNC , D ,
     &                        AUX1 ,    NTX ,     AUX2 ,   KTRUNC     );
          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  ldrho_s(ih,il) = ldrho_s(ih,il) + D(ihl)/pi;
              enddo
          enddo






          Ar = 0.D0;
          Ai = 0.D0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh_nnz(ib1,ib2) ) then
                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);
                  do j = j0 , j0-1+id_spx(ib2)
                      fg2 = fg_spx(j-j0+1,ib2);
                      nz2 = nz_spx(j-j0+1,ib2);
                      nr2 = nr_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);
                      do i = i0 , i0-1+id_spx(ib1)
                          fg1 = fg_spx(i-i0+1,ib1);
                          nz1 = nz_spx(i-i0+1,ib1);
                          nr1 = nr_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1.ne.fg2  .or.
     &                        abs(ml1-ml2).ne.K_multipole ) then
                              CYCLE;
                          endif

                          z = + drho_1(i,j,1) + drho_2(i,j,1)
     &                        + drho_1(i,j,2) + drho_2(i,j,2);
                          if( fg1.eq.'g' .and. fg2.eq.'g' ) then
                              z = -z;
                          endif

                          fac = DBLE(abs(ml1)+abs(ml2))/(b0*b0*bp*bp);
                          fac = fac + DBLE(nz1+nz2+1)/(b0*b0*bz*bz);
                          fac = fac + DBLE(2*(nr1+nr2+1))/(b0*b0*bp*bp);

                          Ar(i,j) = fac * DREAL(z);
                          Ai(i,j) = fac * DIMAG(z);
                      enddo
                  enddo
                  endif
              enddo
          enddo

          call diagUVtAUV( N_total , NCOORD ,    k_PHI , N_blocks ,
     &                          Ar ,    NTX ,       Ai ,      NTX ,
     &                      dh_nnz ,    NBX ,
     &                      ia_spx , id_spx ,
     &                       PHI_U ,    NTX ,   PHI_SVt,   KTRUNC , D ,
     &                        AUX1 ,    NTX ,     AUX2 ,   KTRUNC     );
          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  ldrho_s(ih,il) = ldrho_s(ih,il) + D(ihl)/(-pi);
              enddo
          enddo






          Ar = 0.D0;
          Ai = 0.D0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh_nnz(ib1,ib2) ) then
                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);
                  do j = j0 , j0-1+id_spx(ib2)
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);
                      do i = i0 , i0-1+id_spx(ib1)
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1.ne.fg2  .or.
     &                        abs(ml1-ml2).ne.K_multipole ) then
                              CYCLE;
                          endif

                          z = + drho_1(i,j,1) + drho_2(i,j,1)
     &                        + drho_1(i,j,2) + drho_2(i,j,2);
                          if( fg1.eq.'g' .and. fg2.eq.'g' ) then
                              z = -z;
                          endif

                          fac = 0.5D0 * DBLE( ml1*ml1 + ml2*ml2 );

                          Ar(i,j) = fac * DREAL(z);
                          Ai(i,j) = fac * DIMAG(z);
                      enddo
                  enddo
                  endif
              enddo
          enddo

          call diagUVtAUV( N_total , NCOORD ,   k_PHIr , N_blocks ,
     &                          Ar ,    NTX ,       Ai ,      NTX ,
     &                      dh_nnz ,    NBX ,
     &                      ia_spx , id_spx ,
     &                      PHIr_U ,    NTX ,  PHIr_SVt,   KTRUNC , D ,
     &                        AUX1 ,    NTX ,     AUX2 ,   KTRUNC     );
          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  ldrho_s(ih,il) = ldrho_s(ih,il) + D(ihl)/pi;
              enddo
          enddo






      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_ddensdcurr() *************************';
      write(6,*) '';
      endif

      return;
      end;






c======================================================================c

      subroutine diagUVtAUV(    N ,      M ,      K ,      NB ,
     &                         AR ,   LDAR ,
     &                         AI ,   LDAI ,
     &                       Annz , LDAnnz ,
     &                     ia_blk , id_blk ,
     &                          U ,    LDU ,
     &                          V ,    LDV ,      D ,
     &                       AUX1 , LDAUX1 ,
     &                       AUX2 , LDAUX2                      )

c======================================================================c
      ! Calculates the diagonal D = diag( (U*V)^T * A * (U*V) )
      ! of a complex matrix A whose block structure is given by Annz
      !
      !  [in]: A         is complex  NxN  matrix, A = AR + j*AI
      !        Annz      is logical NBxNB matrix
      !        U         is real     NxK  matrix
      !        V         is real     KxM  matrix
      !        AUX1      is real     NxK  matrix (auxiliary)
      !        AUX2      is real     KxK  matrix (auxiliary)
      !        ia_blk(i) is the first index of i-th block
      !        id_blk(i) is the size of i-th block
      !
      !
      ! [out]: D is complex vector od length M

      IMPLICIT NONE;
      EXTERNAL ddot;
      REAL*8   ddot;

      INTEGER*4 N, M, K, NB;
      INTEGER*4 LDAR, LDAI, LDAnnz, LDU, LDV, LDAUX1, LDAUX2;
      INTEGER*4 ia_blk( NB );
      INTEGER*4 id_blk( NB );
      REAL*8 AR( LDAR , N );
      REAL*8 AI( LDAI , N );
      REAL*8  U(  LDU , K );
      REAL*8  V(  LDV , M );
      REAL*8 AUX1( LDAUX1 , K );
      REAL*8 AUX2( LDAUX2 , K );
      LOGICAL Annz( LDAnnz , NB );
      COMPLEX*16 D( M );



      REAL*8 Y( LDV );
      INTEGER*4 i, ib, jb, i0, j0;

      D = COMPLEX( 0.D0 , 0.D0 );



c-----Real part
      ! Calculation of AUX1 = AR*U
      AUX1 = 0.D0;
      do ib = 1 , NB
          do jb = 1 , NB
              i0 = ia_blk(ib);
              j0 = ia_blk(jb);
              if( Annz(ib,jb) ) then
                  call dgemm(        'N' ,   'N' ,
     &                        id_blk(ib) ,     K , id_blk(jb) ,
     &                              1.D0 ,
     &                         AR(i0,j0) ,  LDAR ,
     &                           U(j0,1) ,   LDU ,
     &                              1.D0 ,
     &                        AUX1(i0,1) , LDAUX1  );

              endif
          enddo
      enddo

      ! Calculation of AUX2 = U^T * AUX1 = U^T * AR * U
      call dgemm(    'T' ,    'N' ,
     &                 K ,      K ,      N ,
     &              1.D0 ,
     &                 U ,    LDU ,
     &              AUX1 , LDAUX1 ,
     &              0.D0 ,
     &              AUX2 , LDAUX2   );

      ! Calculation of < V_.i | AUX2 * V_.i >
      do i = 1 , M
          ! Calculation of Y = AUX2 * V_.i
          call dgemv(  'N' ,
     &                   K ,       K ,
     &                1.D0 ,
     &                AUX2 ,  LDAUX2 ,
     &              V(1,i) ,       1 ,
     &                0.D0 ,
     &                   Y ,       1  );

          ! Calculation of < V_.i | Y >
          D(i) = D(i) + COMPLEX( ddot(K,V(1,i),1,Y,1) , 0.D0 );
      enddo






c-----Imaginary part
      ! Calculation of AUX1 = AI*U
      AUX1 = 0.D0;
      do ib = 1 , NB
          do jb = 1 , NB
              i0 = ia_blk(ib);
              j0 = ia_blk(jb);
              if( Annz(ib,jb) ) then
                  call dgemm(        'N' ,   'N' ,
     &                        id_blk(ib) ,     K , id_blk(jb) ,
     &                              1.D0 ,
     &                         AI(i0,j0) ,  LDAI ,
     &                           U(j0,1) ,   LDU ,
     &                              1.D0 ,
     &                        AUX1(i0,1) , LDAUX1  );

              endif
          enddo
      enddo

      ! Calculation of AUX2 = U^T * AUX1 = U^T * AI * U
      call dgemm(    'T' ,    'N' ,
     &                 K ,      K ,      N ,
     &              1.D0 ,
     &                 U ,    LDU ,
     &              AUX1 , LDAUX1 ,
     &              0.D0 ,
     &              AUX2 , LDAUX2   );

      ! Calculation of < V_.i | AUX2 * V_.i >
      do i = 1 , M
          ! Calculation of Y = AUX2 * V_.i
          call dgemv(  'N' ,
     &                   K ,       K ,
     &                1.D0 ,
     &                AUX2 ,  LDAUX2 ,
     &              V(1,i) ,       1 ,
     &                0.D0 ,
     &                   Y ,       1  );

          ! Calculation of < V_.i | Y >
          D(i) = D(i) + COMPLEX( 0.D0 , ddot(K,V(1,i),1,Y,1) );
      enddo





      return;
      end;
