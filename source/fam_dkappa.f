c======================================================================c

      subroutine fam_dkappa( lpr )

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

      COMPLEX*16 v;
      common /v_matrix/ v( NBSX , NBSX , NBX , 2 );

      COMPLEX*16 u;
      common /u_matrix/ u( NBSX , NBSX , NBX , 2 );

      LOGICAL dh_nnz, dDelta_nnz, dkappa_nnz, f_nnz;
      common /nnz_blocks/ dh_nnz    ( NBX , NBX ),
     &                    dDelta_nnz( NBX , NBX ),
     &                    dkappa_nnz( NBX , NBX ),
     &                    f_nnz     ( NBX , NBX );

      COMPLEX*16 dkappa_pl, dkappa_mi;
      common /delta_kappa/ dkappa_pl( NTX , NTX , 2 ),
     &                     dkappa_mi( NTX , NTX , 2 );

      COMPLEX*16 x_fam, y_fam;
      common /xy_fam/ x_fam( NTX , NTX , 2 ),
     &                y_fam( NTX , NTX , 2 );



      COMPLEX*16 Temp( NTX , NTX );



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_dkappa() **************************';
      write(6,*) '';
      endif


      if( i_pairing .eq. 0 ) then
          dkappa_pl = COMPLEX( 0.D0 , 0.D0 );
          dkappa_mi = COMPLEX( 0.D0 , 0.D0 );
          return;
      endif


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
      write(6,*) '****** END fam_dkappa() ****************************';
      write(6,*) '';
      endif

      return;
      end;
