c======================================================================c

      subroutine fam_h20h02( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

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

      COMPLEX*16 dh_1, dh_2;
      common /delta_h/ dh_1( NTX , NTX , 2 ),
     &                 dh_2( NTX , NTX , 2 );

      COMPLEX*16 dDelta_pl, dDelta_mi;
      common /dDelta/ dDelta_pl( NTX , NTX , 2 ),
     &                dDelta_mi( NTX , NTX , 2 );

      COMPLEX*16 h20, h02;
      common /h20h02/ h20( NTX , NTX , 2 ),
     &                h02( NTX , NTX , 2 );



      COMPLEX*16  Temp( NTX , NTX );



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_h20h02() **************************';
      write(6,*) '';
      endif


c----------------------------------------------------------------------c
c     NOTE:                                                            c
c     In code some matrices are altered by addition of a sign, or by   c
c     transposition/Hermitian conjugation to make some formulae        c
c     simpler, and to avoid many unnecessary transpositions of         c
c     matrices in the code.                                            c
c     All transformations are idempotent and therefore applying it     c
c     twice yields the original matrix.                                c
c                                                                      c
c                                                                      c
c     Matrix in paper:                  Matrix in code:                c
c                                                                      c
c     f_1                               + f1_JK                        c
c     f_2                               + f2_JK                        c
c     \delta \Delta^{(-)}(\omega)       - hermconj(dDelta_mi)          c
c     \delta \kappa^{(-)}(\omega)       - hermconj(dkappa_mi)          c
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
c     Other matrices have consistent definition with ones in the paper c
c----------------------------------------------------------------------c






c-----Calculation of h20 = + hermconj(u) * dh_1         * v
c-----                     + hermconj(v) * transp(dh_2) * u
c-----                     - hermconj(u) * dDelta_pl    * u
c-----                     - hermconj(v) * dDelta_mi    * v
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



                  call zgemm( 'c'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( -1.D0 , 0.D0 )   ,
     &                        v(1,1,ib1,it)             , NBSX        ,
     &                        dDelta_mi(i0,j0,it)       , NTX         ,
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






c-----Calculation of h02 = + hermconj(v) * dh_1         * u
c-----                     + hermconj(u) * transp(dh_2) * v
c-----                     + hermconj(v) * dDelta_pl    * v
c-----                     + hermconj(u) * dDelta_mi    * u
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



                  call zgemm( 'c'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        u(1,1,ib1,it)             , NBSX        ,
     &                        dDelta_mi(i0,j0,it)       , NTX         ,
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






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_h20h02() ****************************';
      write(6,*) '';
      endif

      return;
      end;
