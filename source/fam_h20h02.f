c======================================================================c

      subroutine fam_h20h02( lpr )

c======================================================================c

      USE dirqfampar;
      USE simplex;
      USE u_matrix;
      USE v_matrix;
      USE nnz_blocks;
      USE dh;
      USE dDelta;
      USE h02h20matrix;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      DOUBLE COMPLEX Temp( NTX , NTX );



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
