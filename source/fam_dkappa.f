c======================================================================c

      subroutine fam_dkappa( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE simplex;
      USE u_matrix;
      USE v_matrix;
      USE nnz_blocks;
      USE dkappa;
      USE xyfam;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      DOUBLE COMPLEX Temp( NTX , NTX );



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






c-----Calculation of dkappa_mi = - v * x_fam * hermconj(v)
c-----                           + u * y_fam * hermconj(u)
      dkappa_mi = COMPLEX( 0.D0 , 0.D0 );
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
     &                        v(1,1,ib1,it)             , NBSX        ,
     &                        x_fam(i0,j0,it)           , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         );

                  call zgemm( 'n'         , 'c'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Temp(i0,j0)               , NTX         ,
     &                        v(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        dkappa_mi(i0,j0,it)       , NTX         );



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
