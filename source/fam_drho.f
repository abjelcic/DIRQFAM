c======================================================================c

      subroutine fam_drho( lpr )

c======================================================================c

      USE dirqfampar;
      USE simplex;
      USE u_matrix;
      USE v_matrix;
      USE drho;
      USE xyfam;
      USE nnz_blocks;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      DOUBLE COMPLEX z_tmp;
      DOUBLE COMPLEX Temp( NTX , NTX );



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_drho() ****************************';
      write(6,*) '';
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






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_drho() ******************************';
      write(6,*) '';
      endif

      return;
      end;
