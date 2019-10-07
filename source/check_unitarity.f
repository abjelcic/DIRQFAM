c======================================================================c

      subroutine check_unitarity( lpr , it )

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



      COMPLEX*16 Temp( NBSX , NBSX );



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN check_unitarity() **********************';
      write(6,*) '';
      write(6,*) '---------------------------------------';
      if( it .eq. 1 ) then
          write(6,*) '-------------- Neutrons: --------------';
      else
          write(6,*) '-------------- Protons: ---------------';
      endif
      write(6,*) '---------------------------------------';
      write(6,*) '';
      endif






c----- hermconj(U) * U  +  hermconj(V) * V = I <==>
c----- hermconj(u) * u  +  hermconj(v) * v = I
      error = 0.D0;
      do ib = 1 , N_blocks
          call zgemm( 'c'        ,   'n'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( +1.D0 , 0.D0 )   ,
     &                u(1,1,ib,it)              , NBSX       ,
     &                u(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( 0.D0 , 0.D0 )    ,
     &                Temp(1,1)                 , NBSX         );

          call zgemm( 'c'        ,   'n'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( +1.D0 , 0.D0 )   ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( +1.D0 , 0.D0 )   ,
     &                Temp(1,1)                 , NBSX         );

          do i = 1 , id_spx(ib)
              Temp(i,i) = Temp(i,i) - COMPLEX( +1.D0 , 0.D0 );
              do j = 1 , id_spx(ib)
                  error = MAX( error , ABS(Temp(i,j)) );
              enddo
          enddo

      enddo
      call assert( error.lt.1.D-8 , 'Bogoliubov transf. not unitary' );






c----- U * hermconj(U)  +  V * hermconj(V) = I <==>
c----- u * hermconj(u)  +  v * hermconj(v) = I
      error = 0.D0;
      do ib = 1 , N_blocks
          call zgemm( 'n'        ,   'c'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( +1.D0 , 0.D0 )   ,
     &                u(1,1,ib,it)              , NBSX       ,
     &                u(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( 0.D0 , 0.D0 )    ,
     &                Temp(1,1)                 , NBSX         );

          call zgemm( 'n'        ,   'c'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( +1.D0 , 0.D0 )   ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( +1.D0 , 0.D0 )   ,
     &                Temp(1,1)                 , NBSX         );

          do i = 1 , id_spx(ib)
              Temp(i,i) = Temp(i,i) - COMPLEX( +1.D0 , 0.D0 );
              do j = 1 , id_spx(ib)
                  error = MAX( error , ABS(Temp(i,j)) );
              enddo
          enddo

      enddo
      call assert( error.lt.1.D-8 , 'Bogoliubov transf. not unitary' );






c-----   transp(U) * V  +    transp(V) * U = 0 <==>
c----- hermconj(u) * v  -  hermconj(v) * u = 0
      error = 0.D0;
      do ib = 1 , N_blocks
          call zgemm( 'c'        ,   'n'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( +1.D0 , 0.D0 )   ,
     &                u(1,1,ib,it)              , NBSX       ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( 0.D0 , 0.D0 )    ,
     &                Temp(1,1)                 , NBSX         );

          call zgemm( 'c'        ,   'n'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( +1.D0 , 0.D0 )   ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                u(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( -1.D0 , 0.D0 )   ,
     &                Temp(1,1)                 , NBSX         );

          do i = 1 , id_spx(ib)
              do j = 1 , id_spx(ib)
                  error = MAX( error , ABS(Temp(i,j)) );
              enddo
          enddo

      enddo
      call assert( error.lt.1.D-8 , 'Bogoliubov transf. not unitary' );






c----- U * hermconj(V)  +  conj(V) *   transp(U)   = 0 <==>
c----- u * hermconj(v)  -       v  * hermconj(u)   = 0
      error = 0.D0;
      do ib = 1 , N_blocks
          call zgemm( 'n'        ,   'c'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( +1.D0 , 0.D0 )   ,
     &                u(1,1,ib,it)              , NBSX       ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( 0.D0 , 0.D0 )    ,
     &                Temp(1,1)                 , NBSX         );

          call zgemm( 'n'        ,   'c'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( +1.D0 , 0.D0 )   ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                u(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( -1.D0 , 0.D0 )   ,
     &                Temp(1,1)                 , NBSX         );

          do i = 1 , id_spx(ib)
              do j = 1 , id_spx(ib)
                  error = MAX( error , ABS(Temp(i,j)) );
              enddo
          enddo

      enddo
      call assert( error.lt.1.D-8 , 'Bogoliubov transf. not unitary' );






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END check_unitarity() ************************';
      write(6,*) '';
      endif

      return;
      end;
