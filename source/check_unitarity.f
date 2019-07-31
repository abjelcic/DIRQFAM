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






c-----adj(U)*U + ajd(V)*V = I <==> adj(u)*u + adj(v)*v = I
      error = -1.D0;
      do ib = 1 , N_blocks
          call zgemm( 'c'        ,   'n'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( 1.D0 , 0.D0 )    ,
     &                u(1,1,ib,it)              , NBSX       ,
     &                u(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( 0.D0 , 0.D0 )    ,
     &                Temp(1,1)                 , NBSX         );

          call zgemm( 'c'        ,   'n'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( 1.D0 , 0.D0 )    ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( 1.D0 , 0.D0 )    ,
     &                Temp(1,1)                 , NBSX         );

          do i = 1 , id_spx(ib)
              Temp(i,i) = Temp(i,i) - COMPLEX( 1.D0 , 0.D0 );
              do j = 1 , id_spx(ib)
                  error = max( error , ABS(Temp(i,j)) );
              enddo
          enddo

      enddo
      if( error .gt. 1.D-8 ) then
          write(6,*)'|| ajd(U)*U + adj(V)*V - I ||_max =       ', error;
          stop 'Error: Bogoliubov transformation not unitary!';
      endif
      if( lpr ) then
          write(6,*)'|| ajd(U)*U + adj(V)*V - I ||_max =       ', error;
      endif






c-----U*adj(U) + V*adj(V) = I <==> u*adj(u) + v*adj(v) = I
      error = -1.D0;
      do ib = 1 , N_blocks
          call zgemm( 'n'        ,   'c'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( 1.D0 , 0.D0 )    ,
     &                u(1,1,ib,it)              , NBSX       ,
     &                u(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( 0.D0 , 0.D0 )    ,
     &                Temp(1,1)                 , NBSX         );

          call zgemm( 'n'        ,   'c'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( 1.D0 , 0.D0 )    ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( 1.D0 , 0.D0 )    ,
     &                Temp(1,1)                 , NBSX         );

          do i = 1 , id_spx(ib)
              Temp(i,i) = Temp(i,i) - COMPLEX( 1.D0 , 0.D0 );
              do j = 1 , id_spx(ib)
                  error = max( error , ABS(Temp(i,j)) );
              enddo
          enddo

      enddo
      if( error .gt. 1.D-8 ) then
          write(6,*)'|| U*ajd(U) + V*adj(V) - I ||_max =       ', error;
          stop 'Error: Bogoliubov transformation not unitary!';
      endif
      if( lpr ) then
          write(6,*)'|| U*ajd(U) + V*adj(V) - I ||_max =       ', error;
      endif






c-----transp(U)*V + transp(V)*U = 0 <==> adj(u)*v - adj(v)*u = 0
      error = -1.D0;
      do ib = 1 , N_blocks
          call zgemm( 'c'        ,   'n'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( 1.D0 , 0.D0 )    ,
     &                u(1,1,ib,it)              , NBSX       ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( 0.D0 , 0.D0 )    ,
     &                Temp(1,1)                 , NBSX         );

          call zgemm( 'c'        ,   'n'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( 1.D0 , 0.D0 )    ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                u(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( -1.D0 , 0.D0 )   ,
     &                Temp(1,1)                 , NBSX         );

          do i = 1 , id_spx(ib)
              do j = 1 , id_spx(ib)
                  error = max( error , ABS(Temp(i,j)) );
              enddo
          enddo

      enddo
      if( error .gt. 1.D-8 ) then
          write(6,*)'|| transp(U)*V + transp(V)*U ||_max =     ', error;
          stop 'Error: Bogoliubov transformation not unitary!';
      endif
      if( lpr ) then
          write(6,*)'|| transp(U)*V + transp(V)*U ||_max =     ', error;
      endif






c-----U*adj(V) + conj(V)*transp(U) = 0 <==> u*adj(v) - v*ajd(u) = 0
      error = -1.D0;
      do ib = 1 , N_blocks
          call zgemm( 'n'        ,   'c'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( 1.D0 , 0.D0 )    ,
     &                u(1,1,ib,it)              , NBSX       ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( 0.D0 , 0.D0 )    ,
     &                Temp(1,1)                 , NBSX         );

          call zgemm( 'n'        ,   'c'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( 1.D0 , 0.D0 )    ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                u(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( -1.D0 , 0.D0 )   ,
     &                Temp(1,1)                 , NBSX         );

          do i = 1 , id_spx(ib)
              do j = 1 , id_spx(ib)
                  error = max( error , ABS(Temp(i,j)) );
              enddo
          enddo

      enddo
      if( error .gt. 1.D-8 ) then
          write(6,*)'|| U*adj(V) + conj(V)*transp(U) ||_max =  ', error;
          stop 'Error: Bogoliubov transformation not unitary!';
      endif
      if( lpr ) then
          write(6,*)'|| U*adj(V) + conj(V)*transp(U) ||_max =  ', error;
      endif






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END check_unitarity() ************************';
      write(6,*) '';
      endif

      return;
      end;
