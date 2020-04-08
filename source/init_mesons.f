c======================================================================c

      subroutine init_mesons( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE ddpc1ddme2;
      USE quadrature;
      USE mesmat;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;

      CHARACTER parname*10;
      common /partyp/ parname;
      common /baspar/ hom, hb0, b0;
      common /defbas/ beta0, q, bp, bz;
      common /basnnn/ n0f, n0b;



      INTEGER nz_mes( NHMAX , 0:J_MAX+1 );
      INTEGER nr_mes( NHMAX , 0:J_MAX+1 );

      DOUBLE PRECISION H   ( NHMAX , NHMAX , 0:J_MAX+1 );
      DOUBLE PRECISION Htmp( NHMAX , NHMAX             );

      DOUBLE PRECISION phiz ( -NGH:NGH , 0:(NMESMAX)               );
      DOUBLE PRECISION phirK(    1:NGL , 0:(NMESMAX/2) , 0:J_MAX+1 );

      hbc    = 197.328284D0;
      b0_mes = b0 / DSQRT(2.D0);
      bp_mes = bp;
      bz_mes = bz;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN init_mesons() *************************';
      write(6,*) '';
      endif


      if( parname .ne. 'DD-ME2' ) then
          return;
      endif






c-----Setting dimensions of P and H matrices
      do K = 0 , J_MAX+1
          N = n0b - K;
          nP(K) = ( (N+1)*(N+3) + 1 - MOD(N,2) )/4;

          call assert( nP(K).le.NHMAX , 'NHMAX too small' );
      enddo






c-----Construction of (nz,nr) pairs with nz + 2*nr + K <= NMESMAX
      do K = 0 , J_MAX+1
          il = 0;
          do nr = 0 , ( NMESMAX - K )/2
              do nz = 0 , NMESMAX - K - 2*nr
                  il = il + 1;
                  nz_mes(il,K) = nz;
                  nr_mes(il,K) = nr;
              enddo
          enddo
          call assert( il.eq.nP(K) , 'nP(K) wrong' );
      enddo






c-----Construction of H matrices
      H = 0.D0;
      do K = 0 , J_MAX+1
          do i = 1 , nP(K)
              nz1 = nz_mes(i,K);
              nr1 = nr_mes(i,K);
              do j = 1 , nP(K)
                  nz2 = nz_mes(j,K);
                  nr2 = nr_mes(j,K);

                  facp = 1.D0 / ( b0_mes*bp_mes )**2.D0
                  facz = 1.D0 / ( b0_mes*bz_mes )**2.D0;

                  if( nz1.eq.nz2 .and. nr1.eq.nr2 ) then
                      x = (DBLE(nz1)+0.5D0)*facz + DBLE(2*nr1+K+1)*facp;
                      H(i,j,K) = H(i,j,K) + x;
                  endif

                  if( nz1.eq.nz2 .and. nr1.eq.nr2+1 ) then
                      x = DSQRT( DBLE(nr1*(nr1+K)) ) * facp;
                      H(i,j,K) = H(i,j,K) + x;
                  endif

                  if( nz1.eq.nz2 .and. nr2.eq.nr1+1 ) then
                      x = DSQRT( DBLE(nr2*(nr2+K)) ) * facp;
                      H(i,j,K) = H(i,j,K) + x;
                  endif

                  if( nz1.eq.nz2+2 .and. nr1.eq.nr2 ) then
                      x = - DSQRT(DBLE((nz2+1)*(nz2+2)))/2.D0 * facz;
                      H(i,j,K) = H(i,j,K) + x;
                  endif

                  if( nz2.eq.nz1+2 .and. nr1.eq.nr2 ) then
                      x = - DSQRT(DBLE((nz1+1)*(nz1+2)))/2.D0 * facz;
                      H(i,j,K) = H(i,j,K) + x;
                  endif

              enddo
          enddo
      enddo






c-----Calculation of phiz and phirK
      do nz = 0 , NMESMAX
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              z = DBLE(isign(1,ih)) * zb_fam(abs(ih));
              phiz(ih,nz) = phi_nz(nz,b0_mes*bz_mes,z);

          enddo
      enddo
      do K = 0 , J_MAX+1
          do nr = 0 , NMESMAX/2
              do il = 1 , NGL
                  r = rb_fam(il);
                  phirK(il,nr,K) = phi_nr_ml(nr,K,b0_mes*bp_mes,r);
              enddo
          enddo
      enddo






c-----Calculation of {Psig,Pome,Prho} matrices
      do K = 0 , J_MAX+1


          do i = 1 , nP(K)
              nz = nz_mes(i,K);
              nr = nr_mes(i,K);
              call assert( nz.le.NMESMAX   , 'nz > NMESMAX'   );
              call assert( nr.le.NMESMAX/2 , 'nr > NMESMAX/2' );

              ihl = 0;
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih .eq. 0 ) CYCLE;

                      ihl = ihl + 1;
                      Psig(i,ihl,K) = phiz(ih,nz) * phirK(il,nr,K);
                      Pome(i,ihl,K) = phiz(ih,nz) * phirK(il,nr,K);
                      Prho(i,ihl,K) = phiz(ih,nz) * phirK(il,nr,K);

                  enddo
              enddo
          enddo


          do imes = 1 , 3

              do i = 1 , nP(K)
                  do j = 1 , nP(K)
                      Htmp(i,j) = H(i,j,K);
                  enddo
              enddo

              select case( imes )
                  case( 1 )
                      do i = 1 , nP(K)
                          Htmp(i,i) = Htmp(i,i) + (m_sig/hbc)**2.D0;
                      enddo
                  case( 2 )
                      do i = 1 , nP(K)
                          Htmp(i,i) = Htmp(i,i) + (m_ome/hbc)**2.D0;
                      enddo
                  case( 3 )
                      do i = 1 , nP(K)
                          Htmp(i,i) = Htmp(i,i) + (m_rho/hbc)**2.D0;
                      enddo
                  case default
                      stop 'Error: imes > 3!';
              end select

              call dpotrf( 'L' , nP(K) , Htmp , NHMAX , INFO );
              call assert( INFO.eq.0 , 'dpotrf() failed' );

              select case( imes )
                  case( 1 )
                      call dtrsm( 'L'            , 'L'    ,
     &                            'N'            , 'N'    ,
     &                            nP(K)          , NCOORD , 1.D0 ,
     &                            Htmp           , NHMAX  ,
     &                            Psig(1,1,K)    , NHMAX           );
                  case( 2 )
                      call dtrsm( 'L'            , 'L'    ,
     &                            'N'            , 'N'    ,
     &                            nP(K)          , NCOORD , 1.D0 ,
     &                            Htmp           , NHMAX  ,
     &                            Pome(1,1,K)    , NHMAX           );
                  case( 3 )
                      call dtrsm( 'L'            , 'L'    ,
     &                            'N'            , 'N'    ,
     &                            nP(K)          , NCOORD , 1.D0 ,
     &                            Htmp           , NHMAX  ,
     &                            Prho(1,1,K)    , NHMAX           );
                  case default
                      stop 'Error: imes > 3!';
              end select

          enddo


      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END init_mesons() ***************************';
      write(6,*) '';
      endif

      return;
      end;
