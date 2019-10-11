c======================================================================c

      subroutine init_mesons( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      CHARACTER parname*10;
      common /partyp/ parname;
      common /baspar/ hom, hb0, b0;
      common /defbas/ beta0, q, bp, bz;

      common /fam/ omega_start, omega_end, delta_omega, omega_print,
     &             omega, gamma_smear,
     &             i_calculation_type, i_coulomb, i_pairing,
     &             J_multipole, K_multipole, ISO;

      REAL*8 m_sig, m_ome, m_rho;
      common /DDPC1_DDME2/  a_s , b_s , c_s , d_s ,
     &                      a_v , b_v , c_v , d_v ,
     &                      a_tv, b_tv, c_tv, d_tv,
     &                      del_s,
     &
     &                      a_sig, b_sig, c_sig, d_sig, g0_sig, m_sig,
     &                      a_ome, b_ome, c_ome, d_ome, g0_ome, m_ome,
     &                      a_rho,                      g0_rho, m_rho,
     &
     &                      rho_sat;

      common /quadrature/ zb_fam( 1:NGH ), zb_famK( 1:NGH ),
     &                    rb_fam( 1:NGL ), rb_famK( 1:NGL ),
     &                    wzwr( 1:NGH , 1:NGL ), wzwrK( 1:NGH , 1:NGL );

      common /mesmat/ Qsig( NHSIZE , NCOORD ), QKsig( NHSIZE , NCOORD ),
     &                Qome( NHSIZE , NCOORD ), QKome( NHSIZE , NCOORD ),
     &                Qrho( NHSIZE , NCOORD ), QKrho( NHSIZE , NCOORD ),
     &                   P( NHSIZE , NCOORD );



      INTEGER*4 nz_mes( NHSIZE );
      INTEGER*4 nr_mes( NHSIZE );

      REAL*8 H   ( NHSIZE , NHSIZE );
      REAL*8 Hsig( NHSIZE , NHSIZE );
      REAL*8 Home( NHSIZE , NHSIZE );
      REAL*8 Hrho( NHSIZE , NHSIZE );

      REAL*8 phiz ( -NGH:NGH , 0:(NMESMAX)   );
      REAL*8 phirK(    1:NGL , 0:(NMESMAX/2) );

      hbc    = 197.328284D0;
      b0_mes = b0 / DSQRT(2.D0);
      bp_mes = bp;
      bz_mes = bz;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN init_mesons() **************************';
      write(6,*) '';
      endif


      if( parname .ne. 'DD-ME2' ) then
          return;
      endif






c-----Construction of (nz,nr) pairs with nz + 2*nr + K <= NMESMAX
      il = 0;
      do nr = 0 , ( NMESMAX - K_multipole )/2
          do nz = 0 , NMESMAX - K_multipole - 2*nr
              il = il + 1;
              nz_mes(il) = nz;
              nr_mes(il) = nr;
          enddo
      enddo
      call assert( il.eq.NHSIZE , 'NHSIZE too small' );






c-----Construction of H matrix
      H = 0.D0;
      do i = 1 , NHSIZE
          nz1 = nz_mes(i);
          nr1 = nr_mes(i);
          do j = 1 , NHSIZE
              nz2 = nz_mes(j);
              nr2 = nr_mes(j);

              facp = 1.D0 / ( b0_mes*bp_mes )**2.D0
              facz = 1.D0 / ( b0_mes*bz_mes )**2.D0;

              if( nz1.eq.nz2 .and. nr1.eq.nr2 ) then
                  x = + ( DBLE(nz1) + 0.5D0 )     * facz
     &                + DBLE(2*nr1+K_multipole+1) * facp;

                  H(i,j) = H(i,j) + x;
              endif

              if( nz1.eq.nz2 .and. nr1.eq.nr2+1 ) then
                  x = DSQRT( DBLE(nr1*(nr1+K_multipole)) ) * facp;

                  H(i,j) = H(i,j) + x;
              endif

              if( nz1.eq.nz2 .and. nr2.eq.nr1+1 ) then
                  x = DSQRT( DBLE(nr2*(nr2+K_multipole)) ) * facp;

                  H(i,j) = H(i,j) + x;
              endif

              if( nz1.eq.nz2+2 .and. nr1.eq.nr2 ) then
                  x = - DSQRT(DBLE((nz2+1)*(nz2+2)))/2.D0 * facz;

                  H(i,j) = H(i,j) + x;
              endif

              if( nz2.eq.nz1+2 .and. nr1.eq.nr2 ) then
                  x = - DSQRT(DBLE((nz1+1)*(nz1+2)))/2.D0 * facz;

                  H(i,j) = H(i,j) + x;
              endif

          enddo
      enddo






c-----Construction of P matrix
      do nz = 0 , NMESMAX
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;
              z = DBLE(isign(1,ih)) * zb_fam(abs(ih));
              phiz(ih,nz) = phi_nz( nz , b0_mes*bz_mes , z );
          enddo
      enddo
      do nr = 0 , NMESMAX/2
          do il = 1 , NGL
              r = rb_fam(il);
              K = K_multipole;
              phirK(il,nr) = phi_nr_ml( nr , K , b0_mes*bp_mes , r );
          enddo
      enddo
      do i = 1 , NHSIZE
          nz = nz_mes(i);
          nr = nr_mes(i);

          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  ihl = ihl + 1;
                  P( i , ihl ) = phiz(ih,nz) * phirK(il,nr);

              enddo
          enddo

      enddo






c-----Calculation of {Qsig,Qome,Qrho} matrices
      do nz = 0 , NMESMAX
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;
              z = DBLE(isign(1,ih)) * zb_fam(abs(ih));
              phiz(ih,nz) = phi_nz( nz , b0_mes*bz_mes , z );
          enddo
      enddo
      do nr = 0 , NMESMAX/2
          do il = 1 , NGL
              r = rb_fam(il);
              K = K_multipole;
              phirK(il,nr) = phi_nr_ml( nr , K , b0_mes*bp_mes , r );
          enddo
      enddo
      do i = 1 , NHSIZE
          nz = nz_mes(i);
          nr = nr_mes(i);

          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  ihl = ihl + 1;
                  Qsig( i , ihl ) = phiz(ih,nz) * phirK(il,nr);
                  Qome( i , ihl ) = phiz(ih,nz) * phirK(il,nr);
                  Qrho( i , ihl ) = phiz(ih,nz) * phirK(il,nr);

              enddo
          enddo

      enddo
      do i = 1 , NHSIZE
          do j = 1 , NHSIZE
              Hsig(i,j) = H(i,j);
              Home(i,j) = H(i,j);
              Hrho(i,j) = H(i,j);
          enddo
          Hsig(i,i) = Hsig(i,i) + (m_sig/hbc)**2.D0;
          Home(i,i) = Home(i,i) + (m_ome/hbc)**2.D0;
          Hrho(i,i) = Hrho(i,i) + (m_rho/hbc)**2.D0;
      enddo

      call dposv( 'U'   , NHSIZE , NCOORD ,
     &             Hsig , NHSIZE ,
     &             Qsig , NHSIZE , INFO1   );

      call dposv( 'U'   , NHSIZE , NCOORD ,
     &             Home , NHSIZE ,
     &             Qome , NHSIZE , INFO2   );

      call dposv( 'U'   , NHSIZE , NCOORD ,
     &             Hrho , NHSIZE ,
     &             Qrho , NHSIZE , INFO3   );

      call assert( INFO1.eq.0 , 'INFO =/= 0 in dposv()' );
      call assert( INFO2.eq.0 , 'INFO =/= 0 in dposv()' );
      call assert( INFO3.eq.0 , 'INFO =/= 0 in dposv()' );

      do i = 1 , NHSIZE
          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  Qsig(i,ihl) = Qsig(i,ihl) * wzwr(abs(ih),il);
                  Qome(i,ihl) = Qome(i,ihl) * wzwr(abs(ih),il);
                  Qrho(i,ihl) = Qrho(i,ihl) * wzwr(abs(ih),il);
              enddo
          enddo
      enddo






c-----Calculation of {QKsig,QKome,QKrho} matrices
      do nz = 0 , NMESMAX
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;
              z = DBLE(isign(1,ih)) * zb_famK(abs(ih));
              phiz(ih,nz) = phi_nz( nz , b0_mes*bz_mes , z );
          enddo
      enddo
      do nr = 0 , NMESMAX/2
          do il = 1 , NGL
              r = rb_famK(il);
              K = K_multipole;
              phirK(il,nr) = phi_nr_ml( nr , K , b0_mes*bp_mes , r );
          enddo
      enddo
      do i = 1 , NHSIZE
          nz = nz_mes(i);
          nr = nr_mes(i);

          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  ihl = ihl + 1;
                  QKsig( i , ihl ) = phiz(ih,nz) * phirK(il,nr);
                  QKome( i , ihl ) = phiz(ih,nz) * phirK(il,nr);
                  QKrho( i , ihl ) = phiz(ih,nz) * phirK(il,nr);

              enddo
          enddo

      enddo
      do i = 1 , NHSIZE
          do j = 1 , NHSIZE
              Hsig(i,j) = H(i,j);
              Home(i,j) = H(i,j);
              Hrho(i,j) = H(i,j);
          enddo
          Hsig(i,i) = Hsig(i,i) + (m_sig/hbc)**2.D0;
          Home(i,i) = Home(i,i) + (m_ome/hbc)**2.D0;
          Hrho(i,i) = Hrho(i,i) + (m_rho/hbc)**2.D0;
      enddo

      call dposv( 'U'    , NHSIZE , NCOORD ,
     &             Hsig  , NHSIZE ,
     &             QKsig , NHSIZE , INFO1   );

      call dposv( 'U'    , NHSIZE , NCOORD ,
     &             Home  , NHSIZE ,
     &             QKome , NHSIZE , INFO2   );

      call dposv( 'U'    , NHSIZE , NCOORD ,
     &             Hrho  , NHSIZE ,
     &             QKrho , NHSIZE , INFO3   );

      call assert( INFO1.eq.0 , 'INFO =/= 0 in dposv()' );
      call assert( INFO2.eq.0 , 'INFO =/= 0 in dposv()' );
      call assert( INFO3.eq.0 , 'INFO =/= 0 in dposv()' );

      do i = 1 , NHSIZE
          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  QKsig(i,ihl) = QKsig(i,ihl) * wzwrK(abs(ih),il);
                  QKome(i,ihl) = QKome(i,ihl) * wzwrK(abs(ih),il);
                  QKrho(i,ihl) = QKrho(i,ihl) * wzwrK(abs(ih),il);
              enddo
          enddo
      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END init_mesons() ****************************';
      write(6,*) '';
      endif

      return;
      end;
