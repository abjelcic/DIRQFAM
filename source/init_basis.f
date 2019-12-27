c======================================================================c

      subroutine init_basis( lpr )

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

      CHARACTER fg_spx;
      common /simplex/ N_total         , N_blocks        ,
     &                 ia_spx(NBX)     , id_spx(NBX)     ,
     &                 nf_size(NBX)    , ng_size(NBX)    ,
     &                 nz_spx(NBSX,NBX), nr_spx(NBSX,NBX),
     &                 ml_spx(NBSX,NBX), fg_spx(NBSX,NBX);

      common /quadrature/ zb_fam( 1:NGH ), wz( 1:NGH ),
     &                    rb_fam( 1:NGL ), wr( 1:NGL ),
     &                    wzwr( 1:NGH , 1:NGL );

      common /basis/ phi_z( -NGH:NGH , NTX ),
     &              dphi_z( -NGH:NGH , NTX ),
     &               phi_r(    1:NGL , NTX ),
     &              dphi_r(    1:NGL , NTX );

      common /wbasis/ wPhi( 1:NGH , 1:NGL , NTX );

      common /PHI/ PHI_U  (    NTX , KTRUNC ),
     &             PHI_SVt( KTRUNC , NCOORD ),
     &             k_PHI;



      REAL*8 xgh(2*NGH), wgh(2*NGH);
      REAL*8 xgl(1*NGL), wgl(1*NGL);

      REAL*8 PHI  ( NTX , NCOORD );
      REAL*8 dzPHI( NTX , NCOORD );
      REAL*8 drPHI( NTX , NCOORD );



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN init_basis() **************************';
      write(6,*) '';
      endif






c-----Calculation of quadrature mesh and weights
      alpha = 0.D0;
      call gauher( xgh , wgh , 2*NGH );
      call gaulag( xgl , wgl , 1*NGL , alpha );
      do ih = 1 , NGH
          zb_fam(ih) = b0*bz * xgh(1+NGH-ih);
              wz(ih) = b0*bz * wgh(1+NGH-ih)
     &                       * DEXP(xgh(1+NGH-ih)**2.D0);
      enddo
      do il = 1 , NGL
          rb_fam(il) = b0*bp * DSQRT(xgl(il));
              wr(il) = 0.5D0 * b0*b0*bp*bp * wgl(il)
     &                       * DEXP(xgl(il)) / xgl(il)**alpha;
      enddo
      do il = 1 , NGL
          do ih = 1 , NGH
              wzwr(ih,il) = wz(ih)*wr(il);
          enddo
      enddo






c-----Calculation of phi_z
      do ib = 1 , N_blocks
          do i = 1 , id_spx(ib)
              nz = nz_spx(i,ib);
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  ii = i-1+ia_spx(ib);

                  z = DBLE(isign(1,ih)) * zb_fam(abs(ih));
                  phi_z(ih,ii) = phi_nz(nz,b0*bz,z);

              enddo
          enddo
      enddo






c-----Calculation of phi_r
      do ib = 1 , N_blocks
          do i = 1 , id_spx(ib)
              nr = nr_spx(i,ib);
              ml = ml_spx(i,ib);
              do il = 1 , NGL

                  ii = i-1+ia_spx(ib);

                  r = rb_fam(il);
                  phi_r(il,ii) = phi_nr_ml(nr,abs(ml),b0*bp,r);

              enddo
          enddo
      enddo






c-----Calculation of dphi_z
      do ib = 1 , N_blocks
          do i = 1 , id_spx(ib)
              nz = nz_spx(i,ib);
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  ii = i-1+ia_spx(ib);

                  z = DBLE(isign(1,ih)) * zb_fam(abs(ih));
                  dphi_z(ih,ii) = d_phi_nz(nz,b0*bz,z);

              enddo
          enddo
      enddo






c-----Calculation of dphi_r
      do ib = 1 , N_blocks
          do i = 1 , id_spx(ib)
              nr = nr_spx(i,ib);
              ml = ml_spx(i,ib);
              do il = 1 , NGL

                  ii = i-1+ia_spx(ib);

                  r = rb_fam(il);
                  dphi_r(il,ii) = d_phi_nr_ml(nr,abs(ml),b0*bp,r);

              enddo
          enddo
      enddo






c-----Calculation of wPhi
      do i = 1 , N_total
          do il = 1 , NGL
              do ih = 1 , NGH

                  call assert(wzwr(ih,il).ge.0.D0,'negative weights');

                  wPhi(ih,il,i) =   DSQRT(wzwr(ih,il))
     &                            * phi_z(ih,i)
     &                            * phi_r(il,i);

              enddo
          enddo
      enddo






c-----Calculation of PHI_U, PHI_SVt
      do i = 1 , N_total
          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  ihl = ihl+1;
                  PHI(i,ihl) = phi_z(ih,i)*phi_r(il,i);

              enddo
          enddo
      enddo
      call lowrank_approx(  N_total , NCOORD ,
     &                      1.D-12  ,
     &                      KTRUNC  , k_PHI  ,
     &                      PHI     , NTX    ,
     &                      PHI_U   , NTX    ,
     &                      PHI_SVt , KTRUNC   );






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END init_basis() ****************************';
      write(6,*) '';
      endif

      return;
      end;






c======================================================================c

      subroutine lowrank_approx(     M ,      N  ,
     &                             TOL ,
     &                               K ,    RANK ,
     &                               A ,     LDA ,
     &                             A_U ,    LDAU ,
     &                           A_SVt ,  LDASVt   )

c======================================================================c
      ! Calculates the low-rank approximation of MxN matrix A using
      ! truncated SVD approach with tolerance TOL: A ~= A_U * A_SVt.
      !
      !
      !  [in]:  A     is real M x N matrix (LDA >= M)
      !
      !         TOL   is a threshold below which singular values of
      !               matrix A are neglected (TOL >= 0)
      !
      !         K     is the upper bound on rank of approximation matrix
      !               A_U*A_SVt, routine fails if not large enough
      !
      !
      !
      ! [out]:  A_U   is real M x K matrix (LDAU >= M) with only
      !               M x RANK relevant submatrix
      !
      !         A_SVt is real K x N matrix (LDASVt >= K) with only
      !               RANK x N relevant submatrix
      !
      !         RANK  is the rank of approximation matrix A_U*A_SVt
      !

      IMPLICIT NONE;
      LOGICAL CHECK;

      INTEGER*4 M, N, K, RANK, LDA, LDAU, LDASVt, NWORK, i, j, INFO;
      REAL*8 TOL, normf1, normf2, error;
      REAL*8 A    ( LDA    , N );
      REAL*8 A_cpy( M      , N );
      REAL*8 A_U  ( LDAU   , K );
      REAL*8 A_SVt( LDASVt , N );

      REAL*8 S ( min(M,N) );
      REAL*8 U ( M , min(M,N) );
      REAL*8 VT( min(M,N) , N );

      REAL*8     WORK( 10000*M );
      INTEGER*4 IWORK( 12*min(M,N) );
      NWORK = 10000*M;
      CHECK = .false.;



      if( CHECK .eqv. .true. ) then
          do j = 1 , N
              do i = 1 , M
                  A_cpy(i,j) = A(i,j);
              enddo
          enddo
      endif



      call assert( M.gt.0       , 'M <= 0 in lowrank_approx()'     );
      call assert( N.gt.0       , 'N <= 0 in lowrank_approx()'     );
      call assert( TOL.gt.0.D0  , 'TOL <= 0 in lowrank_approx()'   );
      call assert( LDA.ge.M     , 'LDA < M in lowrank_approx()'    );
      call assert( LDAU.ge.M    , 'LDAU < M in lowrank_approx()'   );
      call assert( LDASVt.ge.K  , 'LDASVt < K in lowrank_approx()' );



      call dgesvdx(   'V' ,       'V' ,
     &                'V' ,
     &                  M ,         N ,
     &                  A ,       LDA ,
     &                TOL ,   1.D+100 ,
     &                  0 ,         0 ,
     &               RANK ,         S ,
     &                  U ,         M ,
     &                  VT,  min(M,N) ,
     &               WORK ,        -1 ,
     &              IWORK ,      INFO   );

      call assert( INT(WORK(1)).le.NWORK , 'NWORK too small' );

      call dgesvdx(   'V' ,       'V' ,
     &                'V' ,
     &                  M ,         N ,
     &                  A ,       LDA ,
     &                TOL ,   1.D+100 ,
     &                  0 ,         0 ,
     &               RANK ,         S ,
     &                  U ,         M ,
     &                  VT,  min(M,N) ,
     &               WORK ,     NWORK ,
     &              IWORK ,      INFO   );

      call assert( INFO.eq.0 , 'INFO =/= 0 in dgesvdx()' );
      call assert( RANK.le.K , 'K too small' );


      do j = 1 , RANK
          do i = 1 , M
              A_U(i,j) = U(i,j);
          enddo
      enddo

      do j = 1 , N
          do i = 1 , RANK
              A_SVt(i,j) = S(i) * VT(i,j);
          enddo
      enddo



      if( CHECK .eqv. .true. ) then

          normf1 = 0.D0;
          do j = 1 , N
              do i = 1 , M
                  normf1 = normf1 + A_cpy(i,j)**2.D0;
              enddo
          enddo
          normf1 = DSQRT( normf1 );

          call dgemm(  'N'    , 'N'      ,
     &                  M     ,  N       ,  RANK  ,
     &                 +1.D0  ,
     &                 A_U    ,  LDAU ,
     &                 A_SVt  ,  LDASVt  ,
     &                 -1.D0  ,
     &                 A_cpy  ,  M                 );

          normf2 = 0.D0;
          do j = 1 , N
              do i = 1 , M
                  normf2 = normf2 + A_cpy(i,j)**2.D0;
              enddo
          enddo
          normf2 = DSQRT( normf2 );

          error = normf2/normf1;

          write(6,*) 'Low-rank approximation error = ' , error;

      endif

      return;
      end;
