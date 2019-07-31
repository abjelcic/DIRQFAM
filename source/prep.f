
      program PREP

      INTEGER*4 tape;
      INTEGER*4 n0f, n0b;
      INTEGER*4 NGH, NGL;
      INTEGER*4 J_multipole, K_multipole;
      INTEGER*4 broyd_n, broyd_m;
      INTEGER*4 NWMAX;
      INTEGER*4 NBSX;
      INTEGER*4 KTRUNC, NCOORD;
      INTEGER*4 NTX, KX, MVX, MVTX;

      parameter( NSIZE = 100000 );
      INTEGER*4 nz_spx(NSIZE), nr_spx(NSIZE), ml_spx(NSIZE);
      CHARACTER fg_spx(NSIZE);
      CHARACTER fg1, fg2;
      INTEGER*4 ml, ml1, ml2;
      INTEGER*4 nz, NZZ, nz1, nz2;
      INTEGER*4 nr, NRR, nr1, nr2;
      INTEGER*4 nf, ng;
      INTEGER*4 i, j;
      INTEGER*4 suma, lam, dimf, dimg, N;






c-----Reading n0f, n0b, NGH, NGL, J_multipole, K_multipole
      tape = 100;
      open( tape , file = 'dirqfam.dat' , status = 'old' );
      read(tape,'(10x,2i5)') n0f, n0b;
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,*);
      read(tape,'(18x,i9)') NGH;
      read(tape,'(18x,i9)') NGL;
      read(tape,*);
      read(tape,'(18x,i9)') J_multipole;
      read(tape,'(18x,i9)') K_multipole;
      close(tape);

      if( MOD(n0f,2) .ne. 0 ) then
          stop 'Error: n0f has to be even number!';
      endif






c-----Construction of quantum numbers
      nf = ((n0f+1)*(n0f+2)*(n0f+3))/6;
      ng = ((n0f+2)*(n0f+3)*(n0f+4))/6;
      if( nf+ng .gt. NSIZE ) then
          stop 'Error: nf+ng>NSIZE!';
      endif
      i = 0;
      do nz = 0 , n0f
          do nr = 0 , n0f/2
              do ml = -n0f , n0f
                  if( nz+2*nr+abs(ml) .le. n0f ) then
                      i = i + 1;
                      fg_spx(i) = 'f';
                      nz_spx(i) = nz;
                      nr_spx(i) = nr;
                      ml_spx(i) = ml;
                  endif
              enddo
          enddo
      enddo
      do nz = 0 , n0f+1
          do nr = 0 , (n0f+1)/2
              do ml = -(n0f+1) , n0f+1
                  if( nz+2*nr+abs(ml) .le. n0f+1 ) then
                      i = i + 1;
                      fg_spx(i) = 'g';
                      nz_spx(i) = nz;
                      nr_spx(i) = nr;
                      ml_spx(i) = ml;
                  endif
              enddo
          enddo
      enddo
      if( nf+ng .ne. i ) then
          stop 'Error: nf+ng!';
      endif






c-----Calculation of maximal block dimension (NBSX) of the U,V matrices
c-----In fact, this is the dimension of the omega^pi = 1/2^+ block.
      NBSX = 0;
      ! f-component
      do nz = 0 , n0f
          do nr = 0 , n0f/2
              if( MOD(nz,2) .eq. 0 ) then
                  if( 2*nr+nz .le. n0f   ) NBSX = NBSX + 1;
              else
                  if( 2*nr+nz .le. n0f-1 ) NBSX = NBSX + 1;
              endif
          enddo
      enddo
      ! g-component
      do nz = 0 , n0f+1
          do nr = 0 , (n0f+1)/2
              if( MOD(nz,2) .eq. 0 ) then
                  if( 2*nr+nz .le. n0f   ) NBSX = NBSX + 1;
              else
                  if( 2*nr+nz .le. n0f+1 ) NBSX = NBSX + 1;
              endif
          enddo
      enddo






c-----Calculating length of the Broyden vector
      broyd_m = 20;
      broyd_n = 0;
      ! dh_1 matrix
      do i = 1 , nf+ng
          fg1 = fg_spx(i);
          ml1 = ml_spx(i);
          do j = i , nf+ng
              fg2 = fg_spx(j);
              ml2 = ml_spx(j);
              if( fg1 .eq. fg2 ) then
                  if( abs(ml1-ml2) .eq. K_multipole ) then
                      broyd_n = broyd_n + 1;
                  endif
              endif
              if( fg1 .ne. fg2 ) then
                  if( abs(ml1-ml2) .eq. K_multipole ) then
                      broyd_n = broyd_n + 1;
                  endif
                  if( abs(ml1+ml2+1) .eq. K_multipole ) then
                      broyd_n = broyd_n + 1;
                  endif
              endif
          enddo
      enddo
      ! dDelta_pl, dDelta_mi
      do i = 1 , nf+ng
          fg1 = fg_spx(i);
          ml1 = ml_spx(i);
          do j = i , nf+ng
              fg2 = fg_spx(j);
              ml2 = ml_spx(j);
              if( fg1.ne.'f' .or. fg2.ne.'f'    ) CYCLE;
              if( abs(ml1-ml2) .ne. K_multipole ) CYCLE;
              broyd_n = broyd_n + 1; !dDelta_pl
              broyd_n = broyd_n + 1; !dDelta_mi
          enddo
      enddo
      broyd_n = 2*broyd_n; !protons/neutrons
      broyd_n = 2*broyd_n; !real/imaginary






c-----Calculating NWMAX
      NWMAX = 0;
      do i = 1 , nf+ng
          fg1 = fg_spx(i);
          nz1 = nz_spx(i);
          nr1 = nr_spx(i);
          ml1 = ml_spx(i);

          do j = i , nf+ng
              fg2 = fg_spx(j);
              nz2 = nz_spx(j);
              nr2 = nr_spx(j);
              ml2 = ml_spx(j);

              if( fg1.ne.'f' .or. fg2.ne.'f'    ) CYCLE;
              if( abs(ml1-ml2) .ne. K_multipole ) CYCLE;

              do NRR = 0 , nr1+nr2+(abs(ml1)+abs(ml2)-K_multipole)/2
                  do NZZ = 0 , nz1+nz2
                      if( mod(NZZ,2) .ne. mod(nz1+nz2,2) ) CYCLE;

                      NWMAX = NWMAX + 1;

                  enddo
              enddo

          enddo
      enddo






c-----Setting KTRUNC and NCOORD
      NCOORD = (2*NGH)*NGL;
      KTRUNC = 260;






c-----Calculating NTX, KX, MVX, MVTX
      NTX = (n0f+2)*(n0f+3)*(2*n0f+5)/6;
      KX  = (n0f+2)*(n0f+3)*(n0f+4)/6;
      ! f-component
      dimf = 0;
      ! Positive parity
      do i = 1 , n0f+1
          suma = 0;
          do lam = i-1 , i
              do N = lam+mod(lam,2) , n0f , 2
                 suma = suma + (N-lam)/2 + 1;
              enddo
          enddo
          dimf = dimf + (suma*(suma+1))/2;
      enddo
      ! Negative parity
      do i = 1 , n0f+1
          suma = 0;
          do lam = i-1 , i
              do N = lam+mod(lam+1,2) , n0f , 2
                  suma = suma + (N-lam)/2 + 1;
              enddo
          enddo
          dimf = dimf + (suma*(suma+1))/2;
      enddo
      MVX = dimf;
      ! g-component
      dimg = 0;
      ! Positive parity
      do i = 1 , (n0f+1) + 1
          suma = 0;
          do lam = i-1 , i
              do N = lam+mod(lam,2) , (n0f+1) , 2
                  suma = suma + (N-lam)/2 + 1;
              enddo
          enddo
          dimg = dimg + (suma*(suma+1))/2;
      enddo
      ! Negative parity
      do i = 1 , (n0f+1)+1
          suma = 0;
          do lam = i-1 , i
              do N = lam+mod(lam+1,2) , (n0f+1) , 2
                  suma = suma + (N-lam)/2 + 1;
              enddo
          enddo
          dimg = dimg + (suma*(suma+1))/2;
      enddo
      MVTX = dimf + dimg;






c-----Generating dirqfam.par file
      tape = 100;
      open( tape , file = './source/dirqfam.par' , status = 'unknown' );

      write(tape,'(a)') 'c-----maximal number for GFV';
      write(tape,'(6x,a)') 'parameter ( IGFV =       100 )';
      write(tape,'(a)') '';
      write(tape,'(a)') 'c-----number of r-meshpoints (4n+1 points)';
      write(tape,'(6x,a)') 'parameter ( RMAX =  15.00000 )';
      write(tape,'(6x,a)') 'parameter ( MR =       301 )';
      write(tape,'(a)') '';
      write(tape,'(a)') 'c-----number of q-meshpoints (4n+1 points)';
      write(tape,'(6x,a)') 'parameter ( QMAX =   6.00000 )';
      write(tape,'(6x,a)') 'parameter ( MQ =        49 )';
      write(tape,'(a)') '';
      write(tape,'(a)') 'c-----number of gauss-meshpoints';
      write(tape,'(6x,a,i2,a)') 'parameter ( NGH =        ',NGH,' )';
      write(tape,'(6x,a,i2,a)') 'parameter ( NGL =        ',NGL,' )';
      write(tape,'(6x,a)') 'parameter ( NGLEG =      18 )';
      write(tape,'(a)') ' ';
      write(tape,'(2a)')'c---- maximal oscillator quantum number ',
     &            'for fermions';
      write(tape,'(6x,a,i2,a)') 'parameter ( N0FX =       ',n0f,' )';
      write(tape,'(6x,a)') 'parameter ( nxx  =       N0FX/2 )';
      write(tape,'(a)') ' ';
      write(tape,'(a)') 'c-----maximal number of (k,parity)-blocks';
      write(tape,'(6x,a)') 'parameter ( NBX =  2*N0FX+3 )';
      write(tape,'(a)') '';
      write(tape,'(2a)') 'c-----max. number of all levels for ',
     &                   'protons or neutrons';
      write(tape,'(6x,a,i8,a)') 'parameter ( NTX = ' , NTX , ' )';
      write(tape,'(1a)') '';
      write(tape,'(2a)') 'c-----max. number of eigenstates',
     &                   ' for protons or neutrons';
      write(tape,'(6x,a,i6,a)') 'parameter ( KX  = ' , KX , ' )';
      write(tape,'(a)') '';
      write(tape,'(a)') 'c-----max. nz-quantum number of fermions';
      write(tape,'(6x,a)') 'parameter ( NZX  =   N0FX+1 )';
      write(tape,'(a)') '';
      write(tape,'(a)') 'c-----max. nr-quantum number of fermions';
      write(tape,'(6x,a)') 'parameter ( NRX  = N0FX/2 )';
      write(tape,'(a)') '';
      write(tape,'(a)') 'c-----max. ml-quantum number of fermions';
      write(tape,'(6x,a)') 'parameter ( MLX  =        N0FX+1 )';
      write(tape,'(a)') '';
      write(tape,'(a)') 'c-----maximal dimension F of one k-block';
      write(tape,'(6x,a)') 'parameter ( NFX  = (N0FX/2+1)**2 )';
      write(tape,'(a)') '';
      write(tape,'(a)') 'c-----maximal dimension G of one k-block';
      write(tape,'(6x,a)') 'parameter ( NGX  = (N0FX/2+1)*(N0FX/2+2) )';
      write(tape,'(a)') '';
      write(tape,'(a)') 'c-----maximum of nf and ng in all blocks';
      write(tape,'(6x,a)') 'parameter ( NDX  =       NGX )';
      write(tape,'(a)') '';
      write(tape,'(a)') 'c-----oscillator quantum number for bosons';
      write(tape,'(6x,a)') 'parameter ( N0BX =        20 )';
      write(tape,'(a)') '';
      write(tape,'(a)') 'c-----number of bosons';
      write(tape,'(6x,a)') 'parameter ( nbxx = N0BX/2 )';
      write(tape,'(6x,a)') 'parameter ( NOX  = (nbxx+1)*(nbxx+2)/2 )';
      write(tape,'(a)') '';
      write(tape,'(a)') 'c-----for the plot';
      write(tape,'(6x,a)') 'parameter ( NOX1 = (N0FX+1)*(N0FX+2)/2 )';
      write(tape,'(a)') '';
      write(tape,'(a)') 'c-----auxiliary';
      write(tape,'(6x,a)') 'parameter ( NGH2   = NGH+NGH )';
      write(tape,'(6x,a)') 'parameter ( NB2X   = NBX+NBX )';
      write(tape,'(6x,a)') 'parameter ( NHX    = NFX+NGX )';
      write(tape,'(6x,a)') 'parameter ( NDDX   = NDX*NDX )';
      write(tape,'(6x,a)') 'parameter ( NHBX   = NHX+NHX )';
      write(tape,'(6x,a)') 'parameter ( NHBQX  = NHBX*NHBX )';
      write(tape,'(6x,a)') 'parameter ( NFFX   = NFX*NFX )';
      write(tape,'(6x,a)') 'parameter ( NFGX   = NFX*NGX )';
      write(tape,'(6x,a)') 'parameter ( NHHX   = NHX*NHX )';
      write(tape,'(6x,a)') 'parameter ( MG     = (NGH+1)*(NGL+1) )';
      write(tape,'(6x,a)') 'parameter ( N02    = 2*N0FX )';
      write(tape,'(6x,a)') 'parameter ( NNNX   = (N0FX+1)*(N0FX+1) )';
      write(tape,'(a)') '';
      write(tape,'(a)') 'c-----working space';
      write(tape,'(6x,a,i8,a)') 'parameter ( MVX  = ' , MVX  , ' )';
      write(tape,'(6x,a,i8,a)') 'parameter ( MVTX = ' , MVTX , ' )';
      write(tape,'(a)') '';
      write(tape,'(a)') 'c-----FAM parameters';
      write(tape,'(6x,a,i9,a)') 'parameter( NBSX       =', NBSX, ' )';
      write(tape,'(6x,a)') 'parameter( J_MAX      =        3 )';
      write(tape,'(6x,a,i9,a)') 'parameter( NFAM_BROYD =',broyd_n,' )';
      write(tape,'(6x,a,i9,a)') 'parameter( MFAM_BROYD =',broyd_m,' )';
      write(tape,'(6x,a,i9,a)') 'parameter( NWMAX      =',NWMAX,' )';
      write(tape,'(6x,a,i9,a)') 'parameter( KTRUNC     =',KTRUNC,' )';
      write(tape,'(6x,a,i9,a)') 'parameter( NCOORD     =',NCOORD,' )';
      write(tape,'(6x,a,i9,a)') 'parameter( JCHECK     =',
     &                           J_multipole,' )';
      write(tape,'(6x,a,i9,a)') 'parameter( KCHECK     =',
     &                           K_multipole,' )';
      close(tape);

      write(6,*) '';
      write(6,'(a)') 'Successfully generated dirqfam.par file.';
      write(6,*) '';

      return;
      end;
