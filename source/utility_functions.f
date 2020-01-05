c======================================================================c

      subroutine assert( statement , error_message )

c======================================================================c

      LOGICAL statement;
      CHARACTER( LEN = * ) error_message;

      if( statement .eqv. .false. ) then
          write(6,'(3a)') 'Error: ', error_message, '!';
          stop;
      endif

      return;
      end;






c======================================================================c

      INTEGER*4 function index_of_vector( fg , nz , nr , ml , ib )

c======================================================================c

      IMPLICIT REAL*8    (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n)
      include 'dirqfam.par'

      CHARACTER fg_spx;
      common /simplex/ N_total         , N_blocks        ,
     &                 ia_spx(NBX)     , id_spx(NBX)     ,
     &                 nf_size(NBX)    , ng_size(NBX)    ,
     &                 nz_spx(NBSX,NBX), nr_spx(NBSX,NBX),
     &                 ml_spx(NBSX,NBX), fg_spx(NBSX,NBX);

      CHARACTER fg;
      INTEGER*4 nz, nr, ml, ib, i;
      LOGICAL bool;

      bool = .false.;
      do i = 1 , id_spx(ib)

          if( ( ml .eq. ml_spx(i,ib) ) .and.
     &        ( nr .eq. nr_spx(i,ib) ) .and.
     &        ( nz .eq. nz_spx(i,ib) ) .and.
     &        ( fg .eq. fg_spx(i,ib) ) ) then
              bool = .true.;
              EXIT;
          endif

      enddo

      call assert( bool .eqv. .true. , 'index_of_vector() wrong' );

      index_of_vector = i;

      return;
      end;






c======================================================================c

      REAL*8 function phi_nz( nz , b , z )

c======================================================================c
c----------------------------------------------------------------------c
c                                                                      c
c  phi_nz(nz,b,z) =   1/sqrt(b) * 1/sqrt(sqrt(pi)*2^nz*nz!)            c
c                   * H_nz(z/b) * exp( -1/2 * (z/b)^2 )                c
c----------------------------------------------------------------------c

      IMPLICIT REAL*8    (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n)

      INTEGER*4 nz;
      REAL*8 b, z;
      pi = 3.14159265358979324D0;

      call assert( nz .ge. 0    , 'nz <  0 in phi_nz()' );
      call assert( b  .gt. 0.D0 , 'b  <= 0 in phi_nz()' );

      x = z/b;
      phi_0 = 1.D0 / DSQRT(b*DSQRT(pi)) * DEXP( -0.5D0 * x*x );
      phi_1 = x * DSQRT(2.D0) * phi_0;

      if( nz .eq. 0 ) then
          phi_nz = phi_0;
          return;
      endif


      do n = 2 , nz
          tmp   = phi_1;
          phi_1 = phi_1 * DSQRT( 2.D0 / DBLE(n) ) * x;
          phi_1 = phi_1 - DSQRT( DBLE(n-1) / DBLE(n) ) * phi_0;
          phi_0 = tmp;
      enddo

      phi_nz = phi_1;

      return;
      end;






c======================================================================c

      REAL*8 function d_phi_nz( nz , b , z )

c======================================================================c
c----------------------------------------------------------------------c
c                                                                      c
c  d_phi_nz(nz,b,z) = d( phi_nz(nz,b,z) )/dz                           c
c                                                                      c
c----------------------------------------------------------------------c

      IMPLICIT REAL*8    (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n)

      INTEGER*4 nz;
      REAL*8 b, z;

      call assert( nz .ge. 0    , 'nz <  0 in d_phi_nz()' );
      call assert( b  .gt. 0.D0 , 'b  <= 0 in d_phi_nz()' );

      x = z/b;
      if( nz .eq. 0 ) then
          d_phi_nz = -x/b * phi_nz(nz,b,z);
          return;
      endif

      d_phi_nz = - (x/b) * phi_nz(nz,b,z)
     &           + DSQRT(DBLE(2*nz))/b * phi_nz(nz-1,b,z);

      return;
      end;






c======================================================================c

      REAL*8 function phi_nr_ml( nr , ml , b , r )

c======================================================================c
c----------------------------------------------------------------------c
c                                                                      c
c  phi_nr_ml(nr,ml,b,r) =   1/b * sqrt( 2*nr! / (nr+|ml|)! )           c
c                         * (r/b)^|ml| * L_{nr}^|ml|( (r/b)^2 )        c
c                         *  exp( -1/2 * (r/b)^2 )                     c
c----------------------------------------------------------------------c

      IMPLICIT REAL*8    (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n)

      parameter( IGFV = 100 )
      common /gfvfak/ fak(0:IGFV);

      INTEGER*4 nr, ml;
      REAL*8 b, r;

      call assert( nr .ge. 0    , 'nr <  0 in phi_nr_ml()' );
      call assert( b  .gt. 0.D0 , 'b  <= 0 in phi_nr_ml()' );
      call assert( r  .gt. 0.D0 , 'r  <= 0 in phi_nr_ml()' );
      call assert( abs(ml).le.IGFV , '|ml| too large in phi_nr_ml()' );

      ml = abs(ml);
      eta = r*r/(b*b);

      phi_0_ml = DSQRT(2.D0) / ( b* DSQRT(fak(ml)) ) *
     &           DEXP( -0.5D0* ( eta - DBLE(ml)*DLOG(eta) ) );
      phi_1_ml =  phi_0_ml * ( -eta + DBLE(ml+1) ) / DSQRT(DBLE(ml+1));

      if( nr .eq. 0 ) then
          phi_nr_ml = phi_0_ml;
          return;
      endif

      do n = 2 , nr
          tmp      = phi_1_ml;
          phi_1_ml = phi_1_ml*(DBLE(2*n-1+ml)-eta)/DSQRT(DBLE(n*(n+ml)))
          phi_0_ml = phi_0_ml*DSQRT(DBLE((n-1)*(n-1+ml))/DBLE(n*(n+ml)))
          phi_1_ml = phi_1_ml - phi_0_ml;
          phi_0_ml = tmp;
      enddo

      phi_nr_ml = phi_1_ml;

      return;
      end;






c======================================================================c

      REAL*8 function d_phi_nr_ml( nr , ml , b , r )

c======================================================================c
c----------------------------------------------------------------------c
c                                                                      c
c  d_phi_nr_ml(nr,ml,b,z) = d( phi_nr_ml(nr,ml,b,r) )/dr               c
c                                                                      c
c----------------------------------------------------------------------c

      IMPLICIT REAL*8    (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n)

      INTEGER*4 nr, ml;
      REAL*8 b, r;

      call assert( nr .ge. 0    , 'nr <  0 in d_phi_nr_ml()' );
      call assert( b  .gt. 0.D0 , 'b  <= 0 in d_phi_nr_ml()' );
      call assert( r  .gt. 0.D0 , 'r  <= 0 in d_phi_nr_ml()' );

      ml = abs(ml);
      eta = r*r/(b*b);

      if( nr .eq. 0 ) then
          d_phi_nr_ml = ( DBLE(ml) - eta )/( DSQRT(eta) * b ) *
     &                  phi_nr_ml(0,ml,b,r);
          return;
      endif

      fac1 = ( DBLE(2*nr+ml) - eta ) / DSQRT(eta);
      fac2 = - 2.D0 * DSQRT( DBLE( nr*(nr+ml) ) / eta );

      d_phi_nr_ml = + fac1/b * phi_nr_ml(nr  ,ml,b,r)
     &              + fac2/b * phi_nr_ml(nr-1,ml,b,r);

      return;
      end;






c======================================================================c

      REAL*8 function I_K( a , K )

c======================================================================c
c----------------------------------------------------------------------c
c Calculates the integral of sqrt(1-a*cos(x)^2)*cos(2*K*x) on [0,pi/2] c
c Parameter a is in [0,1]                                              c
c Parameter K is in {0,1,2,3,4}                                        c
c----------------------------------------------------------------------c
      IMPLICIT NONE;
      REAL*8 a, x, Ec, Kc, pi, fac1, fac2, fac3;
      INTEGER*4 K, n;

      INTEGER*4 LMT;
      parameter( LMT = 30 );

      REAL*8  K0(0:LMT), K1(0:LMT), K2(0:LMT), K3(0:LMT), K4(0:LMT);
      REAL*8  epsi       /1.D-15/;
      LOGICAL first_call /.true./;

      SAVE first_call, epsi;
      SAVE K0, K1, K2, K3, K4;



      call assert( a.ge.0.D0 .and. a.le.1.D0 , 'a not in [0,1]'       );
      call assert( K.ge.0    .and. K.le.4    , 'K not in {0,1,2,3,4}' );



      if( first_call ) then
          first_call = .false.;

          K0 = 0.D0;
          K1 = 0.D0;
          K2 = 0.D0;
          K3 = 0.D0;
          K4 = 0.D0;

          pi = 3.14159265358979324D0;

          K0(0) = + pi *    1.D0 /     2.D0;
          K0(1) = - pi *    1.D0 /     8.D0;
          K0(2) = - pi *    3.D0 /   128.D0;
          K0(3) = - pi *    5.D0 /   512.D0;
          K0(4) = - pi *  175.D0 / 32768.D0;

          K1(1) = - pi *    1.D0 /    16.D0;
          K1(2) = - pi *    1.D0 /    64.D0;
          K1(3) = - pi *   15.D0 /  2048.D0;
          K1(4) = - pi *   35.D0 /  8192.D0;

          K2(2) = - pi *    1.D0 /   256.D0;
          K2(3) = - pi *    3.D0 /  1024.D0;
          K2(4) = - pi *   35.D0 / 16384.D0;

          K3(3) = - pi *    1.D0 /  2048.D0;
          K3(4) = - pi *    5.D0 /  8192.D0;

          K4(4) = - pi *    5.D0 / 65536.D0;

          do n = 5 , LMT
              K0(n) = K0(n-1) * DBLE((2*n-1)*(2*n-3))/DBLE(4*(n*n- 0));
              K1(n) = K1(n-1) * DBLE((2*n-1)*(2*n-3))/DBLE(4*(n*n- 1));
              K2(n) = K2(n-1) * DBLE((2*n-1)*(2*n-3))/DBLE(4*(n*n- 4));
              K3(n) = K3(n-1) * DBLE((2*n-1)*(2*n-3))/DBLE(4*(n*n- 9));
              K4(n) = K4(n-1) * DBLE((2*n-1)*(2*n-3))/DBLE(4*(n*n-16));
          enddo

      endif



      if( DABS( a - 1.D0 ) .lt. 1.D-14 ) then
          select case( K )
              case( 0 )
                  I_K = + 1.D0 /  1.D0;
              case( 1 )
                  I_K = - 1.D0 /  3.D0;
              case( 2 )
                  I_K = - 2.D0 / 30.D0;
              case( 3 )
                  I_K = - 1.D0 / 35.D0;
              case( 4 )
                  I_K = - 1.D0 / 63.D0;
              case default
                  stop 'Error: K > 4 in I_K()!';
          end select

          return;
      endif



      if( a .gt. 0.465D0 ) then

          call CElliptic( epsi , a , Kc , Ec , n );

          select case( K )
              case( 0 )
                  fac1 = 0.D0;
                  fac2 = 1.D0;
                  fac3 = 1.D0;
              case( 1 )
                  fac1 = 2*(1-a);
                  fac2 = a-2;
                  fac3 = 3*a;
              case( 2 )
                  fac1 = (8*(1-a)*(2-a));
                  fac2 = -(a*(a-16)+16);
                  fac3 = 15*a*a;
              case( 3 )
                  fac1 = 2*(1-a)*(-5*a*a+32*(a-2)*(a-2));
                  fac2 = (2-a)*(5*a*a-8*(a*(a-16)+16));
                  fac3 = 105*a*a*a;
              case( 4 )
                  fac1 = 32*(1-a)*(2-a)*(a*(5*a-32)+32);
                  fac2 = -( a*(a*(a*(5*a-448)+2496)-4096)+2048);
                  fac3 = 315*a*a*a*a;
              case default
                  stop 'Error: K > 4 in I_K()!';
          end select

          I_K = ( fac1*Kc + fac2*Ec ) / fac3;
          return;

      else

          I_K = 0.D0;
          x   = 1.D0;
          select case( K )
              case( 0 )
                  do n = 0 , LMT
                      I_K = I_K + x*K0(n);
                      x = x*a;
                  enddo
              case( 1 )
                  do n = 0 , LMT
                      I_K = I_K + x*K1(n);
                      x = x*a;
                  enddo
              case( 2 )
                  do n = 0 , LMT
                      I_K = I_K + x*K2(n);
                      x = x*a;
                  enddo
              case( 3 )
                  do n = 0 , LMT
                      I_K = I_K + x*K3(n);
                      x = x*a;
                  enddo
              case( 4 )
                  do n = 0 , LMT
                      I_K = I_K + x*K4(n);
                      x = x*a;
                  enddo
              case default
                  stop 'Error: K > 4 in I_K()!';
          end select

          return;

      endif



      return;
      end;






c======================================================================c

      REAL*8 function TalmiMoshinsky_1d( nz1 , nz2 ,
     &                                   NZZ , nz   )

c======================================================================c

      IMPLICIT REAL*8    (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n)
      include 'dirqfam.par'

      common /gfviv / iv ( -IGFV : IGFV ); ! iv(n)  = (-1)^n
      common /gfvfak/ fak(   0   : IGFV ); ! fak(n) = n!
      common /gfvfi / fi (   0   : IGFV ); ! fi(n)  = 1/n!
      common /gfvwf / wf (   0   : IGFV ); ! wf(n)  = sqrt(n!)
      common /gfvwfi/ wfi(   0   : IGFV ); ! wfi(n) = 1/sqrt(n!)

      REAL*8 fac, acc, bin1, bin2;
      INTEGER*4 nz1, nz2,
     &          NZZ, nz;
      INTEGER*4 p;



      TalmiMoshinsky_1d = 0.D0;



      if( nz1.lt.0 .or. nz2.lt.0 .or. NZZ.lt.0 .or. nz.lt.0 ) then
          return;
      endif
      if( nz1+nz2 .ne. NZZ+nz ) then
          return;
      endif

      call assert( nz1.le.IGFV , 'IGFV too small' );
      call assert( nz2.le.IGFV , 'IGFV too small' );
      call assert( NZZ.le.IGFV , 'IGFV too small' );
      call assert(  nz.le.IGFV , 'IGFV too small' );

      fac = wf(nz1)*wf(nz2)*wfi(NZZ)*wfi(nz);
      fac = fac / DSQRT(2.D0**(DBLE( nz1+nz2 )));

      acc = 0.D0;
      do p = max(0,nz1-NZZ) , min(nz,nz1)
          bin1 = fak(nz) *fi(p)    *fi(nz-p);
          bin2 = fak(NZZ)*fi(nz1-p)*fi(NZZ-nz1+p);

          acc = acc + iv(p)*bin1*bin2;
      enddo

      TalmiMoshinsky_1d = fac * acc;

      return;
      end;






c======================================================================c

      REAL*8 function TalmiMoshinsky_2d( nr1 , ml1 , nr2 , ml2 ,
     &                                   NRR , MLL , nr  , ml   )

c======================================================================c

      IMPLICIT REAL*8    (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n)
      include 'dirqfam.par'

      common /gfviv / iv ( -IGFV : IGFV ); ! iv(n)  = (-1)^n
      common /gfvfak/ fak(   0   : IGFV ); ! fak(n) = n!
      common /gfvfi / fi (   0   : IGFV ); ! fi(n)  = 1/n!
      common /gfvwf / wf (   0   : IGFV ); ! wf(n)  = sqrt(n!)
      common /gfvwfi/ wfi(   0   : IGFV ); ! wfi(n) = 1/sqrt(n!)

      REAL*8 fac, acc, mult1, mult2, bin1, bin2;
      INTEGER*4 nr1, ml1, nr2, ml2,
     &          NRR, MLL, nr , ml ;
      INTEGER*4 Mmin, Mmax, Nrmax;
      INTEGER*4 cond1, cond2;
      INTEGER*4 p , q , r , s ,
     &          PP, QQ, RR, SS,
     &          t , TT;
      INTEGER*4 sgn_ml, sgn_MLL;



      TalmiMoshinsky_2d = 0.D0;



      if( nr1.lt.0 .or. nr2.lt.0 .or. NRR.lt.0 .or. nr.lt.0 ) then
          return;
      endif
      if( ml1+ml2 .ne. MLL+ml ) then
          return;
      endif
      if( 2*nr1+abs(ml1) + 2*nr2+abs(ml2) .ne.
     &    2*NRR+abs(MLL) + 2*nr +abs(ml)        ) then
          return;
      endif
      Mmin = ( ml1+ml2 - (2*nr1+abs(ml1)+2*nr2+abs(ml2)) )/2;
      Mmax = ( ml1+ml2 + (2*nr1+abs(ml1)+2*nr2+abs(ml2)) )/2;
      if( MLL.lt.Mmin .or. MLL.gt.Mmax ) then
          return;
      endif
      Nrmax = ( 2*nr1+abs(ml1)+2*nr2+abs(ml2)-abs(MLL)-abs(ml) )/2;
      if( NRR.gt.Nrmax ) then
          return;
      endif



      call assert( nr1+nr2+NRR+nr.le.IGFV , 'IGFV too small' );
      call assert( nr1+abs(ml1)  .le.IGFV , 'IGFV too small' );
      call assert( nr2+abs(ml2)  .le.IGFV , 'IGFV too small' );
      call assert( NRR+abs(MLL)  .le.IGFV , 'IGFV too small' );
      call assert(  nr+abs(ml)   .le.IGFV , 'IGFV too small' );

      fac = iv(nr1+nr2+NRR+nr);
      fac = fac / DSQRT( 2.D0**(DBLE( 2*NRR+abs(MLL)+2*nr+abs(ml) )) );
      fac = fac * wf(nr1)*wf(nr1+abs(ml1))*wf(nr2)*wf(nr2+abs(ml2));
      fac = fac * wfi(NRR)*wfi(NRR+abs(MLL))*wfi(nr)*wfi(nr+abs(ml));

      sgn_ml  = 0;
      sgn_MLL = 0;
      if( ml  .ne. 0 ) sgn_ml  = isign(1,ml );
      if( MLL .ne. 0 ) sgn_MLL = isign(1,MLL);

      cond1 = (2*nr1+abs(ml1)) - (2*nr2+abs(ml2)) + abs(MLL) + abs(ml);
      cond2 = ml1;

      acc = 0.D0;
      do p = 0 , nr
      do q = 0 , nr-p
      do r = 0 , nr-p-q
      s = nr-p-q-r;

          mult1 = fak(nr)*fi(p)*fi(q)*fi(r)*fi(s);

          do PP = 0 , NRR
          do QQ = 0 , NRR-PP
          do RR = 0 , NRR-PP-QQ
          SS = NRR-PP-QQ-RR;

              mult2 = fak(NRR)*fi(PP)*fi(QQ)*fi(RR)*fi(SS);

              do t  = 0 , abs(ml)
              do TT = 0 , abs(MLL)

                  if( 2*(p+PP) - 2*(q+QQ) + 2*(t+TT) .ne. cond1 ) CYCLE;
                  if( r+RR-s-SS+sgn_ml*t+sgn_MLL*TT  .ne. cond2 ) CYCLE;


                  bin1 = fak(abs(ml) )*fi(t) *fi(abs(ml) -t );
                  bin2 = fak(abs(MLL))*fi(TT)*fi(abs(MLL)-TT);

                  acc = acc + iv(t+r+s)*mult1*mult2*bin1*bin2;

              enddo
              enddo

          enddo
          enddo
          enddo

      enddo
      enddo
      enddo

      TalmiMoshinsky_2d = fac * acc;

      return;
      end;






c======================================================================c

      subroutine CElliptic( epsi , xk , e1 , e2 , n )

c======================================================================c
      !******************************************************
      !* Complete elliptic integral of the first and second *
      !* kind. The input parameter is xk, which should be   *
      !* between 0 and 1. Technique uses Gauss' formula for *
      !* the arithmogeometrical mean. epsi is a measure of  *
      !* convergence accuracy. The returned values are e1,  *
      !* the elliptic integral of the first kind, and e2,   *
      !* the elliptic integral of the second kind.          *
      !* -------------------------------------------------- *
      !* Reference: Ball, algorithms for RPN calculators.   *
      !******************************************************
      IMPLICIT REAL*8    (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n)
      parameter( MAXN = 100 )

      REAL*8 epsi, x, xk, e1, e2, pi;
      REAL*8 A(0:MAXN), B(0:MAXN);
      INTEGER*4 j, m, n;


      pi = 3.14159265358979324D0;
      x  = DSQRT(xk);

      A(0) = 1.D0 + x ;
      B(0) = 1.D0 - x ;

      n = 0;
  100 n = n + 1;
      call assert( n.le.MAXN , 'Arrays too small in CElliptic()' );
      A(n) = ( A(n-1) + B(n-1) ) / 2.D0;
      B(n) = DSQRT( A(n-1) * B(n-1) );
      if( DABS( A(n) - B(n) ) > epsi ) GOTO 100;


      call assert( n.lt.30 , 'Possible overflow in CElliptic()' );

      e1 = pi/2.D0/A(n);
      e2 = 2.D0;
      m = 1;
      do j = 1 , n
          e2 = e2 - DBLE(m)*( A(j)*A(j) - B(j)*B(j) );
          m = m*2;
      enddo
      e2 = e2*e1/2.D0;

      return;
      end;
