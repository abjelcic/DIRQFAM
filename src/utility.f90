!======================================================================!

      subroutine assert( statement , error_message )

!======================================================================!
      implicit none;
      logical          , intent(in) :: statement;
      character(len=*) , intent(in) :: error_message;

          if( statement .eqv. .false. ) then
              write(6,'(a,a,a)') 'Error: [' , error_message,  '] !';
              stop 'Assertion failed, stopping the program.';
          end if

      return;
      end subroutine assert






!======================================================================!

      pure double precision function factorial(n) result(ans)

!======================================================================!
      implicit none;
      integer , intent(in) :: n;
      integer              :: i;

          ans = 1.d0;
          do i = 1 , n
              ans = ans * i;
          end do

          return;
      end function factorial






!======================================================================!

      double precision function tau3( it )

!======================================================================!
      implicit none;
      integer , intent(in) :: it;

          call assert( any(it==[1,2]) , 'tau3 function error.' );

          if( it == 1 ) tau3 = -1.d0; ! Neutron (coded as it=1) isospin convention.
          if( it == 2 ) tau3 = +1.d0; ! Proton  (coded as it=2) isospin convention.

      return;
      end function tau3






!======================================================================!

      integer function index_of_vector( fg , nz , nr , ml , ib )

!======================================================================!
      use simplex;
      implicit none;
      character , intent(in) :: fg;
      integer   , intent(in) :: nz;
      integer   , intent(in) :: nr;
      integer   , intent(in) :: ml;
      integer   , intent(in) :: ib;
      integer                :: i;
      logical                :: bool;

          bool = .false.;
          do i = 1 , id_spx(ib)
              if( ( fg == fg_spx(ib)%index(i) ) .and. &
                  ( nz == nz_spx(ib)%index(i) ) .and. &
                  ( nr == nr_spx(ib)%index(i) ) .and. &
                  ( ml == ml_spx(ib)%index(i) ) ) then
                  bool = .true.;
                  exit;
              end if
          end do

          call assert( bool .eqv. .true. , 'index_of_vector() wrong!' );

          index_of_vector = i;

      return;
      end function index_of_vector






!======================================================================!

      subroutine allocqfam()

!======================================================================!
      use v_matrix;
      use v_energy;
      use u_matrix;
      use u_energy;
      use fam_energies;
      use quadrature;
      use basis;
      use wbasis;
      use PHI;
      use gs_dens;
      use gs_mesons;
      use TempBlockMatrix;
      use fmatrix;
      use f02f20matrix;
      use spurious;
      use mesmat;
      use fam_green;
      use Wpairing;
      use dh;
      use dDelta;
      use ddens;
      use dh02dh20matrix;
      use dh11matrix;
      use xyfam;
      use drho;
      use dkappa;
      use dcurr;
      use dlaplace;
      use famlaplacianmod;
      use dmesons;
      use dcoulomb;
      use dpotentials;
      use nuclearlocfunc;
      use gmres;
      use KPMdata;
      implicit none;

          call alloc_v_matrix();
          call alloc_v_energy();
          call alloc_u_matrix();
          call alloc_u_energy();
          call alloc_fam_energies();
          call alloc_quadrature();
          call alloc_basis();
          call alloc_wbasis();
          call alloc_PHI();
          call alloc_gs_dens();
          call alloc_gs_mesons();
          call alloc_TempBlockMatrix();
          call alloc_fmatrix();
          call alloc_f02f20matrix();
          call alloc_spurious();
          call alloc_mesmat();
          call alloc_fam_green();
          call alloc_Wpairing();
          call alloc_dh();
          call alloc_dDelta();
          call alloc_ddens();
          call alloc_dh02dh20matrix();
          call alloc_dh11matrix();
          call alloc_xyfam();
          call alloc_drho();
          call alloc_dkappa();
          call alloc_dcurr();
          call alloc_dlaplace();
          call alloc_famlaplacianmod();
          call alloc_dmesons();
          call alloc_dcoulomb();
          call alloc_dpotentials();
          call alloc_nuclearlocfunc();
          call alloc_gmres();
          call alloc_KPMdata();

      return;
      end subroutine allocqfam






!======================================================================!

      subroutine deallocqfam()

!======================================================================!
      use v_matrix;
      use v_energy;
      use u_matrix;
      use u_energy;
      use fam_energies;
      use quadrature;
      use basis;
      use wbasis;
      use PHI;
      use gs_dens;
      use gs_mesons;
      use TempBlockMatrix;
      use fmatrix;
      use f02f20matrix;
      use spurious;
      use mesmat;
      use fam_green;
      use Wpairing;
      use dh;
      use dDelta;
      use ddens;
      use dh02dh20matrix;
      use dh11matrix;
      use xyfam;
      use drho;
      use dkappa;
      use dcurr;
      use dlaplace;
      use famlaplacianmod;
      use dmesons;
      use dcoulomb;
      use dpotentials;
      use nuclearlocfunc;
      use gmres;
      use KPMdata;
      implicit none;

          call dealloc_v_matrix();
          call dealloc_v_energy();
          call dealloc_u_matrix();
          call dealloc_u_energy();
          call dealloc_fam_energies();
          call dealloc_quadrature();
          call dealloc_basis();
          call dealloc_wbasis();
          call dealloc_PHI();
          call dealloc_gs_dens();
          call dealloc_gs_mesons();
          call dealloc_TempBlockMatrix();
          call dealloc_fmatrix();
          call dealloc_f02f20matrix();
          call dealloc_spurious();
          call dealloc_mesmat();
          call dealloc_fam_green();
          call dealloc_Wpairing();
          call dealloc_dh();
          call dealloc_dDelta();
          call dealloc_ddens();
          call dealloc_dh02dh20matrix();
          call dealloc_dh11matrix();
          call dealloc_xyfam();
          call dealloc_drho();
          call dealloc_dkappa();
          call dealloc_dcurr();
          call dealloc_dlaplace();
          call dealloc_famlaplacianmod();
          call dealloc_dmesons();
          call dealloc_dcoulomb();
          call dealloc_dpotentials();
          call dealloc_nuclearlocfunc();
          call dealloc_gmres();
          call dealloc_KPMdata();

      return;
      end subroutine deallocqfam






!======================================================================!

      double precision function phi_nz( nz , b , z )

!======================================================================!
!----------------------------------------------------------------------!
!                                                                      !
!  phi_nz(nz,b,z) =   1/sqrt(b) * 1/sqrt(sqrt(pi)*2^nz*nz!)            !
!                   * H_nz(z/b) * exp( -1/2 * (z/b)^2 )                !
!----------------------------------------------------------------------!
      implicit none;
      integer          , intent(in) :: nz;
      double precision , intent(in) :: b;
      double precision , intent(in) :: z;

      double precision , parameter  :: pi = 4.d0*atan(1.d0);
      double precision              :: x;
      double precision              :: phi_0, phi_1, tmp;
      integer                       :: n;


          call assert( nz >= 0    , 'nz <  0 in phi_nz()' );
          call assert( b  >  0.D0 , 'b  <= 0 in phi_nz()' );

          x = z/b;

          phi_0 = 1.d0/sqrt(b*sqrt(pi)) * exp( -0.5d0 * x**2 );
          phi_1 = x * sqrt(2.d0) * phi_0;

          if( nz == 0 ) then
              phi_nz = phi_0;
              return;
          end if

          if( nz == 1 ) then
              phi_nz = phi_1;
              return;
          end if


          do n = 2 , nz
              tmp   = phi_1;
              phi_1 = sqrt(2.d0/n)*x * phi_1 - sqrt((n-1.d0)/n) * phi_0;
              phi_0 = tmp;
          end do

          phi_nz = phi_1;


      return;
      end function phi_nz






!======================================================================!

      double precision function d_phi_nz( nz , b , z )

!======================================================================!
!----------------------------------------------------------------------!
!                                                                      !
!  d_phi_nz(nz,b,z) = d( phi_nz(nz,b,z) )/dz                           !
!                                                                      !
!----------------------------------------------------------------------!
      implicit none;
      integer          , intent(in) :: nz;
      double precision , intent(in) :: b;
      double precision , intent(in) :: z;

      double precision              :: x;
      double precision , external   :: phi_nz;


          call assert( nz >= 0    , 'nz <  0 in d_phi_nz()' );
          call assert( b  >  0.D0 , 'b  <= 0 in d_phi_nz()' );

          x = z/b;

          if( nz == 0 ) then
              d_phi_nz = -x/b * phi_nz(0,b,z);
              return;
          end if

          d_phi_nz = -(x/b) * phi_nz(nz,b,z) + sqrt(2.d0*nz)/b * phi_nz(nz-1,b,z);


      return;
      end function d_phi_nz






!======================================================================!

      double precision function phi_nr_ml( nr , ml , b , r )

!======================================================================!
!----------------------------------------------------------------------!
!                                                                      !
!  phi_nr_ml(nr,ml,b,r) =   1/b * sqrt( 2*(nr)! / (nr+|ml|)! )         !
!                         * (r/b)^|ml| * L_{nr}^|ml|( (r/b)^2 )        !
!                         *  exp( -1/2 * (r/b)^2 )                     !
!----------------------------------------------------------------------!
      implicit none;
      integer          , intent(in) :: nr;
      integer          , intent(in) :: ml;
      double precision , intent(in) :: b;
      double precision , intent(in) :: r;

      double precision , external   :: factorial;
      double precision              :: eta;
      double precision              :: phi_0_ml, phi_1_ml, tmp;
      integer                       :: n;


          call assert( nr >= 0    , 'nr <  0 in phi_nr_ml()' );
          call assert( b  >  0.d0 , 'b  <= 0 in phi_nr_ml()' );
          call assert( r  >  0.d0 , 'r  <= 0 in phi_nr_ml()' );

          eta = (r/b)**2;

          phi_0_ml = 1.d0/b * sqrt(2.d0)/sqrt(factorial(abs(ml))) * exp( -0.5d0 * ( eta - abs(ml)*log(eta) ) );
          phi_1_ml = phi_0_ml * ( -eta + (abs(ml)+1) ) / sqrt(abs(ml)+1.d0);

          if( nr == 0 ) then
              phi_nr_ml = phi_0_ml;
              return;
          end if

          if( nr == 1 ) then
              phi_nr_ml = phi_1_ml;
              return;
          end if

          do n = 2 , nr
              tmp      = phi_1_ml;
              phi_1_ml = (2*n-1+abs(ml)-eta)/sqrt(1.d0*n*(n+abs(ml))) * phi_1_ml - sqrt((n-1.d0)*(n-1.d0+abs(ml))/(n*(n+abs(ml)))) * phi_0_ml;
              phi_0_ml = tmp;
          end do

          phi_nr_ml = phi_1_ml;


      return;
      endfunction phi_nr_ml






!======================================================================!

      double precision function d_phi_nr_ml( nr , ml , b , r )

!======================================================================!
!----------------------------------------------------------------------!
!                                                                      !
!  d_phi_nr_ml(nr,ml,b,z) = d( phi_nr_ml(nr,ml,b,r) )/dr               !
!                                                                      !
!----------------------------------------------------------------------!
      implicit none;
      integer          , intent(in) :: nr;
      integer          , intent(in) :: ml;
      double precision , intent(in) :: b;
      double precision , intent(in) :: r;

      double precision              :: eta;
      double precision , external   :: phi_nr_ml;

          call assert( nr >= 0    , 'nr <  0 in d_phi_nr_ml()' );
          call assert( b  >  0.d0 , 'b  <= 0 in d_phi_nr_ml()' );
          call assert( r  >  0.d0 , 'r  <= 0 in d_phi_nr_ml()' );

          eta = (r/b)**2;

          if( nr == 0 ) then
              d_phi_nr_ml = (abs(ml)-eta)/(sqrt(eta)*b) * phi_nr_ml(0,abs(ml),b,r);
              return;
          end if

          d_phi_nr_ml = + (2*nr+abs(ml)-eta)/sqrt(eta*b**2)         * phi_nr_ml(nr  ,abs(ml),b,r) &
                        - 2.d0*sqrt( (nr*(nr+abs(ml)))/(eta*b**2) ) * phi_nr_ml(nr-1,abs(ml),b,r);


      return;
      end function d_phi_nr_ml






!======================================================================!

      double precision function I_K( a , K )

!======================================================================!
!
! Calculates I_K(a) = integral of sqrt(1-a*cos(x)^2)*cos(2*K*x) on [0,pi/2] interval.
! Parameter a is in [0,1].
! Parameter K is in {0,1,2,3,4,5,6}.
!
! There holds:
!     I_K(a) = ((4*K-4)/(2*K+1))*((2-a)/a) * I_{K-1}(a) - ((2*K-5)/(2*K+1)) * I_{K-2}(a),
!     I_0 = E(a),
!     I_1 = ((2-2*a)/(3*a)) * K(a) - ((2-a)/(3*a)) * E(a),
! where K(a) and E(a) are complete elliptic integrals of the first and second half respectively.
!
! Calculating I_K(a) via that recursion is inaccurate for small values of a,
! however for small values of a we can expand the integral I_K(a) around a=0.
!
! Taylor expansion of I_K(a) around a=0 is given by:
!     I_K(a)   = sum_{n=K}^{+oo} b_{n}(K) a^n,
!     b_{K}(K) = -pi/16^K * (1/(4*K-2)) * ( (2*K)! / (K!)^2 )
!     b_{n}(K) = ((2*n-1)*(2*n-3))/(4*(n^2-K^2)) * b_{n-1}(K), for n>K.
!
      implicit none;

      double precision ,                           intent(in) :: a;
      integer          ,                           intent(in) :: K;

      double precision ,                           external   :: factorial;
      double precision                                        :: Kc, Ec, x;
      double precision                                        :: I0, I1, Itmp;
      integer                                                 :: KK, n;
      double precision ,                           parameter  :: pi = 4.d0*atan(1.d0);
      integer          ,                           parameter  :: Kmax = 6;
      integer          ,                           parameter  :: LMT = 64;
      double precision , dimension(0:LMT,0:Kmax) , save       :: b;
      logical          ,                           save       :: first_call = .true.;


          call assert( a>=0.d0 .and. a<=1.d0 , 'a not in [0,1]'           );
          call assert( K>=0    .and. K<=Kmax , 'K not in {0,1,2,3,4,5,6}' );

          if( first_call ) then
              first_call = .false.;

              do KK = 0 , Kmax

                  b( 0:KK-1 , KK ) = 0.d0;
                  b(   KK   , KK ) = (-pi/(16**KK)) * (1.d0/(4*KK-2)) * ( factorial(2*KK) / factorial(KK)**2 );

                  do n = KK+1 , LMT
                      b(n,KK) = b(n-1,KK) * ((2*n-1)*(2*n-3))/(4.d0*(n**2-KK**2));
                  end do

              end do

          end if


          if( a > 0.5d0 ) then

              call CElliptic( a , Kc , Ec );
              ! CElliptic is inaccurate if a=1 is given, however one
              ! can easily see that there holds: I_K(1) = -1/(4*K^2-1).
              if( abs(a-1) < 1.d-15  ) then
                  I_K = -1.d0/(4*K**2-1);
                  return;
              end if


              I0 = Ec;
              I1 = ((2-2*a)/(3*a))*Kc - ((2-a)/(3*a))*Ec;

              if( K == 0 ) then
                  I_K = I0;
                  return;
              end if

              if( K == 1 ) then
                  I_K = I1;
                  return;
              end if

              do KK = 2 , K
                  Itmp = I1;
                  I1   = ((4*KK-4)*(2-a)/((2*KK+1)*a))*I1 - ((2.d0*KK-5)/(2*KK+1))*I0;
                  I0   = Itmp;
              end do
              I_K = I1;
              return;

          else

              I_K = 0.d0;
              x   = 1.d0;
              do n = 0 , LMT
                  I_K = I_K + b(n,K) * x;
                  x   = x*a;
              end do
              return;

          end if


          return;

      contains

          subroutine CElliptic( x , K , E )
          ! Subroutine uses Gauss's formula for the arithmogeometrical mean.
          ! Reference: Ball, Algorithms for RPN calculators.
          implicit none;
          double precision , intent(in)       :: x; ! Must be in [0,1).
          double precision , intent(out)      :: K; ! Returns K(x).
          double precision , intent(out)      :: E; ! Returns E(x).

          integer                             :: n, i;
          double precision , parameter        :: pi   = 4.d0*atan(1.d0);
          double precision , parameter        :: epsi = 1.d-12;
          integer          , parameter        :: LMT  = 2048;
          double precision , dimension(0:LMT) :: A;
          double precision , dimension(0:LMT) :: B;


              A(0) = 1.d0 + sqrt(x);
              B(0) = 1.d0 - sqrt(x);

              do n = 1 , LMT
                  A(n) = 0.5d0 * ( A(n-1) + B(n-1) );
                  B(n) = sqrt(abs( A(n-1) * B(n-1) ));

                  if( abs(A(n)-B(n)) < epsi ) exit;
              end do
              if( n == LMT+1 ) write(*,'(a)') 'Possible accuracy problem in CElliptic() subroutine!';


              K = pi/2 * 1/A(n);
              E = 2.d0;
              do i = 1 , n
                  E = E + 2.d0**(i-1) * ( B(i)**2 - A(i)**2 );
              end do
              E = K * E / 2;


              return;
          end subroutine CElliptic

      end function I_K






!======================================================================!

      double precision function TalmiMoshinsky_z( nz1 , nz2 , N_z , nz )

!======================================================================!
      implicit none;
      integer          ,                    intent(in) :: nz1;
      integer          ,                    intent(in) :: nz2;
      integer          ,                    intent(in) :: N_z;
      integer          ,                    intent(in) :: nz;

      integer                                          :: n, p;
      double precision                                 :: bin1, bin2;
      integer          ,                    parameter  :: LMT = 127;
      double precision ,                    external   :: factorial;
      double precision , dimension(0:LMT) , save       :: fact;
      double precision , dimension(0:LMT) , save       :: factinv;
      logical          ,                    save       :: first_call = .true.;


          if( first_call ) then
              first_call = .false.;
              do n = 0 , LMT
                  fact(n)    = factorial(n);
                  factinv(n) = 1.d0 / factorial(n);
              end do
          end if

          ! Implementation of Eq. (F.2).
          TalmiMoshinsky_z = 0.d0;

          if( .not.( nz1>=0 .and. nz2>=0 .and. N_z>=0 .and. nz>=0 .and. nz1+nz2==N_z+nz ) ) return;

          call assert( nz1 <= LMT , 'LMT too small in TalmiMoshinsky_z' );
          call assert( nz2 <= LMT , 'LMT too small in TalmiMoshinsky_z' );
          call assert( N_z <= LMT , 'LMT too small in TalmiMoshinsky_z' );
          call assert( nz  <= LMT , 'LMT too small in TalmiMoshinsky_z' );

          do p = max(0,nz1-N_z) , min(nz,nz1)

              bin1 = fact(nz)  * factinv(p)     * factinv(nz-p);
              bin2 = fact(N_z) * factinv(nz1-p) * factinv(N_z-nz1+p);

              TalmiMoshinsky_z = TalmiMoshinsky_z + (-1)**p * bin1 * bin2;
          end do

          TalmiMoshinsky_z = sqrt(fact(nz1)*fact(nz2)*factinv(N_z)*factinv(nz)) / sqrt(2.d0**(nz1+nz2)) * TalmiMoshinsky_z;


          return;
      end function TalmiMoshinsky_z






!======================================================================!

      double precision function TalmiMoshinsky_r( nr1 , ml1 , nr2 , ml2 , &
                                                  N_r , M_l , nr  , ml    )

!======================================================================!
      implicit none;
      integer          ,                    intent(in) :: nr1;
      integer          ,                    intent(in) :: ml1;
      integer          ,                    intent(in) :: nr2;
      integer          ,                    intent(in) :: ml2;
      integer          ,                    intent(in) :: N_r;
      integer          ,                    intent(in) :: M_l;
      integer          ,                    intent(in) :: nr;
      integer          ,                    intent(in) :: ml;

      integer                                          :: n, p, q, r, s, PP, QQ, RR, SS, t, TT;
      integer                                          :: Mmin, Mmax, Nrmax;
      integer                                          :: conditionC3, conditionC4;
      double precision                                 :: mult1, mult2, bin1, bin2;
      integer          ,                    parameter  :: LMT = 127;
      double precision ,                    external   :: factorial;
      double precision , dimension(0:LMT) , save       :: fact;
      double precision , dimension(0:LMT) , save       :: factinv;
      logical          ,                    save       :: first_call = .true.;


          if( first_call ) then
              first_call = .false.;
              do n = 0 , LMT
                  fact(n)    = factorial(n);
                  factinv(n) = 1.d0 / factorial(n);
              end do
          end if

          ! Implementation of Eq. (F.6).
          TalmiMoshinsky_r = 0.D0;

          Mmin  = ( ml1+ml2 - ( 2*nr1+abs(ml1) + 2*nr2+abs(ml2) ) )/2;
          Mmax  = ( ml1+ml2 + ( 2*nr1+abs(ml1) + 2*nr2+abs(ml2) ) )/2;
          Nrmax = ( 2*nr1+abs(ml1) + 2*nr2+abs(ml2) - abs(M_l) - abs(ml) )/2;

          if( .not.( nr1>=0 .and. nr2>=0 .and. N_r>=0 .and. nr>=0                     ) ) return;
          if( .not.( ml1+ml2 == M_l+ml                                                ) ) return;
          if( .not.( 2*nr1+abs(ml1) + 2*nr2+abs(ml2) == 2*N_r+abs(M_l) + 2*nr+abs(ml) ) ) return;
          if( .not.( M_l>=Mmin .and. M_l<=Mmax                                        ) ) return;
          if( .not.( N_r <= Nrmax                                                     ) ) return;

          call assert( nr1 + abs(ml1) <= LMT , 'LMT too small in TalmiMoshinsky_r' );
          call assert( nr2 + abs(ml2) <= LMT , 'LMT too small in TalmiMoshinsky_r' );
          call assert( N_r + abs(M_l) <= LMT , 'LMT too small in TalmiMoshinsky_r' );
          call assert( nr  + abs(ml)  <= LMT , 'LMT too small in TalmiMoshinsky_r' );


          conditionC3 = (2*nr1+abs(ml1)) - (2*nr2+abs(ml2)) + abs(M_l) + abs(ml);
          conditionC4 = ml1;

          do p = 0 , nr
          do q = 0 , nr-p
          do r = 0 , nr-p-q
              s = nr-p-q-r; ! Condition (C1).

              mult1 = fact(nr) * factinv(p)*factinv(q)*factinv(r)*factinv(s);

              do PP = 0 , N_r
              do QQ = 0 , N_r-PP
              do RR = 0 , N_r-PP-QQ
                  SS = N_r-PP-QQ-RR; ! Condition (C2).

                  mult2 = fact(N_r) * factinv(PP)*factinv(QQ)*factinv(RR)*factinv(SS);

                  do t  = 0 , abs(ml)
                  do TT = 0 , abs(M_l)

                      if( 2*(p+PP) - 2*(q+QQ) + 2*(t+TT) == conditionC3 ) then                      ! Condition (C3).
                          if( (r+RR) - (s+SS) + signum(ml)*t + signum(M_l)*TT == conditionC4 ) then ! Condition (C4).

                              bin1 = fact(abs(ml))  * factinv(t)  * factinv(abs(ml)-t);
                              bin2 = fact(abs(M_l)) * factinv(TT) * factinv(abs(M_l)-TT);

                              TalmiMoshinsky_r = TalmiMoshinsky_r + (-1)**(t+r+s) * mult1 * mult2 * bin1 * bin2;

                          end if
                      end if

                  end do
                  end do

              end do
              end do
              end do

          end do
          end do
          end do

          TalmiMoshinsky_r = TalmiMoshinsky_r * (-1)**(nr1+nr2+N_r+nr)                                                   &
                                              * 1.d0 / sqrt( 2.d0**( 2*N_r+abs(M_l) + 2*nr+abs(ml) ) )                   &
                                              * sqrt( fact(nr1) * fact(nr1+abs(ml1)) * fact(nr2) * fact(nr2+abs(ml2)) )  &
                                              / sqrt( fact(N_r) * fact(N_r+abs(M_l)) * fact(nr)  * fact(nr+abs(ml))   );


          return;

      contains

          pure integer function signum(n) result(ans)
              implicit none;
              integer , intent(in) :: n;

                  ans = 0;
                  if( n /= 0 ) ans = sign(1,n);

              return;
          end function signum

      end function TalmiMoshinsky_r






!======================================================================!

      module module_SVD

!======================================================================!
          implicit none;

          public :: calculateSVD;
          public :: getSVDerror;

          private :: allocate_1Dreal;
          private :: allocate_2Dreal;
          private :: deallocate_1Dreal;
          private :: deallocate_2Dreal;

      contains

          subroutine calculateSVD( A , rank , U , V )
              implicit none;

              double precision , dimension(:,:) ,               intent(in)  :: A;
              integer          ,                                intent(in)  :: rank;
              double precision , dimension(:,:) , allocatable , intent(out) :: U;
              double precision , dimension(:,:) , allocatable , intent(out) :: V;

              double precision , dimension(:,:) , allocatable               :: Acpy;
              integer                                                       :: i, j;
              integer                                                       :: m, n;
              integer                                                       :: mx, mn;
              integer                                                       :: INFO;
              double precision , dimension(:)   , allocatable               :: SingVals;
              double precision , dimension(:,:) , allocatable               :: Ufull, Vfull;
              double precision , dimension(1)                               :: WORKquery;
              double precision , dimension(:)   , allocatable               :: WORK;
              external                                                      :: dgesvd;


                  m = size(A,1);
                  n = size(A,2);

                  mx = max(m,n);
                  mn = min(m,n);

                  if( rank > mn ) then
                      stop 'Requested rank in calculateSVD is larger than min(m,n).';
                  end if

                  if( mn == 0 ) then
                      call allocate_2Dreal( U , [0,0] );
                      call allocate_2Dreal( V , [0,0] );
                      return;
                  end if

                  call allocate_2Dreal( Acpy     , [ m , n ]        );
                  call allocate_1Dreal( SingVals , [ min(m,n) ]     );
                  call allocate_2Dreal( Ufull    , [ m , min(m,n) ] );
                  call allocate_2Dreal( Vfull    , [ min(m,n) , n ] );


                  Acpy(1:m,1:n) = A(1:m,1:n);

                  call dgesvd( 'S'       , 'S'           , &
                               m         , n             , &
                               Acpy      , size(Acpy,1)  , &
                               SingVals  ,                 &
                               Ufull     , size(Ufull,1) , &
                               Vfull     , size(Vfull,1) , &
                               WORKquery , -1            , &
                               INFO                        );
                  if( INFO/=0 ) stop 'INFO/=0 in dgesvd workspace query.';

                  call allocate_1Dreal( WORK , [ int(WORKquery(1)+0.5d0) ] );

                  call dgesvd( 'S'       , 'S'           , &
                               m         , n             , &
                               Acpy      , size(Acpy,1)  , &
                               SingVals  ,                 &
                               Ufull     , size(Ufull,1) , &
                               Vfull     , size(Vfull,1) , &
                               WORK      , size(WORK)    , &
                               INFO                        );
                  if( INFO/=0 ) stop 'INFO/=0 in dgesvd!';



                  call allocate_2Dreal( U , [ m , rank ] );
                  call allocate_2Dreal( V , [ rank , n ] );

                  if( any( SingVals < 0.0 ) ) stop 'dgesvd calculated negative sungular value!';

                  do j = 1 , rank
                      do i = 1 , m
                          U(i,j) = sqrt(SingVals(j)) * Ufull(i,j);
                      end do
                  end do

                  do j = 1 , n
                      do i = 1 , rank
                          V(i,j) = sqrt(SingVals(i)) * Vfull(i,j);
                      end do
                  end do




                  call deallocate_2Dreal( Acpy     );
                  call deallocate_1Dreal( SingVals );
                  call deallocate_2Dreal( Ufull    );
                  call deallocate_2Dreal( Vfull    );
                  call deallocate_1Dreal( WORK     );


              return;
          end subroutine calculateSVD

          double precision function getSVDerror( E , U , V ) result(relerr) ! relerr = || U*V - E ||_F / || E ||_F.
              implicit none;

              double precision , dimension(:,:) , intent(in)  :: E;
              double precision , dimension(:,:) , intent(in)  :: U;
              double precision , dimension(:,:) , intent(in)  :: V;

              double precision , dimension(:,:) , allocatable :: Eapprox;
              external                                        :: dgemm;


                  if( .not. ( size(E,1)==size(U,1) .and. size(U,2)==size(V,1) .and. size(V,2)==size(E,2) ) ) stop 'Incompatible sizes in getLowRankApproxError.';

                  call allocate_2Dreal( Eapprox , [size(E,1),size(E,2)] );

                  call dgemm( 'N' , 'N' , size(U,1) , size(V,2) , size(U,2) , 1.d0 , U,size(U,1) , V,size(V,1) , 0.d0 , Eapprox,size(Eapprox,1) );

                  relerr =  sqrt(sum(abs(Eapprox-E)**2)) / sqrt(sum(abs(E)**2));

                  call deallocate_2Dreal( Eapprox );


              return;
          end function getSVDerror

          subroutine allocate_1Dreal( A , size , Message )
              implicit none;
              double precision , dimension(:) , allocatable , intent(inout) :: A;
              integer          , dimension(1) ,               intent(in)    :: size;
              character(len=*) ,                optional    , intent(in)    :: Message;
              integer                                                       :: n;
              integer                                                       :: allocStat, deallocStat;


                  if( allocated(A) ) then
                      deallocate( A , stat=deallocStat );
                      if( deallocStat /= 0 ) then
                          if( present(Message) ) then
                              write(*,'(a,a,a,i0,a)') 'Allocation failed [' // Message // '] with status code: ' , deallocStat , '.';
                          else
                              write(*,'(a,i0,a)') 'Allocation failed with status code: ' , deallocStat , '.';
                          end if
                          stop;
                      end if
                  end if

                  n = size(1);

                  allocate( A( 1:n ) , stat=allocStat );
                  if( allocStat /= 0 ) then
                      if( present(Message) ) then
                          write(*,'(a,a,a,i0,a)') 'Allocation failed [' // Message // '] with status code: ' , allocStat , '.';
                      else
                          write(*,'(a,i0,a)') 'Allocation failed with status code: ' , allocStat , '.';
                      end if
                      stop;
                  end if


              return;
          end subroutine allocate_1Dreal

          subroutine allocate_2Dreal( A , size , Message )
              implicit none;
              double precision , dimension(:,:) , allocatable , intent(inout) :: A;
              integer          , dimension(2)   ,               intent(in)    :: size;
              character(len=*) ,                  optional    , intent(in)    :: Message;
              integer                                                         :: m;
              integer                                                         :: n;
              integer                                                         :: allocStat, deallocStat;


                  if( allocated(A) ) then
                      deallocate( A , stat=deallocStat );
                      if( deallocStat /= 0 ) then
                          if( present(Message) ) then
                              write(*,'(a,a,a,i0,a)') 'Allocation failed [' // Message // '] with status code: ' , deallocStat , '.';
                          else
                              write(*,'(a,i0,a)') 'Allocation failed with status code: ' , deallocStat , '.';
                          end if
                          stop;
                      end if
                  end if

                  m = size(1);
                  n = size(2);

                  allocate( A( 1:m , 1:n ) , stat=allocStat );
                  if( allocStat /= 0 ) then
                      if( present(Message) ) then
                          write(*,'(a,a,a,i0,a)') 'Allocation failed [' // Message // '] with status code: ' , allocStat , '.';
                      else
                          write(*,'(a,i0,a)') 'Allocation failed with status code: ' , allocStat , '.';
                      end if
                      stop;
                  end if


              return;
          end subroutine allocate_2Dreal

          subroutine deallocate_1Dreal( A , Message )
              implicit none;
              double precision , dimension(:) , allocatable , intent(inout) :: A;
              character(len=*) ,                optional    , intent(in)    :: Message;
              integer                                                       :: deallocStat;


                  if( allocated(A) ) then
                      deallocate( A , stat=deallocStat );
                      if( deallocStat /= 0 ) then
                          if( present(Message) ) then
                              write(*,'(a,a,a,i0,a)') 'Deallocation failed [' // Message // '] with status code: ' , deallocStat , '.';
                          else
                              write(*,'(a,i0,a)') 'Deallocation failed with status code: ' , deallocStat , '.';
                          end if
                          stop;
                      end if
                  end if


              return;
          end subroutine deallocate_1Dreal

          subroutine deallocate_2Dreal( A , Message )
              implicit none;
              double precision , dimension(:,:) , allocatable , intent(inout) :: A;
              character(len=*) ,                  optional    , intent(in)    :: Message;
              integer                                                         :: deallocStat;


                  if( allocated(A) ) then
                      deallocate( A , stat=deallocStat );
                      if( deallocStat /= 0 ) then
                          if( present(Message) ) then
                              write(*,'(a,a,a,i0,a)') 'Deallocation failed [' // Message // '] with status code: ' , deallocStat , '.';
                          else
                              write(*,'(a,i0,a)') 'Deallocation failed with status code: ' , deallocStat , '.';
                          end if
                          stop;
                      end if
                  end if


              return;
          end subroutine deallocate_2Dreal

      end module module_SVD






!======================================================================!

      subroutine transposeBlockMatrix( A )

!======================================================================!
          use dataTypes;
          implicit none;
          type(complexBlockMatrix) , intent(inout) :: A;
          integer                                  :: ib1, ib2, i1, i2;
          double complex                           :: tmp;

              do ib2 = 1 , size(A%nnzblocks,2)
                  do ib1 = 1 , ib2-1
                      if( A%nnzblocks(ib1,ib2) ) then
                          do i2 = 1 , size(A%blocks(ib1,ib2)%mat,2)
                              do i1 = 1 , size(A%blocks(ib1,ib2)%mat,1)
                                  tmp                          = A%blocks(ib1,ib2)%mat(i1,i2)
                                  A%blocks(ib1,ib2)%mat(i1,i2) = A%blocks(ib2,ib1)%mat(i2,i1);
                                  A%blocks(ib2,ib1)%mat(i2,i1) = tmp;
                              end do
                          end do
                      end if
                  end do
                  if( A%nnzblocks(ib2,ib2) ) then
                      do i2 = 1 , size(A%blocks(ib2,ib2)%mat,2)
                          do i1 = 1 , i2
                              tmp                          = A%blocks(ib2,ib2)%mat(i1,i2);
                              A%blocks(ib2,ib2)%mat(i1,i2) = A%blocks(ib2,ib2)%mat(i2,i1);
                              A%blocks(ib2,ib2)%mat(i2,i1) = tmp;
                          end do
                      end do
                  end if
              end do

          return;
      end subroutine transposeBlockMatrix






!======================================================================!

      subroutine hermitianTransposeBlockMatrix( A )

!======================================================================!
          use dataTypes;
          implicit none;
          type(complexBlockMatrix) , intent(inout) :: A;
          integer                                  :: ib1, ib2, i1, i2;
          double complex                           :: tmp;

              do ib2 = 1 , size(A%nnzblocks,2)
                  do ib1 = 1 , ib2-1
                      if( A%nnzblocks(ib1,ib2) ) then
                          do i2 = 1 , size(A%blocks(ib1,ib2)%mat,2)
                              do i1 = 1 , size(A%blocks(ib1,ib2)%mat,1)
                                  tmp                          = A%blocks(ib1,ib2)%mat(i1,i2)
                                  A%blocks(ib1,ib2)%mat(i1,i2) = conjg(A%blocks(ib2,ib1)%mat(i2,i1));
                                  A%blocks(ib2,ib1)%mat(i2,i1) = conjg(tmp);
                              end do
                          end do
                      end if
                  end do
                  if( A%nnzblocks(ib2,ib2) ) then
                      do i2 = 1 , size(A%blocks(ib2,ib2)%mat,2)
                          do i1 = 1 , i2
                              tmp                          = A%blocks(ib2,ib2)%mat(i1,i2);
                              A%blocks(ib2,ib2)%mat(i1,i2) = conjg(A%blocks(ib2,ib2)%mat(i2,i1));
                              A%blocks(ib2,ib2)%mat(i2,i1) = conjg(tmp);
                          end do
                      end do
                  end if
              end do

          return;
      end subroutine hermitianTransposeBlockMatrix






!======================================================================!

      subroutine setToZeroBlockMatrix( A )

!======================================================================!
          use dataTypes;
          implicit none;
          type(complexBlockMatrix) , intent(inout) :: A;
          integer                                  :: ib1, ib2, i1, i2;

              do ib2 = 1 , size(A%blocks,2)
                  do ib1 = 1 , size(A%blocks,1)
                      if( A%nnzblocks(ib1,ib2) ) then
                          do i2 = 1 , size(A%blocks(ib1,ib2)%mat,2)
                              do i1 = 1 , size(A%blocks(ib1,ib2)%mat,1)
                                  A%blocks(ib1,ib2)%mat(i1,i2) = cmplx(0.d0,0.d0,kind=8);
                              end do
                          end do
                      end if
                  end do
              end do

          return;
      end subroutine setToZeroBlockMatrix






!======================================================================!

      double precision function frobeniousNormOfBlockMatrix( A ) result(ans)

!======================================================================!
          use dataTypes;
          implicit none;
          type(complexBlockMatrix) , intent(in) :: A;
          integer                               :: ib1, ib2, i1, i2;

              ans = 0.d0;
              do ib2 = 1 , size(A%blocks,2)
                  do ib1 = 1 , size(A%blocks,1)
                      if( A%nnzblocks(ib1,ib2) ) then
                          do i2 = 1 , size(A%blocks(ib1,ib2)%mat,2)
                              do i1 = 1 , size(A%blocks(ib1,ib2)%mat,1)
                                  ans = ans + abs(A%blocks(ib1,ib2)%mat(i1,i2))**2;
                              end do
                          end do
                      end if
                  end do
              end do
              ans = sqrt(ans);

          return;
      end function frobeniousNormOfBlockMatrix
