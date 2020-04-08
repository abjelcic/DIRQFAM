      module blokap
      implicit none;

      !common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
      integer nb;
      integer,     dimension(:), allocatable :: kb;
      integer,     dimension(:), allocatable :: mb;
      character*6, dimension(:), allocatable :: tb;

      contains
          subroutine alloc_blokap()
          use dirqfampar
          implicit none;

          allocate( kb(NBX) );
          allocate( mb(NBX) );
          allocate( tb(NBX) );

          end

          subroutine dealloc_blokap()
          implicit none;

          deallocate( kb );
          deallocate( mb );
          deallocate( tb );

          end

      end module blokap

      module bloosc
      implicit none;

      !common /bloosc/ ia(NBX,2),id(NBX,2)
      integer, dimension(:,:), allocatable :: ia;
      integer, dimension(:,:), allocatable :: id;

      contains
          subroutine alloc_bloosc()
          use dirqfampar
          implicit none;

          allocate( ia(NBX,2) );
          allocate( id(NBX,2) );

          end

          subroutine dealloc_bloosc()
          implicit none;

          deallocate( ia );
          deallocate( id );

          end

      end module bloosc

      module quaosc
      implicit none;

      !common /quaosc/ nt,nz(NTX),nr(NTX),ml(NTX),ms(NTX),np(NTX),tt(NTX)
      integer nt;
      integer,     dimension(:), allocatable :: nz;
      integer,     dimension(:), allocatable :: nr;
      integer,     dimension(:), allocatable :: ml;
      integer,     dimension(:), allocatable :: ms;
      integer,     dimension(:), allocatable :: np;
      character*8, dimension(:), allocatable :: tt;

      contains
          subroutine alloc_quaosc()
          use dirqfampar
          implicit none;

          allocate( nz(NTX) );
          allocate( nr(NTX) );
          allocate( ml(NTX) );
          allocate( ms(NTX) );
          allocate( np(NTX) );
          allocate( tt(NTX) );

          end

          subroutine dealloc_quaosc()
          implicit none;

          deallocate( nz );
          deallocate( nr );
          deallocate( ml );
          deallocate( ms );
          deallocate( np );
          deallocate( tt );

          end

      end module quaosc

      module bosqua
      implicit none;

      !common /bosqua/ nzb(NOX),nrb(NOX),NO
      integer NO;
      integer, dimension(:), allocatable :: nzb;
      integer, dimension(:), allocatable :: nrb;

      contains
          subroutine alloc_bosqua()
          use dirqfampar
          implicit none;

          allocate( nzb(NOX) );
          allocate( nrb(NOX) );

          end

          subroutine dealloc_bosqua()
          implicit none;

          deallocate( nzb );
          deallocate( nrb );

          end

      end module bosqua

      module vvvikf
      implicit none;

      !common /vvvikf/ mv,ipos(NBX),nib(MVX),nni(2,MVX)
      integer mv;
      integer, dimension(:),   allocatable :: ipos;
      integer, dimension(:),   allocatable :: nib;
      integer, dimension(:,:), allocatable :: nni;

      contains
          subroutine alloc_vvvikf()
          use dirqfampar
          implicit none;

          allocate( ipos(NBX)  );
          allocate( nib(MVX)   );
          allocate( nni(2,MVX) );

          end

          subroutine dealloc_vvvikf()
          implicit none;

          deallocate( ipos );
          deallocate( nib  );
          deallocate( nni  );

          end

      end module vvvikf

      module deldel
      implicit none;

      !common /deldel/ de(NHHX,NB2X)
      double precision, dimension(:,:), allocatable :: de;

      contains
          subroutine alloc_deldel()
          use dirqfampar
          implicit none;

          allocate( de(NHHX,NB2X)  );

          end

          subroutine dealloc_deldel()
          implicit none;

          deallocate( de );

          end

      end module deldel

      module gamgam
      implicit none;

      !common /gamgam/ hh(NHHX,NB2X)
      double precision, dimension(:,:), allocatable :: hh;

      contains
          subroutine alloc_gamgam()
          use dirqfampar
          implicit none;

          allocate( hh(NHHX,NB2X)  );

          end

          subroutine dealloc_gamgam()
          implicit none;

          deallocate( hh );

          end

      end module gamgam

      module broyde
      implicit none;

      !parameter (nn = 2*MVTX+2*MVX+2+1)
      !parameter (mm = 7)
      !common /broyde/ bbeta(mm,mm),df(nn,mm),dv(nn,mm), bwork(mm),curv(nn),bw0,ibwork(mm)
      !common /broyde1/ vin(nn)
      !common /broyde2/ ibroyd
      !dimension vou(nn)
      integer nn;
      integer mm;
      integer ibroyd;
      double precision bw0;
      double precision, dimension(:,:), allocatable :: bbeta;
      double precision, dimension(:,:), allocatable :: df;
      double precision, dimension(:,:), allocatable :: dv;
      double precision, dimension(:),   allocatable :: bwork;
      double precision, dimension(:),   allocatable :: curv;
      double precision, dimension(:),   allocatable :: vin;
      double precision, dimension(:),   allocatable :: vou;
      integer         , dimension(:),   allocatable :: ibwork;

      contains
          subroutine alloc_broyde()
          use dirqfampar
          implicit none;

          nn = 2*MVTX + 2*MVX + 2 + 1;
          mm = 7;

          allocate( bbeta(mm,mm) );
          allocate(    df(nn,mm) );
          allocate(    dv(nn,mm) );
          allocate(    bwork(mm) );
          allocate(     curv(nn) );
          allocate(      vin(nn) );
          allocate(      vou(nn) );
          allocate(   ibwork(mm) );

          end

          subroutine dealloc_broyde()
          implicit none;

          deallocate( bbeta  );
          deallocate( df     );
          deallocate( dv     );
          deallocate( bwork  );
          deallocate( curv   );
          deallocate( vin    );
          deallocate( vou    );
          deallocate( ibwork );

          end

      end module broyde

      module canonmod
      implicit none;

      !dimension aa(NHHX),dd(NHHX),v2(NHX),z(NHX),eb(NHX),h(NHX),d(NHX)
      double precision, dimension(:), allocatable :: aa;
      double precision, dimension(:), allocatable :: dd;
      double precision, dimension(:), allocatable :: v2;
      double precision, dimension(:), allocatable :: z;
      double precision, dimension(:), allocatable :: eb;
      double precision, dimension(:), allocatable :: h;
      double precision, dimension(:), allocatable :: d;

      contains
          subroutine alloc_canonmod()
          use dirqfampar
          implicit none;

          allocate( aa(NHHX) );
          allocate( dd(NHHX) );
          allocate( v2(NHX ) );
          allocate(  z(NHX ) );
          allocate( eb(NHX ) );
          allocate(  h(NHX ) );
          allocate(  d(NHX ) );

          end

          subroutine dealloc_canonmod()
          implicit none;

          deallocate( aa );
          deallocate( dd );
          deallocate( v2 );
          deallocate( z  );
          deallocate( eb );
          deallocate( h  );
          deallocate( d  );

          end

      end module canonmod

      module blodir
      implicit none;

      !common /blodir/ ka(NBX,4),kd(NBX,4)
      integer, dimension(:,:), allocatable :: ka;
      integer, dimension(:,:), allocatable :: kd;

      contains
          subroutine alloc_blodir()
          use dirqfampar
          implicit none;

          allocate( ka(NBX,4) );
          allocate( kd(NBX,4) );

          end

          subroutine dealloc_blodir()
          implicit none;

          deallocate( ka );
          deallocate( kd );

          end

      end module blodir

      module eeecan
      implicit none;

      !common /eeecan/ eecan(KX,4),decan(KX,4),vvcan(KX,4),fgcan(NHX,KX,4),ibkcan(KX,4)
      double precision, dimension(:,:),   allocatable :: eecan;
      double precision, dimension(:,:),   allocatable :: decan;
      double precision, dimension(:,:),   allocatable :: vvcan;
      double precision, dimension(:,:,:), allocatable :: fgcan;
      integer,          dimension(:,:),   allocatable :: ibkcan;

      contains
          subroutine alloc_eeecan()
          use dirqfampar
          implicit none;

          allocate(     eecan(KX,4) );
          allocate(     decan(KX,4) );
          allocate(     vvcan(KX,4) );
          allocate( fgcan(NHX,KX,4) );
          allocate(    ibkcan(KX,4) );

          end

          subroutine dealloc_eeecan()
          implicit none;

          deallocate( eecan  );
          deallocate( decan  );
          deallocate( vvcan  );
          deallocate( fgcan  );
          deallocate( ibkcan );

          end

      end module eeecan

      module blocan
      implicit none;

      !common /blocan/ kacan(nbx,4),kdcan(nbx,4),nkcan(4)
      integer, dimension(:,:), allocatable :: kacan;
      integer, dimension(:,:), allocatable :: kdcan;
      integer, dimension(:),   allocatable :: nkcan;

      contains
          subroutine alloc_blocan()
          use dirqfampar
          implicit none;

          allocate( kacan(NBX,4) );
          allocate( kdcan(NBX,4) );
          allocate(     nkcan(4) );

          end

          subroutine dealloc_blocan()
          implicit none;

          deallocate( kacan );
          deallocate( kdcan );
          deallocate( nkcan );

          end

      end module blocan

      module waveuv
      implicit none;

      !common /waveuv/ fguv(NHBX,KX,4),equ(KX,4)
      double precision, dimension(:,:,:), allocatable :: fguv;
      double precision, dimension(:,:),   allocatable :: equ;

      contains
          subroutine alloc_waveuv()
          use dirqfampar
          implicit none;

          allocate( fguv(NHBX,KX,4) );
          allocate(       equ(KX,4) );

          end

          subroutine dealloc_waveuv()
          implicit none;

          deallocate( fguv );
          deallocate( equ  );

          end

      end module waveuv

      module centmasmod
      implicit none;

      !parameter (MG4 = 4*MG)
      !dimension wc(MG4,KX),dz(MG4,KX),dr(MG4,KX)
      integer MG4;
      double precision, dimension(:,:), allocatable :: wc;
      double precision, dimension(:,:), allocatable :: dz;
      double precision, dimension(:,:), allocatable :: dr;

      contains
          subroutine alloc_centmasmod()
          use dirqfampar
          implicit none;

          MG4 = 4*MG;
          allocate( wc(MG4,KX) );
          allocate( dz(MG4,KX) );
          allocate( dr(MG4,KX) );

          end

          subroutine dealloc_centmasmod()
          implicit none;

          deallocate( wc );
          deallocate( dz );
          deallocate( dr );

          end

      end module centmasmod

      module gaussl
      implicit none;

      !common /gaussl/ xl(0:NGL),wl(0:NGL),sxl(0:NGL),rb(0:NGL)
      double precision, dimension(:), allocatable :: xl;
      double precision, dimension(:), allocatable :: wl;
      double precision, dimension(:), allocatable :: sxl;
      double precision, dimension(:), allocatable :: rb;

      contains
          subroutine alloc_gaussl()
          use dirqfampar
          implicit none;

          allocate(  xl(0:NGL) );
          allocate(  wl(0:NGL) );
          allocate( sxl(0:NGL) );
          allocate(  rb(0:NGL) );

          end

          subroutine dealloc_gaussl()
          implicit none;

          deallocate( xl  );
          deallocate( wl  );
          deallocate( sxl );
          deallocate( rb  );

          end

      end module gaussl

      module herpol
      implicit none;

      !common /herpol/ qh(0:NZX,0:NGH),qh1(0:NZX,0:NGH)
      double precision, dimension(:,:), allocatable :: qh;
      double precision, dimension(:,:), allocatable :: qh1;

      contains
          subroutine alloc_herpol()
          use dirqfampar
          implicit none;

          allocate(  qh(0:NZX,0:NGH) );
          allocate( qh1(0:NZX,0:NGH) );

          end

          subroutine dealloc_herpol()
          implicit none;

          deallocate( qh  );
          deallocate( qh1 );

          end

      end module herpol

      module lagpol
      implicit none;

      !common /lagpol/ ql(0:2*NRX,0:MLX,0:NGL),ql1(0:2*NRX,0:MLX,0:NGL)
      double precision, dimension(:,:,:), allocatable :: ql;
      double precision, dimension(:,:,:), allocatable :: ql1;

      contains
          subroutine alloc_lagpol()
          use dirqfampar
          implicit none;

          allocate(  ql(0:2*NRX,0:MLX,0:NGL) );
          allocate( ql1(0:2*NRX,0:MLX,0:NGL) );

          end

          subroutine dealloc_lagpol()
          implicit none;

          deallocate( ql  );
          deallocate( ql1 );

          end

      end module lagpol

      module coulmb
      implicit none;

      !common /coulmb/ cou(MG),drvp(MG)
      double precision, dimension(:), allocatable :: cou;
      double precision, dimension(:), allocatable :: drvp;

      contains
          subroutine alloc_coulmb()
          use dirqfampar
          implicit none;

          allocate(  cou(MG) );
          allocate( drvp(MG) );

          end

          subroutine dealloc_coulmb()
          implicit none;

          deallocate( cou  );
          deallocate( drvp );

          end

      end module coulmb

      module procou
      implicit none;

      !common /procou/ ggc(MG,MG)
      double precision, dimension(:,:), allocatable :: ggc;

      contains
          subroutine alloc_procou()
          use dirqfampar
          implicit none;

          allocate( ggc(MG,MG) );

          end

          subroutine dealloc_procou()
          implicit none;

          deallocate( ggc  );

          end

      end module procou

      module gaussh
      implicit none;

      !common /gaussh/ xh(0:NGH),wh(0:NGH),zb(0:NGH)
      double precision, dimension(:), allocatable :: xh;
      double precision, dimension(:), allocatable :: wh;
      double precision, dimension(:), allocatable :: zb;

      contains
          subroutine alloc_gaussh()
          use dirqfampar
          implicit none;

          allocate( xh(0:NGH) );
          allocate( wh(0:NGH) );
          allocate( zb(0:NGH) );

          end

          subroutine dealloc_gaussh()
          implicit none;

          deallocate( xh );
          deallocate( wh );
          deallocate( zb );

          end

      end module gaussh

      module constr
      implicit none;

      !common /constr/ vc(MG,2)
      double precision, dimension(:,:), allocatable :: vc;

      contains
          subroutine alloc_constr()
          use dirqfampar
          implicit none;

          allocate( vc(MG,2) );

          end

          subroutine dealloc_constr()
          implicit none;

          deallocate( vc );

          end

      end module constr

      module couplf
      implicit none;

      !common /couplf/ ff(MG,4,2)
      double precision, dimension(:,:,:), allocatable :: ff;

      contains
          subroutine alloc_couplf()
          use dirqfampar
          implicit none;

          allocate( ff(MG,4,2) );

          end

          subroutine dealloc_couplf()
          implicit none;

          deallocate( ff );

          end

      end module couplf

      module deltamod
      implicit none;

      !dimension pnn(NNNX)
      double precision, dimension(:), allocatable :: pnn;

      contains
          subroutine alloc_deltamod()
          use dirqfampar
          implicit none;

          allocate( pnn(NNNX) );

          end

          subroutine dealloc_deltamod()
          implicit none;

          deallocate( pnn );

          end

      end module deltamod

      module rokaos
      implicit none;

      !common /rokaos/ rosh(NHHX,NB2X),aka(MVX,2)
      double precision, dimension(:,:), allocatable :: rosh;
      double precision, dimension(:,:), allocatable :: aka;

      contains
          subroutine alloc_rokaos()
          use dirqfampar
          implicit none;

          allocate( rosh(NHHX,NB2X) );
          allocate( aka(MVX,2)      );

          end

          subroutine dealloc_rokaos()
          implicit none;

          deallocate( rosh );
          deallocate( aka );

          end

      end module rokaos

      module tmrwnn
      implicit none;

      !common /tmrwnn/ wnn(MVX,NNNX),nnmax
      integer nnmax;
      double precision, dimension(:,:), allocatable :: wnn;

      contains
          subroutine alloc_tmrwnn()
          use dirqfampar
          implicit none;

          allocate( wnn(MVX,NNNX) );

          end

          subroutine dealloc_tmrwnn()
          implicit none;

          deallocate( wnn );

          end

      end module tmrwnn

      module densitmod
      implicit none;

      !dimension drs(MG,2),drv(MG,2)
      double precision, dimension(:,:), allocatable :: drs;
      double precision, dimension(:,:), allocatable :: drv;

      contains
          subroutine alloc_densitmod()
          use dirqfampar
          implicit none;

          allocate( drs(MG,2) );
          allocate( drv(MG,2) );

          end

          subroutine dealloc_densitmod()
          implicit none;

          deallocate( drs );
          deallocate( drv );

          end

      end module densitmod

      module dens
      implicit none;

      !common /dens  / ro(MG,4),dro(MG,4)
      double precision, dimension(:,:), allocatable :: ro;
      double precision, dimension(:,:), allocatable :: dro;

      contains
          subroutine alloc_dens()
          use dirqfampar
          implicit none;

          allocate(  ro(MG,4) );
          allocate( dro(MG,4) );

          end

          subroutine dealloc_dens()
          implicit none;

          deallocate( ro  );
          deallocate( dro );

          end

      end module dens

      module rhorho
      implicit none;

      !common /rhorho/ rs(MG,2),rv(MG,2)
      double precision, dimension(:,:), allocatable :: rs;
      double precision, dimension(:,:), allocatable :: rv;

      contains
          subroutine alloc_rhorho()
          use dirqfampar
          implicit none;

          allocate( rs(MG,2) );
          allocate( rv(MG,2) );

          end

          subroutine dealloc_rhorho()
          implicit none;

          deallocate( rs );
          deallocate( rv );

          end

      end module rhorho

      module gaucor
      implicit none;

      !common /gaucor/ ww(MG)
      double precision, dimension(:), allocatable :: ww;

      contains
          subroutine alloc_gaucor()
          use dirqfampar
          implicit none;

          allocate( ww(MG) );

          end

          subroutine dealloc_gaucor()
          implicit none;

          deallocate( ww );

          end

      end module gaucor

      module ptenso
      implicit none;

      !common /ptenso/ aka(MVX,2),dkd(MVX,2)
      double precision, dimension(:,:), allocatable :: aka;
      double precision, dimension(:,:), allocatable :: dkd;

      contains
          subroutine alloc_ptenso()
          use dirqfampar
          implicit none;

          allocate( aka(MVX,2) );
          allocate( dkd(MVX,2) );

          end

          subroutine dealloc_ptenso()
          implicit none;

          deallocate( aka );
          deallocate( dkd );

          end

      end module ptenso

      module dirhbmod
      implicit none;

      !character*8 tbb(NHBX)
      !dimension hb(NHBQX),e(NHBX),ez(NHBX)
      double precision, dimension(:), allocatable :: hb;
      double precision, dimension(:), allocatable :: e;
      double precision, dimension(:), allocatable :: ez;
      character*8,      dimension(:), allocatable :: tbb;

      contains
          subroutine alloc_dirhbmod()
          use dirqfampar
          implicit none;

          allocate(  hb(NHBQX) );
          allocate(   e(NHBX)  );
          allocate(  ez(NHBX)  );
          allocate( tbb(NHBX)  );

          end

          subroutine dealloc_dirhbmod()
          implicit none;

          deallocate( hb  );
          deallocate( e   );
          deallocate( ez  );
          deallocate( tbb );

          end

      end module dirhbmod

      module gfvmod
      implicit none;

      !common /gfviv / iv(-IGFV:IGFV)
      !common /gfvsq / sq(0:IGFV)
      !common /gfvsqi/ sqi(0:IGFV)
      !common /gfvsqh/ sqh(0:IGFV)
      !common /gfvshi/ shi(0:IGFV)
      !common /gfvfak/ fak(0:IGFV)
      !common /gfvfad/ fad(0:IGFV)
      !common /gfvfi / fi(0:IGFV)
      !common /gfvfdi/ fdi(0:IGFV)
      !common /gfvwf / wf(0:IGFV)
      !common /gfvwfi/ wfi(0:IGFV)
      !common /gfvwfd/ wfd(0:IGFV)
      !common /gfvgm2/ gm2(0:IGFV)
      !common /gfvgmi/ gmi(0:IGFV)
      !common /gfvwg / wg(0:IGFV)
      !common /gfvwgi/ wgi(0:IGFV)
      !common /bin0/   bin(0:IGFV,0:IGFV)
      integer,          dimension(:),   allocatable :: iv;
      double precision, dimension(:),   allocatable :: sq;
      double precision, dimension(:),   allocatable :: sqi;
      double precision, dimension(:),   allocatable :: sqh;
      double precision, dimension(:),   allocatable :: shi;
      double precision, dimension(:),   allocatable :: fak;
      double precision, dimension(:),   allocatable :: fad;
      double precision, dimension(:),   allocatable :: fi;
      double precision, dimension(:),   allocatable :: fdi;
      double precision, dimension(:),   allocatable :: wf;
      double precision, dimension(:),   allocatable :: wfi;
      double precision, dimension(:),   allocatable :: wfd;
      double precision, dimension(:),   allocatable :: gm2;
      double precision, dimension(:),   allocatable :: gmi;
      double precision, dimension(:),   allocatable :: wg;
      double precision, dimension(:),   allocatable :: wgi;
      double precision, dimension(:,:), allocatable :: bin;

      contains
          subroutine alloc_gfvmod()
          use dirqfampar
          implicit none;

          allocate(  iv(   -IGFV:IGFV) );
          allocate(  sq(       0:IGFV) );
          allocate( sqi(       0:IGFV) );
          allocate( sqh(       0:IGFV) );
          allocate( shi(       0:IGFV) );
          allocate( fak(       0:IGFV) );
          allocate( fad(       0:IGFV) );
          allocate(  fi(       0:IGFV) );
          allocate( fdi(       0:IGFV) );
          allocate(  wf(       0:IGFV) );
          allocate( wfi(       0:IGFV) );
          allocate( wfd(       0:IGFV) );
          allocate( gm2(       0:IGFV) );
          allocate( gmi(       0:IGFV) );
          allocate(  wg(       0:IGFV) );
          allocate( wgi(       0:IGFV) );
          allocate( bin(0:IGFV,0:IGFV) );

          end

          subroutine dealloc_gfvmod()
          implicit none;

          deallocate( iv  );
          deallocate( sq  );
          deallocate( sqi );
          deallocate( sqh );
          deallocate( shi );
          deallocate( fak );
          deallocate( fad );
          deallocate( fi  );
          deallocate( fdi );
          deallocate( wf  );
          deallocate( wfi );
          deallocate( wfd );
          deallocate( gm2 );
          deallocate( gmi );
          deallocate( wg  );
          deallocate( wgi );
          deallocate( bin );

          end

      end module gfvmod

      module gaushmod
      implicit none;

      !dimension x(2*NGH),w(2*NGH)
      double precision, dimension(:), allocatable :: x;
      double precision, dimension(:), allocatable :: w;

      contains
          subroutine alloc_gaushmod()
          use dirqfampar
          implicit none;

          allocate( x(2*NGH) );
          allocate( w(2*NGH) );

          end

          subroutine dealloc_gaushmod()
          implicit none;

          deallocate( x );
          deallocate( w );

          end

      end module gaushmod

      module ekinmod
      implicit none;

      !dimension h0(NHHX)
      double precision, dimension(:), allocatable :: h0;

      contains
          subroutine alloc_ekinmod()
          use dirqfampar
          implicit none;

          allocate( h0(NHHX) );

          end

          subroutine dealloc_ekinmod()
          implicit none;

          deallocate( h0 );

          end

      end module ekinmod

      module single
      implicit none;

      !common /single/ sp(NFGX,NBX)
      double precision, dimension(:,:), allocatable :: sp;

      contains
          subroutine alloc_single()
          use dirqfampar
          implicit none;

          allocate( sp(NFGX,NBX) );

          end

          subroutine dealloc_single()
          implicit none;

          deallocate( sp );

          end

      end module single

      module eparmod
      implicit none;

      !dimension ro(NHHX)
      double precision, dimension(:), allocatable :: ro;

      contains
          subroutine alloc_eparmod()
          use dirqfampar
          implicit none;

          allocate( ro(NHHX) );

          end

          subroutine dealloc_eparmod()
          implicit none;

          deallocate( ro );

          end

      end module eparmod

      module fields
      implicit none;

      !common /fields/ phi(MG,4)
      double precision, dimension(:,:), allocatable :: phi;

      contains
          subroutine alloc_fields()
          use dirqfampar
          implicit none;

          allocate( phi(MG,4) );

          end

          subroutine dealloc_fields()
          implicit none;

          deallocate( phi );

          end

      end module fields

      module fieldmod
      implicit none;

      !dimension gm(4),so(MG),ph(MG)
      double precision, dimension(:), allocatable :: gm;
      double precision, dimension(:), allocatable :: so;
      double precision, dimension(:), allocatable :: ph;

      contains
          subroutine alloc_fieldmod()
          use dirqfampar
          implicit none;

          allocate( gm(4)  );
          allocate( so(MG) );
          allocate( ph(MG) );

          end

          subroutine dealloc_fieldmod()
          implicit none;

          deallocate( gm );
          deallocate( so );
          deallocate( ph );

          end

      end module fieldmod

      module propag
      implicit none;

      !common /propag/ gg(NOX,NOX,4),psi(NOX,MG)
      double precision, dimension(:,:,:), allocatable :: gg;
      double precision, dimension(:,:),   allocatable :: psi;

      contains
          subroutine alloc_propag()
          use dirqfampar
          implicit none;

          allocate(  gg(NOX,NOX,4) );
          allocate( psi(NOX,MG)    );

          end

          subroutine dealloc_propag()
          implicit none;

          deallocate( gg  );
          deallocate( psi );

          end

      end module propag

      module gordonmod
      implicit none;

      !dimension rn(NOX),sn(NOX)
      double precision, dimension(:), allocatable :: rn;
      double precision, dimension(:), allocatable :: sn;

      contains
          subroutine alloc_gordonmod()
          use dirqfampar
          implicit none;

          allocate( rn(NOX) );
          allocate( sn(NOX) );

          end

          subroutine dealloc_gordonmod()
          implicit none;

          deallocate( rn  );
          deallocate( sn  );

          end

      end module gordonmod

      module greemesmod
      implicit none;

      !dimension dd(NOX,NOX),gi(NOX,NOX)
      double precision, dimension(:,:), allocatable :: dd;
      double precision, dimension(:,:), allocatable :: gi;

      contains
          subroutine alloc_greemesmod()
          use dirqfampar
          implicit none;

          allocate( dd(NOX,NOX) );
          allocate( gi(NOX,NOX) );

          end

          subroutine dealloc_greemesmod()
          implicit none;

          deallocate( dd );
          deallocate( gi );

          end

      end module greemesmod

      module potpot
      implicit none;

      !common /potpot/ vps(MG,2),vms(MG,2)
      double precision, dimension(:,:), allocatable :: vps;
      double precision, dimension(:,:), allocatable :: vms;

      contains
          subroutine alloc_potpot()
          use dirqfampar
          implicit none;

          allocate( vps(MG,2) );
          allocate( vms(MG,2) );

          end

          subroutine dealloc_potpot()
          implicit none;

          deallocate( vps );
          deallocate( vms );

          end

      end module potpot

      module bospol
      implicit none;

      !common /bospol/  qhb(0:N0BX,0:NGH),qlb(0:N0BX,0:NGL))
      double precision, dimension(:,:), allocatable :: qhb;
      double precision, dimension(:,:), allocatable :: qlb;

      contains
          subroutine alloc_bospol()
          use dirqfampar
          implicit none;

          allocate( qhb(0:N0BX,0:NGH) );
          allocate( qlb(0:N0BX,0:NGL) );

          end

          subroutine dealloc_bospol()
          implicit none;

          deallocate( qhb );
          deallocate( qlb );

          end

      end module bospol

      module plotmod
      implicit none;

      !dimension pn(NOX1),zp(0:NGL,0:NGH,3)
      double precision, dimension(:,:,:), allocatable :: zp;
      double precision, dimension(:),     allocatable :: pn;

      contains
          subroutine alloc_plotmod()
          use dirqfampar
          implicit none;

          allocate( zp(0:NGL,0:NGH,3) );
          allocate( pn(NOX1)          );

          end

          subroutine dealloc_plotmod()
          implicit none;

          deallocate( zp );
          deallocate( pn );

          end

      end module plotmod

      module singdmod
      implicit none;

      !dimension pnosc(0:N0FX),vnz(1:MVX,0:N0FX),vnr(1:MVX,0:N0FX)
      !dimension wn(MVX)
      double precision, dimension(:),   allocatable :: pnosc;
      double precision, dimension(:,:), allocatable :: vnz;
      double precision, dimension(:,:), allocatable :: vnr;
      double precision, dimension(:),   allocatable :: wn;

      contains
          subroutine alloc_singdmod()
          use dirqfampar
          implicit none;

          allocate(     pnosc(0:N0FX) );
          allocate( vnz(1:MVX,0:N0FX) );
          allocate( vnr(1:MVX,0:N0FX) );
          allocate(    wn(MVX)        );

          end

          subroutine dealloc_singdmod()
          implicit none;

          deallocate( pnosc );
          deallocate( vnz   );
          deallocate( vnr   );
          deallocate( wn    );

          end

      end module singdmod

      module vvnpmod
      implicit none;

      !dimension pnosc(0:N0FX)
      double precision, dimension(:), allocatable :: pnosc;

      contains
          subroutine alloc_vvnpmod()
          use dirqfampar
          implicit none;

          allocate( pnosc(0:N0FX) );

          end

          subroutine dealloc_vvnpmod()
          implicit none;

          deallocate( pnosc );

          end

      end module vvnpmod

      module single2
      implicit none;

      !common /single2/ sprr(NDDX,2*NBX)
      double precision, dimension(:,:), allocatable :: sprr;

      contains
          subroutine alloc_single2()
          use dirqfampar
          implicit none;

          allocate( sprr(NDDX,2*NBX) );

          end

          subroutine dealloc_single2()
          implicit none;

          deallocate( sprr );

          end

      end module single2
