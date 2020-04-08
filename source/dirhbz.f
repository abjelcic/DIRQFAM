c----------------------------------------------------------------------c
c     DIRHBZ code downloaded directly from CPC program library.        c
c                                                                      c
c     File 'dirhb.par' is casted to 'dirqfampar.mod' Fortran module.   c
c                                                                      c
c     Input file 'dirhb.dat' is renamed to 'dirqfam.dat'.              c
c                                                                      c
c     Some unused subroutines in DIRHBZ code are deleted.              c
c                                                                      c
c     Minor changes are labeled with 'ab' mark.                        c
c                                                                      c
c     All static commons and local arrays with sizes previously        c
c     statically determined by 'dirhb.par' are recasted to             c
c     appropriate Fortran modules which are dynamically allocated      c
c     in the beginning of program execution.                           c
c                                                                      c
c     Plotting subroutine in DIRHBZ code is not called since QFAM      c
c     submodule will print the density anyway.                         c
c----------------------------------------------------------------------c




c======================================================================c
c
c     PROGRAM DIRHB-axial
c
c======================================================================c
c     Relativistic Hartree-Bogoliubov theory in a axially symmetric basis
c     Main part of the code
c----------------------------------------------------------------------c
c
c-------------------------------------
      program DIRHBZ

      implicit real*8 (a-h,o-z)

      !ab
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      !ab


      !ab
      call calcdim();
      call allocdirhbz();
      !ab
c
c---- sets data
      call default(.false.)
c
c---- reads in data
      call reader(.true.)
c
c---- force-parameters
      call forces(.true.)
c
c---- Gauss-Hermite mesh points
      call gaush( 2.D0 ,.false.)
      call gausl(.false.)
c
c---- oscillator basis for single particle states
      call base(.false.)
c
c---- preparations
      call prep(.false.)
c
c---- wavefunctions at Gauss-Meshpoints
      call gaupol(.false.)
c
c---- initialization of the potentials
      call inout(1,.false.)
c
c---- initialization of the pairing field
      call dinout(1,.false.)
      call start(.false.)
c
c---- single-particle matix elements
      call singf(.false.)
c
c---- preparation of pairing matrix elements
      call singd(.false.)
c
c---- coulomb and meson propagators
      call greecou(.false.)
      call greemes(.false.)
c
c---- iteration
      call iter(.true.)
c
c---- transformation to the canonical basis
      call canon(.true.)
c
c---- center of mass correction
      call centmas(.false.)
c
c---- results
      call resu(.true.)
c

!ab   call plot(.false.)

c
c---- punching of potentials to tape  dis.wel
      call inout(2,.false.)
      call dinout(2,.false.)

      !ab
      close(l6);
      !ab (close dirhb.out file)


      !ab
c-----FAM submodule
      call main_fam( .false. );
      !ab

      !ab
      call deallocdirhbz();
      !ab

!ab   stop ' FINAL STOP OF DIRHBZ'
c-end-DIZ
      end





      !ab
c======================================================================c

      subroutine calcdim()

c======================================================================c

      use dirqfampar;
      implicit none;

      integer tape;
      integer nn0f, nn0b;
      integer NNGH, NNGL;
      integer JJ_multipole, KK_multipole;
      integer nf, ng;
      integer i, j;
      integer nz, nr, ml;
      character, dimension(:), allocatable :: fg_spx;
      integer,   dimension(:), allocatable :: nz_spx;
      integer,   dimension(:), allocatable :: nr_spx;
      integer,   dimension(:), allocatable :: ml_spx;
      integer NNBSX;
      integer broyd_m, broyd_n;
      character fg1, fg2;
      integer   nz1, nz2, nr1, nr2, ml1, ml2;
      integer NNWMAX;
      integer NRR, NZZ;
      integer NNCOORD, KKTRUNC;
      integer NNMESMAX, JJ_MAX, NNHMAX;
      integer NNTX, KKX, MMVX, MMVTX;
      integer dimf, dimg, suma, lam, N;






c-----Reading parameters which determine dimensions of arrays/matrices
      tape = 100;
      open( tape , file = 'dirqfam.dat' , status = 'old' );
      read(tape,'(10x,2i5)') nn0f, nn0b;
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
      read(tape,'(18x,i9)') NNGH;
      read(tape,'(18x,i9)') NNGL;
      read(tape,*);
      read(tape,'(18x,i9)') JJ_multipole;
      read(tape,'(18x,i9)') KK_multipole;
      close(tape);

      call assert( MOD(nn0f,2).eq.0  , 'n0f has to be even number'   );
      call assert( MOD(nn0b,2).eq.0  , 'n0b has to be even number'   );
      call assert( NNGH.ge.30        , 'NGH has to be at least 30'   );
      call assert( NNGL.ge.30        , 'NGL has to be at least 30'   );
      call assert( NNGH.lt.100       , 'NGH has to be less than 100' );
      call assert( NNGL.lt.100       , 'NGL has to be less than 100' );






c-----Construction of quantum numbers
      nf = ((nn0f+1)*(nn0f+2)*(nn0f+3))/6;
      ng = ((nn0f+2)*(nn0f+3)*(nn0f+4))/6;
      allocate( fg_spx(nf+ng) );
      allocate( nz_spx(nf+ng) );
      allocate( nr_spx(nf+ng) );
      allocate( ml_spx(nf+ng) );
      i = 0;
      do nz = 0 , nn0f
          do nr = 0 , nn0f/2
              do ml = -nn0f , nn0f
                  if( nz+2*nr+abs(ml) .le. nn0f ) then
                      i = i + 1;
                      fg_spx(i) = 'f';
                      nz_spx(i) = nz;
                      nr_spx(i) = nr;
                      ml_spx(i) = ml;
                  endif
              enddo
          enddo
      enddo
      do nz = 0 , nn0f+1
          do nr = 0 , (nn0f+1)/2
              do ml = -(nn0f+1) , nn0f+1
                  if( nz+2*nr+abs(ml) .le. nn0f+1 ) then
                      i = i + 1;
                      fg_spx(i) = 'g';
                      nz_spx(i) = nz;
                      nr_spx(i) = nr;
                      ml_spx(i) = ml;
                  endif
              enddo
          enddo
      enddo
      call assert( nf+ng.eq.i , 'nf+ng =/= i' );






c-----Calculation of maximal block dimension (NBSX) of the U,V matrices
c-----(in fact, this is the dimension of the omega^pi = 1/2^+ block)
      NNBSX = ((nn0f+2)*(nn0f+3))/2;






c-----Calculating Broyden's vector length
      broyd_m = 35;
      broyd_n = 0;
      ! dh_1 matrix
      do i = 1 , nf+ng
          fg1 = fg_spx(i);
          ml1 = ml_spx(i);
          do j = i , nf+ng
              fg2 = fg_spx(j);
              ml2 = ml_spx(j);
              if( fg1 .eq. fg2 ) then
                  if( abs(ml1-ml2) .eq. KK_multipole ) then
                      broyd_n = broyd_n + 1;
                  endif
              endif
              if( fg1 .ne. fg2 ) then
                  if( abs(ml1-ml2) .eq. KK_multipole ) then
                      broyd_n = broyd_n + 1;
                  endif
                  if( abs(ml1+ml2+1) .eq. KK_multipole ) then
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
              if( abs(ml1-ml2) .ne. KK_multipole ) CYCLE;
              broyd_n = broyd_n + 1; !dDelta_pl
              broyd_n = broyd_n + 1; !dDelta_mi
          enddo
      enddo
      broyd_n = 2*broyd_n; !protons/neutrons
      broyd_n = 2*broyd_n; !real/imaginary






c-----Calculating NWMAX
      NNWMAX = 0;
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
              if( abs(ml1-ml2) .ne. KK_multipole ) CYCLE;

              do NRR = 0 , nr1+nr2+(abs(ml1)+abs(ml2)-KK_multipole)/2
                  do NZZ = 0 , nz1+nz2
                      if( mod(NZZ,2) .ne. mod(nz1+nz2,2) ) CYCLE;

                      NNWMAX = NNWMAX + 1;

                  enddo
              enddo

          enddo
      enddo






c-----Setting KTRUNC and NCOORD
      NNCOORD = (2*NNGH)*NNGL;
      KKTRUNC = NNBSX + 10;






c-----Setting NMESMAX shells and number of nz+2*nr+K <= NMESMAX pairs
      NNMESMAX = nn0b;
      JJ_MAX = 3;
      call assert( NNMESMAX-(JJ_MAX+1).gt.0 , 'n0b too small' );
      NNHMAX = ( (NNMESMAX+1)*(NNMESMAX+3) + 1 - MOD(NNMESMAX,2) )/4;






c-----Calculating NTX, KX, MVX, MVTX
      NNTX = (nn0f+2)*(nn0f+3)*(2*nn0f+5)/6;
      KKX  = (nn0f+2)*(nn0f+3)*(nn0f+4)/6;
      ! f-component
      dimf = 0;
      ! Positive parity
      do i = 1 , nn0f+1
          suma = 0;
          do lam = i-1 , i
              do N = lam+mod(lam,2) , nn0f , 2
                 suma = suma + (N-lam)/2 + 1;
              enddo
          enddo
          dimf = dimf + (suma*(suma+1))/2;
      enddo
      ! Negative parity
      do i = 1 , nn0f+1
          suma = 0;
          do lam = i-1 , i
              do N = lam+mod(lam+1,2) , nn0f , 2
                  suma = suma + (N-lam)/2 + 1;
              enddo
          enddo
          dimf = dimf + (suma*(suma+1))/2;
      enddo
      MMVX = dimf;
      ! g-component
      dimg = 0;
      ! Positive parity
      do i = 1 , (nn0f+1) + 1
          suma = 0;
          do lam = i-1 , i
              do N = lam+mod(lam,2) , (nn0f+1) , 2
                  suma = suma + (N-lam)/2 + 1;
              enddo
          enddo
          dimg = dimg + (suma*(suma+1))/2;
      enddo
      ! Negative parity
      do i = 1 , (nn0f+1)+1
          suma = 0;
          do lam = i-1 , i
              do N = lam+mod(lam+1,2) , (nn0f+1) , 2
                  suma = suma + (N-lam)/2 + 1;
              enddo
          enddo
          dimg = dimg + (suma*(suma+1))/2;
      enddo
      MMVTX = dimf + dimg;






c-----Setting dirqfampar module
      IGFV       = 100;

      RMAX       = 15.D0;;
      MR         = 301;

      QMAX       = 6.D0;
      MQ         = 49;

      NGH        = NNGH;
      NGL        = NNGL;
      NGLEG      = 18;

      N0FX       = nn0f;
      nxx        = N0FX/2;

      NBX        = 2*N0FX+3;

      NTX        = NNTX;

      KX         = KKX;

      NZX        =  N0FX+1;

      NRX        = N0FX/2;

      MLX        = N0FX+1;

      NFX        = (N0FX/2+1)**2;

      NGX        = (N0FX/2+1)*(N0FX/2+2);

      NDX        = NGX;

      N0BX       = nn0b;

      nbxx       = N0BX/2;
      NOX        = ((nbxx+1)*(nbxx+2))/2;

      NOX1       = ((N0FX+1)*(N0FX+2))/2;

      NGH2       = NGH+NGH;
      NB2X       = NBX+NBX;
      NHX        = NFX+NGX;
      NDDX       = NDX*NDX;
      NHBX       = NHX+NHX;
      NHBQX      = NHBX*NHBX;
      NFFX       = NFX*NFX;
      NFGX       = NFX*NGX;
      NHHX       = NHX*NHX;
      MG         = (NGH+1)*(NGL+1);
      N02        = 2*N0FX;
      NNNX       = (N0FX+1)*(N0FX+1);

      MVX        = MMVX;
      MVTX       = MMVTX;

      NBSX       = NNBSX;
      J_MAX      = JJ_MAX;
      NFAM_BROYD = broyd_n;
      MFAM_BROYD = broyd_m;
      NWMAX      = NNWMAX;
      KTRUNC     = KKTRUNC;
      NCOORD     = NNCOORD;
      NMESMAX    = NNMESMAX;
      NHMAX      = NNHMAX;






      deallocate( fg_spx );
      deallocate( nz_spx );
      deallocate( nr_spx );
      deallocate( ml_spx );

      end

c======================================================================c

      subroutine allocdirhbz()

c======================================================================c

      use blokap;
      use bloosc;
      use quaosc;
      use bosqua;
      use vvvikf;
      use deldel;
      use gamgam;
      use broyde;
      use canonmod;
      use blodir;
      use eeecan;
      use blocan;
      use waveuv;
      use centmasmod;
      use gaussl;
      use herpol;
      use lagpol;
      use coulmb;
      use procou;
      use gaussh;
      use constr;
      use couplf;
      use deltamod;
      use rokaos;
      use tmrwnn;
      use densitmod;
      use dens;
      use rhorho;
      use gaucor;
      use ptenso;
      use dirhbmod;
      use gfvmod;
      use gaushmod;
      use ekinmod;
      use single;
      use eparmod;
      use fields;
      use fieldmod;
      use propag;
      use gordonmod;
      use greemesmod;
      use potpot;
      use bospol;
      use plotmod;
      use singdmod;
      use vvnpmod;
      use single2;

      implicit none;

      call alloc_blokap();
      call alloc_bloosc();
      call alloc_quaosc();
      call alloc_bosqua();
      call alloc_vvvikf();
      call alloc_deldel();
      call alloc_gamgam();
      call alloc_broyde();
      call alloc_canonmod();
      call alloc_blodir();
      call alloc_eeecan();
      call alloc_blocan();
      call alloc_waveuv();
      call alloc_centmasmod();
      call alloc_gaussl();
      call alloc_herpol();
      call alloc_lagpol();
      call alloc_coulmb();
      call alloc_procou();
      call alloc_gaussh();
      call alloc_constr();
      call alloc_couplf();
      call alloc_deltamod();
      call alloc_rokaos();
      call alloc_tmrwnn();
      call alloc_densitmod();
      call alloc_dens();
      call alloc_rhorho();
      call alloc_gaucor();
      call alloc_ptenso();
      call alloc_dirhbmod();
      call alloc_gfvmod();
      call alloc_gaushmod();
      call alloc_ekinmod();
      call alloc_single();
      call alloc_eparmod();
      call alloc_fields();
      call alloc_fieldmod();
      call alloc_propag();
      call alloc_gordonmod();
      call alloc_greemesmod();
      call alloc_potpot();
      call alloc_bospol();
      call alloc_plotmod();
      call alloc_singdmod();
      call alloc_vvnpmod();
      call alloc_single2();

      end

c======================================================================c

      subroutine deallocdirhbz()

c======================================================================c

      use blokap;
      use bloosc;
      use quaosc;
      use bosqua;
      use vvvikf;
      use deldel;
      use gamgam;
      use broyde;
      use canonmod;
      use blodir;
      use eeecan;
      use blocan;
      use waveuv;
      use centmasmod;
      use gaussl;
      use herpol;
      use lagpol;
      use coulmb;
      use procou;
      use gaussh;
      use constr;
      use couplf;
      use deltamod;
      use rokaos;
      use tmrwnn;
      use densitmod;
      use dens;
      use rhorho;
      use gaucor;
      use ptenso;
      use dirhbmod;
      use gfvmod;
      use gaushmod;
      use ekinmod;
      use single;
      use eparmod;
      use fields;
      use fieldmod;
      use propag;
      use gordonmod;
      use greemesmod;
      use potpot;
      use bospol;
      use plotmod;
      use singdmod;
      use vvnpmod;
      use single2;

      implicit none;

      call dealloc_blokap();
      call dealloc_bloosc();
      call dealloc_quaosc();
      call dealloc_bosqua();
      call dealloc_vvvikf();
      call dealloc_deldel();
      call dealloc_gamgam();
      call dealloc_broyde();
      call dealloc_canonmod();
      call dealloc_blodir();
      call dealloc_eeecan();
      call dealloc_blocan();
      call dealloc_waveuv();
      call dealloc_centmasmod();
      call dealloc_gaussl();
      call dealloc_herpol();
      call dealloc_lagpol();
      call dealloc_coulmb();
      call dealloc_procou();
      call dealloc_gaussh();
      call dealloc_constr();
      call dealloc_couplf();
      call dealloc_deltamod();
      call dealloc_rokaos();
      call dealloc_tmrwnn();
      call dealloc_densitmod();
      call dealloc_dens();
      call dealloc_rhorho();
      call dealloc_gaucor();
      call dealloc_ptenso();
      call dealloc_dirhbmod();
      call dealloc_gfvmod();
      call dealloc_gaushmod();
      call dealloc_ekinmod();
      call dealloc_single();
      call dealloc_eparmod();
      call dealloc_fields();
      call dealloc_fieldmod();
      call dealloc_propag();
      call dealloc_gordonmod();
      call dealloc_greemesmod();
      call dealloc_potpot();
      call dealloc_bospol();
      call dealloc_plotmod();
      call dealloc_singdmod();
      call dealloc_vvnpmod();
      call dealloc_single2();

      end
      !ab






c======================================================================c

      subroutine base(lpr)

c======================================================================c
c
c     determines the basis for Dirac equation in zylindric coordinates
c
c     NB        number of K-parity-blocks
c     abs(KB)   K-quantum number of each block (K+1/2)
c     sign(KB)  parity of each block ( +1 or  -)
c     IA(ib,1): begin of the large components of block b is IA(b,1)+1
c     IA(ib,2): begin of the small components of block b is IA(b,2)+1
c     ID(ib,1): dimension large components of block b
c     ID(ib,2): dimension small components of block b
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use blokap;
      use bloosc;
      use quaosc;
      use bosqua;
      use vvvikf;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
      character tp*1,tis*1,tit*8,tl*1                           ! textex
!ab   character nucnam*2 (not even used...)                     ! nucnuc
!ab   character tb*6                                            ! blokap
!ab   character tt*8                                            ! quaosc
c
c
      common /basnnn/ n0f,n0b
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
!ab   common /quaosc/ nt,nz(NTX),nr(NTX),ml(NTX),ms(NTX),np(NTX),tt(NTX)
!ab   common /bosqua/ nzb(NOX),nrb(NOX),NO
      common /zdimos/ nzm,nrm,mlm
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
!ab   common /vvvikf/ mv,ipos(NBX),nib(MVX),nni(2,MVX)
c
      if (lpr) then
      write(l6,*) '****** BEGIN BASE **********************************'
      endif
c
c
c
c======================================================================c
c     Oscillator-Base for Fermions:
c======================================================================c
c
c
c----------------------------------------------------------------------c
c     Construction of the different K - Parity blocks.
c----------------------------------------------------------------------c
      nzm = 0
      nrm = 0
      mlm = 0
      nfx0 = 0
      ngx0 = 0
c
      ib = 0
      il = 0
      ik = 0
c     loop over K-quantum number
      do k  = 1,n0f+1
c
c     loop over the parity
      do 10 ipk = 1,2
	     if (k.eq.n0f+1.and.mod(k+ipk,2).ne.0) goto 10
c
         if (ib+1.gt.NBX) stop ' in BASE: NBX too small'
	     ib = ib + 1
	     ipar   = 2*mod(ipk,2) - 1
         kb(ib) = ipar*k
	     mb(ib) = 2
         write(tb(ib),'(i2,2h/2,a1)') k+k-1,tp(ipk)
         ipf = ipk
         ipg = 3 - ipk

c        loop over large and small components
         do ifg = 1,2
	        ia(ib,ifg) = il
	        !ab
	        ip = 0;
            !ab
            if (ifg.eq.1) ip = ipf
            if (ifg.eq.2) ip = ipg
c
c           loop over major quantum number nn
            if (mod(k+ip,2).eq.0) then
               n0 = k-1
            else
               n0 = k
            endif
            do nn = n0,n0f+ifg-1,2
c
c           loop over quantum number nz
            do iz = 0,nn
c
c           loop over quantum number ms
            do is = 0,1
               im = k-is
               ir = nn - iz - im
               if (ir.ge.0.and.mod(ir,2).eq.0) then
                  ir = ir/2
c
                  il = il + 1
	          if (il.gt.NTX) stop 'in BASE: NTX too small'
                  nz(il) = iz
                  nr(il) = ir
                  ml(il) = im
                  ms(il) = k-im
	              np(il) = ifg
                  write(tt(il),100) iz,ir,im
c                 write(6,100) iz,ir,im               ! remove
c                 read*                               ! remove
                  nzm = max0(nzm,nz(il))
                  nrm = max0(nrm,nr(il))
                  mlm = max0(mlm,ml(il))
               endif
            enddo   ! is
            enddo   ! iz
            enddo   ! nn
            id(ib,ifg) = il - ia(ib,ifg)
	     enddo   ! ifg
         nf = id(ib,1)
         ng = id(ib,2)
         nh = nf + ng
	     nfx0 = max(nfx0,id(ib,1))
	     ngx0 = max(ngx0,id(ib,2))
	     ik   = ik + max0(id(ib,1),id(ib,2))
   10 continue   ! ipf
      enddo   ! k




      !ab
c=======================================================================c
c=======================================================================c
c=======================================================================c

      il = il + 1;

      ib     = ib + 1;
      kb(ib) = -( n0f + 1 );
      mb(ib) = 2;
      write(tb(ib),'(i2,2h/2,a1)') 2*abs(kb(ib))-1 , '-';

      id(ib,1) = 0;
      id(ib,2) = 1;
      ia(ib,1) = il - 1;
      ia(ib,2) = il - 1;
      nz(il)   = 0;
      nr(il)   = 0;
      ml(il)   = n0f + mod(n0f,2);
      ms(il)   = 2*mod(n0f,2) - 1;
      np(il)   = 2;
      write(tt(il),100) nz(il) , nr(il) , ml(il);


      il = il + 1;

      ib     = ib + 1;
      kb(ib) = -( 2*mod(n0f,2) - 1 )*( n0f + 2 );
      mb(ib) = 2;
      if ( kb(ib) .gt. 0 ) then
          write(tb(ib),'(i2,2h/2,a1)') 2*abs(kb(ib))-1 , '+';
      else
          write(tb(ib),'(i2,2h/2,a1)') 2*abs(kb(ib))-1 , '-';
      endif


      id(ib,1) = 0;
      id(ib,2) = 1;
      ia(ib,1) = il - 1;
      ia(ib,2) = il - 1;
      nz(il)   = 0;
      nr(il)   = 0;
      ml(il)   = n0f + 1;
      ms(il)   = 1
      np(il)   = 2;
      write(tt(il),100) nz(il) , nr(il) , ml(il);

c=======================================================================c
c=======================================================================c
c=======================================================================c
      !ab



      nb = ib
      nt = il
      nk = ik
c
c----------------------------------------------------------------------c
c     Construction of two-body pairs (i1,i2)
c----------------------------------------------------------------------c
c---- only f-pairs
      il = 0
      do ib = 1,nb
         ipos(ib) = il
         nf = id(ib,1)
         do n2 =  1,nf
         do n1 = n2,nf
            il = il + 1
            nib(il)   = ib
            nni(1,il) = n1
            nni(2,il) = n2
         enddo  ! n1
         enddo  ! n2
      enddo  ! ib
      mv = il
      if (lpr) write(l6,105) ' number of pairs f  mv = ',mv,MVX
      if (mv.gt.MVX) then
	 write(6,*) 'mv =',mv,' MVX = ',MVX
	 stop ' in BASE: MVX too small'
      endif
c
c
      if (lpr) then
         write(l6,*) ' '
         write(l6,*)   'Maximal values:             needed     given'
         write(l6,103) ' Number of blocks: nb  = ',nb,NBX
         write(l6,103) ' Number of levels  nt  = ',nt,NTX
         write(l6,103) ' Number of levels  nk  = ',nk,KX
         write(l6,103) ' Maximal nz:       nzm = ',nzm,NZX
         write(l6,103) ' Maximal nr:       nrm = ',nrm,NRX
         write(l6,103) ' Maximal ml:       mlm = ',mlm,MLX
         write(l6,103) ' Maximal nf        nfm = ',nfx0,NFX
         write(l6,103) ' Maximal ng        ngm = ',ngx0,NGX
         write(l6,103) ' Maximal mv        mv  = ',mv,MVX
         write(l6,103) ' Maximal mv1       mv1 = ',mv1,MVTX
      endif
c
c----------------------------------------------------------------------c
c     Printout
c----------------------------------------------------------------------c
      if (lpr) then
         do ib = 1,nb
	        kap = kb(ib)
	        nf  = id(ib,1)
	        ng  = id(ib,2)
	        nh  = nf + ng
            i0f = ia(ib,1)
            i0g = ia(ib,2)
            write(l6,'(/,a,6i4)') tb(ib),nf,ng,i0f,i0g
            do i = i0f+1,i0f+nh
               nn = nz(i) + 2*nr(i) + ml(i)
	           mx = 2*(iabs(kap) - ml(i)) - 1
               write(l6,102) i,'   NN = ',nn,'   nz = ',nz(i),
     &         '   nr = ',nr(i),'   ml =',ml(i),'   ms = ',mx,tt(i)
	       if (i.eq.i0f+nf) write(l6,'(3x,61(1h-))')
      read*
            enddo  ! i
         enddo   ! ib
      endif
c
c
c
c======================================================================c
c     Oscillator Base for Bosons:
c======================================================================c
      if (mod(N0B,2).ne.0) stop ' in BASE: N0B must be even'
      il   = 0
c
c---- loop over nr-quantum number
      do ir = 0,N0B/2
c
c---- loop over nz-quantum number
      do iz = 0,N0B,2
         if (iz + 2*ir.le.N0B) then
            il = il + 1
            nrb(il) = ir
            nzb(il) = iz
            call icheck_n1len2(il,NOX,'in BASE: il too large')
         endif
      enddo   ! iz
      enddo   ! ir
      NO=il
      if (lpr) then
         write(l6,'(//,a,/)') ' Boson base: '
         do i = 1,NO
            nn = nzb(i) + 2*nrb(i)
            write(l6,104) i,'  nn = ',nn,'  nz = ',nzb(i),
     &                                   '  nr = ',nrb(i)
         enddo
         write(l6,*) ' '
         write(l6,*)   'Maximal values:             needed     given'
         write(l6,103) ' Number of levels NO   = ',NO,NOX
      endif
c
      if (lpr) then
      write(l6,*) '****** END BASE ************************************'
      endif
c
  100 format('[',3i2,']')
!ab 101 format(4i2)
  102 format(i4,a,i2,a,i2,a,i2,a,i2,a,i2,3h/2 ,a)
  103 format(a,2i10)
  104 format(i5,a,i3,a,i3,a,i3)
  105 format(a,2i6)
c
      return
c-end-BASE
      end

c======================================================================c

      subroutine broyden(lpr)

c======================================================================c
      use dirqfampar;
      use bloosc;
      use blokap;
      use deldel;
      use gamgam;
      use broyde;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
!ab   character tb*6                                           ! blokap
c
      ! broyden iteration sizes
!ab   parameter (nn = 2*MVTX+2*MVX+2+1)
!ab   parameter (mm = 7)
c
      common /baspar/ hom,hb0,b0
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
!ab   common /bloosc/ ia(nbx,2),id(nbx,2)
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
      common /physco/ hbc,alphi,r0
!ab   common /deldel/ de(nhhx,nb2x)
!ab   common /gamgam/ hh(nhhx,nb2x)
      common /fermi / ala(2),tz(2)
      common /con_b2/ betac,q0c,cquad,c0,alaq,calcq0,icstr
      common /pair  / del(2),spk(2),spk0(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
!ab   common /broyde/ bbeta(mm,mm),df(nn,mm),dv(nn,mm),
!ab  &                bwork(mm),curv(nn),bw0,ibwork(mm)
!ab   common /broyde1/ vin(nn)
!ab   common /broyde2/ ibroyd
!ab   dimension vou(nn)
      data bmix /0.7d0/
c
      if (lpr) then
      write(l6,*) '****** BEGIN BROYDEN *************************'
      endif

      xmi = xmix
      si  = zero
      do i=1,nn
         vou(i)=zero
      enddo
c
      if (ii.eq.0) then
         do i=1,nn
            vin(i)=zero
         enddo
         ipos=0
         do it = 1,2
            do ib=1,nb
               nf=id(ib,1)
               ng=id(ib,2)
               nh=nf+ng
               m = ib+(it-1)*nbx
               do i2 = 1,nf
                  do i1 = i2,nf
                     ipos=ipos+1
                     vin(ipos) = hh(i1+(i2-1)*nh,m)
                  enddo
               enddo
               do i2 = nf+1,nh
                  do i1 = i2,nh
                     ipos=ipos+1
                     vin(ipos) = hh(i1+(i2-1)*nh,m)
                  enddo
               enddo
               do i2 = 1,nf
                  do i1 = i2,nf
                     ipos=ipos+1
                     vin(ipos) = de(i1+(i2-1)*nh,m)
                  enddo
               enddo
            enddo  !ib
            ipos=ipos+1
            vin(ipos)=ala(it)
         enddo  !it
         ipos=ipos+1
         vin(ipos)=alaq
         return
      endif
      ipos=0
      do it = 1,2
         do ib=1,nb
            nf=id(ib,1)
            ng=id(ib,2)
            nh=nf+ng
            m = ib+(it-1)*nbx
            do i2 = 1,nf
               do i1 = i2,nf
                  ipos=ipos+1
                  vou(ipos) = hh(i1+(i2-1)*nh,m) -vin(ipos)
               enddo
            enddo
            do i2 = nf+1,nh
               do i1 = i2,nh
                  ipos=ipos+1
                  vou(ipos) = hh(i1+(i2-1)*nh,m) -vin(ipos)
               enddo
            enddo
            do i2 = 1,nf
               do i1 = i2,nf
                  ipos=ipos+1
                  vou(ipos) = de(i1+(i2-1)*nh,m) -vin(ipos)
                  if(del(it).lt.1.d-5) then
                     vou(ipos)=zero
                     vin(ipos)=zero
                  endif
               enddo
            enddo
         enddo
         ipos=ipos+1
         vou(ipos)=ala(it)-vin(ipos)
      enddo
      ipos=ipos+1
      vou(ipos)=alaq-vin(ipos)
      nnn = ipos
      ! calculate si
      do i=1,nnn
         si = max( si, abs(vou(i)) )
      enddo
      ! broyden's mixing procedure starts here...
      if ( mm .eq. 0 .or. ii.eq. 1 .or. ibroyd.eq.0) then ! linear mixing
         do i = 1, nnn
            vin(i) = vin(i) + xmi*vou(i)
         enddo
         ilast = 0
      else                                                ! broyden mixing

         iuse = min( ii-1 - 1, mm )
         ipos = ii - 2 - ( (ii-3)/mm )*mm
         inex = ii - 1 - ( (ii-2)/mm )*mm

         if( ii .eq. 2 ) then
            do j = 1, mm
               do i = 1, nnn
                  df(i,j) = zero
                  dv(i,j) = zero
               enddo
            enddo
            bw0 = 0.01d0
         else
            do i = 1, nnn
               df(i,ipos) = vou(i) - df(i,ipos)
               dv(i,ipos) = vin(i) - dv(i,ipos)
            enddo

            dnorm = sqrt( dnrm2(nnn,df(1,ipos),1)**2.0d0 )
            call dscal( nnn, 1.0d0 / dnorm, df(1,ipos), 1 )
            call dscal( nnn, 1.0d0 / dnorm, dv(1,ipos), 1 )
         endif

         do i = 1, iuse
            do j = i+1, iuse
               bbeta(i,j) = ddot( nnn, df(1,j), 1, df(1,i), 1 )
            enddo
            bbeta(i,i) = 1.0d0 + bw0*bw0
         enddo

         call dsytrf( 'U', iuse, bbeta, mm, ibwork, bwork, mm, info )
         if( info .ne. 0 ) stop '   in broyden: info at DSYTRF V+S '

         call dsytri( 'U', iuse, bbeta, mm, ibwork, bwork, info )
         if( info .ne. 0 ) stop '   in broyden: info at DSYTRI V+S '

         do i = 1, iuse
            do j = i+1, iuse
               bbeta(j,i) = bbeta(i,j)
            enddo
            bwork(i) = ddot( nnn, df(1,i), 1, vou, 1 )
         enddo

         do i = 1, nnn
            curv(i) = bmix * vou(i)
         enddo
         do i = 1, iuse
            gamma = 0.0d0
            do j = 1, iuse
               gamma = gamma + bbeta(j,i) * bwork(j)
            enddo
            do k = 1, nnn
               curv(k) = curv(k) - gamma * ( dv(k,i) + bmix * df(k,i) )
            enddo
         enddo

         call dcopy( nnn, vou, 1, df(1,inex), 1 )
         call dcopy( nnn, vin, 1, dv(1,inex), 1 )

         curvature = ddot( nnn, vou, 1, curv, 1 )
         if( curvature .gt. 0.0d0 ) then
            ilast = 1
            do i = 1, nnn
               vin(i) = vin(i) + curv(i)
            enddo
c            write(*,100) 'broyden mixing: mm =', iuse, 'c=',curvature
         else
            ilast = 0
            do i = 1, nnn
               vin(i) = vin(i) + xmi*vou(i)
            enddo
c            write(*,100) 'linear  mixing: mm =', iuse, 'c=',curvature
         endif
      endif
!ab 100  format(10x,a,i2,2x,a,f16.8)
      ! broyden's mixing procedure ends here

      ! set the new matrix elements
      ipos=0
      do it = 1,2
         do ib=1,nb
            nf=id(ib,1)
            ng=id(ib,2)
            nh=nf+ng
            m = ib+(it-1)*nbx
            do i2 = 1,nf
               do i1 = i2,nf
                  ipos=ipos+1
                  hh(i1+(i2-1)*nh,m) = vin(ipos)
                  hh(i2+(i1-1)*nh,m) = vin(ipos)
               enddo
            enddo
            do i2 = nf+1,nh
               do i1 = i2,nh
                  ipos=ipos+1
                  hh(i1+(i2-1)*nh,m) = vin(ipos)
                  hh(i2+(i1-1)*nh,m) = vin(ipos)
               enddo
            enddo
            do i2 = 1,nf
               do i1 = i2,nf
                  ipos=ipos+1
                  de(i1+(i2-1)*nh,m) = vin(ipos)
                  de(i2+(i1-1)*nh,m) = vin(ipos)
               enddo
            enddo
         enddo  !ib
         ipos=ipos+1
         ala(it) =vin(ipos)
      enddo  !it
      ipos=ipos+1
      alaq=vin(ipos)
c
c
      if (lpr) then
      write(l6,*) '****** END BROYDEN *************************'
      endif
c
      return
c-end-BROYDEN
      end

c======================================================================cz

      subroutine canon(lpr)

c======================================================================c
c
c     transforms to the canonical basis
c     version for RHB
c----------------------------------------------------------------------c
      use dirqfampar;
      use canonmod;
      use blodir;
      use blokap;
      use bloosc;
      use gamgam;
      use deldel;
      use eeecan;
      use blocan;
      use quaosc;
      use waveuv;

      implicit real*8 (a-h,o-z)
!ab   include 'dirqfam.par'
c
!ab   logical lpr,lpr1
      logical lpr;
c
!ab   character tb*6                                            ! blokap
!ab   character tt*8                                            ! quaosc
      character tp*1,tis*1,tit*8,tl*1                           ! textex
c
!ab   dimension aa(NHHX),dd(NHHX),v2(NHX),z(NHX),eb(NHX),h(NHX),d(NHX)
c
!ab   common /blodir/ ka(NBX,4),kd(NBX,4)
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
!ab   common /gamgam/ hh(NHHX,NB2X)
!ab   common /deldel/ de(NHHX,NB2X)
!ab   common /eeecan/ eecan(KX,4),decan(KX,4),vvcan(KX,4),
!ab  &                fgcan(NHX,KX,4),ibkcan(KX,4)
!ab   common /blocan/ kacan(nbx,4),kdcan(nbx,4),nkcan(4)
      common /mathco/ zero,one,two,half,third,pi
      common /fermi / ala(2),tz(2)
      common /optopt/ itx,icm,icou,ipc,inl,idd
!ab   common /quaosc/ nt,nz(NTX),nr(NTX),ml(NTX),ms(NTX),np(NTX),tt(NTX)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
!ab   common /waveuv/ fguv(NHBX,KX,4),equ(KX,4)

      data ash/100.d0/
c
      if (lpr) then
      write(l6,*) '****** BEGIN CANON *********************************'
      endif
c
c
c======================================================================c
      do it = 1,2   ! loop over neutrons and protons
c======================================================================c
         if (lpr) then
            write(l6,100) tit(it)
            write(l6,101) 'K pi','[nz,nr,ml]','smax',
     &                    'eecan','vvcan','decan'
         endif
         klp = 0
         kla = 0
c======================================================================c
      do ib = 1,nb   ! loop over the blocks
c======================================================================c
         nf  = id(ib,1)
         ng  = id(ib,2)
         nh  = nf + ng
         nhb = 2*nh
         i0f = ia(ib,1)
         i0g = ia(ib,2)
         m   = ib + (it-1)*NBX
         kf  = kd(ib,it)
         kg  = kd(ib,it+2)
         k0f = ka(ib,it)
         k0g = ka(ib,it+2)
c
c------- transformation to the canonical basis
c------- calculation of the generalized density V*VT
         do n2 = 1,nh
         do n1 = 1,nh
            s = zero
            do k = k0f+1,k0f+kf
               s = s + fguv(nh+n1,k,it)*fguv(nh+n2,k,it)
            enddo
            s1 = zero
            do k = k0g+1,k0g+kg
               s1 = s1 + fguv(nh+n1,k,it+2)*fguv(nh+n2,k,it+2)
            enddo
            aa(n1+(n2-1)*nh) = s + ash*s1
         enddo   ! n1
         enddo   ! n2
c
c------- diagonalizaton
         call sdiag(nh,nh,aa,v2,dd,z,1)
         eps=1.0e-6
         call degen(nh,nh,v2,dd,hh(1,m),eb,eps,aa,z)
c
c        major component of the wave function should be > 0
         do k = 1,nh
            cmax = zero
            do n = 1,nh
               s = abs(dd(n+(k-1)*nh))
               if (s.gt.abs(cmax)) cmax = dd(n+(k-1)*nh)
            enddo   ! n
            if (cmax.lt.zero) then
               do n = 1,nh
                  dd(n+(k-1)*nh) = - dd(n+(k-1)*nh)
               enddo   ! n
            endif
         enddo   ! k
c------- diagonal matrix elements of HH and DE in the canonical basis
         do k = 1,nh
            hk = zero
            dk = zero
            do n2 = 1,nh
               h2 = zero
               d2 = zero
               do n1 = 1,nh
                  h2 = h2 + dd(n1+(k-1)*nh)*hh(n1+(n2-1)*nh,m)
                  d2 = d2 + dd(n1+(k-1)*nh)*de(n1+(n2-1)*nh,m)
               enddo
               hk = hk + h2*dd(n2+(k-1)*nh)
               dk = dk + d2*dd(n2+(k-1)*nh)
            enddo
            h(k) = hk
            d(k) = dk
         enddo   ! k
c
c------- reordering according to the energy h(k)
         call ordx(nh,h,d,v2,dd)
c
         do k = 1,nh
            if (v2(k).lt.zero .or. v2(k).gt.2.d0) v2(k) = zero
            if (v2(k).gt.one) v2(k) = one
         enddo ! k
c
         kacan(ib,it)   = klp
         do k=1,nf
            klp=klp+1
            eecan(klp,it)=h(ng+k)
            decan(klp,it)=d(ng+k)
            vvcan(klp,it)=v2(ng+k)
            ibkcan(klp,it)=ib
            do n = 1,nh
               fgcan(n,klp,it)=dd(n+(ng+k-1)*nh)
            enddo
         enddo
         kdcan(ib,it)   = klp - kacan(ib,it)
c
         kacan(ib,it+2) = kla
         do k=1,ng
            kla=kla+1
            eecan(kla,it+2)=h(k)
            decan(kla,it+2)=d(k)
            vvcan(kla,it+2)=v2(k)
            ibkcan(kla,it+2)=ib
            do n=1,nh
               fgcan(n,kla,it+2)=dd(n+(k-1)*nh)
            enddo
         enddo
         kdcan(ib,it+2) = kla - kacan(ib,it+2)

c------- printout for particles
         if (lpr) then
         if (ib.eq.1) e0 = h(ng+1)
         k1 = kacan(ib,it)+1
         k2 = kacan(ib,it)+kdcan(ib,it)
         if (kb(ib) .gt. 0) ip=1
         if (kb(ib) .lt. 0) ip=2
         do k = k1,k2
            e1 = eecan(k,it)
            v1 = vvcan(k,it)
            d1 = decan(k,it)
c           search for the main oscillator component
            smax = zero
            !ab
            imax = 1;
            !ab
            do i = 1,nf
               s = abs(fgcan(i,k,it))
               if (s.gt.smax) then
                  smax = s
                  imax = i
               endif
            enddo
            fx = fgcan(imax,k,it)**2
c
            write(l6,102) k,tb(ib),tt(i0f+imax),fx,e1,v1,d1
            enddo  ! k
          endif   ! lpr
c======================================================================c
   10 enddo      ! ib      end loop over the blocks
c======================================================================c
         nkcan(it)   = klp
         nkcan(it+2) = kla
c======================================================================c
      enddo      ! it      end loop over neutrons and protons --------c
c======================================================================c
c
      if (lpr) then
      write(l6,*) '****** END CANON ***********************************'
      endif
c
  100 format(' single-particle energies and gaps ',1x,
     &       'in the canonical basis: ',a,/1x,66(1h-))

  101 format(7x,a,a,2x,a,5x,a,5x,a,5x,a)
  102 format(i4,2x,a6,a8,f7.2,5f10.4)
c
      return
C-end-CANON
      end
c======================================================================c
c
      subroutine degen(na,n,ea,dd,bb,eb,eps,zz,z)
c
c======================================================================c
c
c     EA    is a set of partially degenerate eigenvalues of some matrix
c     DD    are the corresponding eigenvectors DD.
c     BB    is a matrix, which is diagonalized in the subspaces of
c           degenerate eigenvales EA.
c     EB    contains the eigenvalues of BB in these subspaces.
c     EPS   determines, to which accuracy the eigenvalues EA are
c           considered to be degenerate
c     ZZ,Z  are auxiliary arrays
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      dimension bb(na,n),dd(na,n),ea(n),eb(n)
      dimension zz(na,n),z(n)
c
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
c---- check for degeneracies
      k1 = 0
      do i1 = 1,n
         k1 = k1 + 1
         if (k1.ge.n) goto 20
         do k2 = k1+1,n
            if (ea(k2)-ea(k1).gt.eps) goto 10
         enddo
   10    me = k2 - 1
         ma = k1
c
c----    diagonalize together with bb
         if (me.gt.k1) then
            m0 = ma - 1
            mm = me - m0
c           write(6,'(3i3,a,/,(10f12.8))') mm,ma,me,
c     &      ' Eigenvalues degenerate:',(ea(i),i=ma,me)
c          write(l6,'(3i3,a,/,(10f12.8))') mm,ma,me,
c    &      ' Eigenvalues degenerate:',(ea(i),i=ma,me)

            do m1 = ma,me
               do k = 1,n
                  s = zero
                  do i = 1,n
                     s = s + dd(i,m1) * bb(i,k)
                  enddo   ! i
                  z(k) = s
               enddo   ! k
               do m2 = ma,me
                  s = zero
                  do k = 1,n
                     s = s + z(k) * dd(k,m2)
                  enddo   ! k
                  zz(m1-m0,m2-m0) = s
               enddo   ! m2
            enddo   ! m1
            call sdiag(na,mm,zz,eb(ma),zz,z,+1)
c           call aprint(1,1,6,1,1,mm,eb,' ',' ','H_x')
c	    call aprint(1,1,6,1,1,mm,eb(ma),' ',' ','H_x')
c           call aprint(1,1,1,na,mm,mm,zz,' ',' ','DD_x')
            do i = 1,n
               do m = ma,me
                  s = zero
                  do l = ma,me
                     s = s + dd(i,l) * zz(l-m0,m-m0)
                  enddo   ! l
                  z(m) = s
               enddo   ! m
               do m = ma,me
                  dd(i,m) = z(m)
               enddo   ! m
            enddo   ! i
            k1 = me
         endif
      enddo   ! i1
c
   20 return
c-end-DEGEN
      end

c======================================================================c
c
      subroutine centmas(lpr)
c
c======================================================================c
      use dirqfampar;
      use centmasmod;
      use blodir;
      use bloosc;
      use blokap;
      use gaussl;
      use gfvmod;
      use herpol;
      use lagpol;
      use quaosc;
      use eeecan;
      use blocan;

      implicit real*8(a-h,o-z)
c
!ab   include 'dirqfam.par'
!ab   parameter (MG4 = 4*MG)
c
      logical lpr,lpar
c
!ab   character tb*6                                            ! blokap
!ab   character tt*8                                            ! quaosc
      character nucnam*2                                        ! nucnuc

!ab   dimension wc(MG4,KX),dz(MG4,KX),dr(MG4,KX)
c
      common /baspar/ hom,hb0,b0
!ab   common /blodir/ ka(NBX,4),kd(NBX,4)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
      common /centma/ cmas(3)
      common /defbas/ beta0,q,bp,bz
!ab   common /gaussl/ xl(0:NGL),wl(0:NGL),sxl(0:NGL),rb(0:NGL)
!ab   common /gfviv / iv(-IGFV:IGFV)
!ab   common /herpol/ qh(0:NZX,0:NGH),qh1(0:NZX,0:NGH)
!ab   common /lagpol/ ql(0:2*NRX,0:mlx,0:NGL),ql1(0:2*NRX,0:MLX,0:NGL)
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ hbc,alphi,r0
      common /masses/ amu,ames(4)
      common /nucnuc/ amas,npr(3),nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
!ab   common /quaosc/ nt,nz(NTX),nr(NTX),ml(NTX),ms(NTX),np(NTX),tt(NTX)
!ab   common /eeecan/ eecan(KX,4),decan(KX,4),vvcan(KX,4),
!ab  &                fgcan(NHX,KX,4),ibkcan(KX,4)
!ab   common /blocan/ kacan(nbx,4),kdcan(nbx,4),nkcan(4)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp

c
      ihl(ih,il) = 1+ ih + il*(NGH+1)
c
      if (icm.ne.2) return
c
c
      if (lpr) then
      write(l6,*) '****** BEGIN CENTMAS *******************************'
      endif
c
c======================================================================c
c     calculation of the wave-functions in r-space
c======================================================================c
      az = one/(b0*bz)
      ap = one/(b0*bp)
      n1 = 0
      n2 = MG
      n3 = 2*MG
      n4 = 3*MG
c
      do it = 1,2
      do k  = 1,nkcan(it)
         ib = ibkcan(k,it)
         kk = iabs(kb(ib))
         nf = id(ib,1)
         ng = id(ib,2)
	     nd = nf + ng
         i0 = ia(ib,1)
c
         do il = 1,NGL
         do ih = 1,NGH
            i  = ihl(ih,il)
c
c           large components
            sp  = zero
            srp = zero
            szp = zero
            sm  = zero
            srm = zero
            szm = zero
            do n = 1,nf
               nzz = nz(n+i0)
               nrr = nr(n+i0)
               mll = ml(n+i0)
               mss = iv(kk-mll+1)
c
               fg=fgcan(n,k,it)
               if (mss.eq.1) then
                  sp  = sp  +fg*qh(nzz,ih)*ql(nrr,mll,il)
                  szp = szp +fg*qh1(nzz,ih)*ql(nrr,mll,il)
                  srp = srp +fg*qh(nzz,ih)*ql1(nrr,mll,il)
               else
                  sm  = sm  +fg*qh(nzz,ih)*ql(nrr,mll,il)
                  szm = szm +fg*qh1(nzz,ih)*ql(nrr,mll,il)
                  srm = srm +fg*qh(nzz,ih)*ql1(nrr,mll,il)
               endif
            enddo  ! n
            wc(i+n1,k) = sp
            dz(i+n1,k) = szp*az
            dr(i+n1,k) = srp*ap
            wc(i+n2,k) = sm
            dz(i+n2,k) = szm*az
            dr(i+n2,k) = srm*ap
c
c           small components
            sp  = zero
            srp = zero
            szp = zero
            sm  = zero
            srm = zero
            szm = zero
            do n = nf+1,nd
               nzz = nz(n+i0)
               nrr = nr(n+i0)
               mll = ml(n+i0)
               mss = iv(kk-mll+1)
c
               fg=fgcan(n,k,it)
               if (mss.eq.1) then
                  sp  = sp  +fg*qh(nzz,ih)*ql(nrr,mll,il)
                  szp = szp +fg*qh1(nzz,ih)*ql(nrr,mll,il)
                  srp = srp +fg*qh(nzz,ih)*ql1(nrr,mll,il)
               else
                  sm  = sm  +fg*qh(nzz,ih)*ql(nrr,mll,il)
                  szm = szm +fg*qh1(nzz,ih)*ql(nrr,mll,il)
                  srm = srm +fg*qh(nzz,ih)*ql1(nrr,mll,il)
               endif
            enddo  ! n
            wc(i+n3,k) = sp
            dz(i+n3,k) = szp*az
            dr(i+n3,k) = srp*ap
            wc(i+n4,k) = sm
            dz(i+n4,k) = szm*az
            dr(i+n4,k) = srm*ap
c
         enddo ! il
         enddo ! ih
      enddo ! k

      fac = hbc/(two*amu*amas)
c
c======================================================================c
c     direct term
c======================================================================c
      sz = zero
      sr = zero
      do k = 1,nkcan(it)
         ib = ibkcan(k,it)
         kk = iabs(kb(ib))
c
c        loop over large and small components
         do ifg = 1,2
c
c        loop over spin
         do mss = 1,2
            if (mss.eq.1) then
               mll = kk-1
            else
               mll = kk
            endif
            ipos = (mss-1 + 2*(ifg-1))*MG
c
            sz1 = zero
            sr1 = zero
            do i = 1,MG
               sz1 = sz1 - dz(ipos+i,k)**2
               sr1 = sr1 - dr(ipos+i,k)**2
            enddo
c
            s = zero
            do il = 1,NGL
            do ih = 1,NGH
               i  = ihl(ih,il)
               s = s - wc(ipos+i,k)**2/rb(il)**2
            enddo
            enddo
            sr1 = sr1 + s*mll**2
c
            sz = sz + sz1*vvcan(k,it)
            sr = sr + sr1*vvcan(k,it)
         enddo   ! mss
         enddo   ! ifg
      enddo   ! k
      sdz = two * sz * fac
      sdr = two * sr * fac
      sd  = sdz + sdr
c
c
c======================================================================c
c     exchange term
c======================================================================c
      sz = zero
      sr = zero
c
c---- loop block 1
c---- loop block 2
      do ib1 = 1,nb
         kk1 = iabs(kb(ib1))
         ip1 = kb(ib1)/kk1
         ka1 = kacan(ib1,it) + 1
         ke1 = kacan(ib1,it) + kdcan(ib1,it)
      do ib2 = 1,nb
         kk2 = iabs(kb(ib2))
         ip2 = kb(ib2)/kk2
         ka2 = kacan(ib2,it) + 1
         ke2 = kacan(ib2,it) + kdcan(ib2,it)
         lpar = (ip1.ne.ip2)
c
         if (kk1 .eq. kk2) then
            do k1 = ka1,ke1
               v1 = sqrt(vvcan(k1,it))
               u1 = sqrt(1.d0-v1**2)
            do k2 = ka2,ke2
               v2 = sqrt(vvcan(k2,it))
               u2 = sqrt(1.d0-v2**2)
               uv = v1*v2*(v1*v2 + u1*u2)
                  if(lpar) then
                  s = zero
                  do ifg = 1,2
                  do mss  = 1,2
                     ipos = (mss-1 + 2*(ifg-1))*MG
                     do i = 1,MG
                        s = s + wc(ipos+i,k1)*dz(ipos+i,k2)
                     enddo
                  enddo   ! mss
                  enddo   ! ifg
                  sz = sz + s**2 * uv
               endif
            enddo  ! k2
            enddo  ! k1
         elseif (kk1 .eq. kk2+1) then
            do k1 = ka1,ke1
               v1 = sqrt(vvcan(k1,it))
               u1 = sqrt(1.d0-v1**2)
            do k2 = ka2,ke2
               v2 = sqrt(vvcan(k2,it))
               u2 = sqrt(1.d0-v2**2)
               uv = v1*v2*(v1*v2 + u1*u2)
                if(lpar) then
                  s = zero
                  do ifg = 1,2
                  do mss  = 1,2
                     if (mss.eq.1) then
                        ml2 = kk2 - 1
                     else
                        ml2 = kk2
                     endif
                     ipos = (mss-1 + 2*(ifg-1))*MG
                     do il = 1,NGL
                     do ih = 1,NGH
                        i = ihl(ih,il)
                        s = s + wc(ipos+i,k1)*dr(ipos+i,k2)
                        s = s - wc(ipos+i,k1)*wc(ipos+i,k2) *
     &                          ml2/rb(il)
                     enddo   ! ih
                     enddo   ! il
                  enddo   ! mss
                  enddo   ! ifg
                  sr = sr + half*s**2 * uv
               endif
            enddo   ! k2
            enddo   ! k1
         elseif (kk1 .eq.  kk2-1) then
            do k1 = ka1,ke1
               v1 = sqrt(vvcan(k1,it))
               u1 = sqrt(1.d0-v1**2)
            do k2 = ka2,ke2
               v2 = sqrt(vvcan(k2,it))
               u2 = sqrt(1.d0-v2**2)
               uv = v1*v2*(v1*v2 + u1*u2)
                if(lpar) then
                  s = zero
                  do ifg = 1,2
                  do mss  = 1,2
                     if (mss.eq.1) then
                        ml2 = kk2 - 1
                     else
                        ml2 = kk2
                     endif
                     ipos = (mss-1 + 2*(ifg-1))*MG
                     do il = 1,NGL
                     do ih = 1,NGH
                        i  = ihl(ih,il)
                        s = s+wc(ipos+i,k1)*dr(ipos+i,k2)
                        s = s+wc(ipos+i,k1)*wc(ipos+i,k2) *
     &                      ml2/rb(il)
                     enddo   ! ih
                     enddo   ! il
                  enddo   ! mss
                  enddo   ! ifg
                  sr = sr + half*s**2 * uv
               endif
            enddo   ! k2
            enddo   ! k1
         endif
         if ((kk1.eq.1) .and. (kk2.eq.1)) then
            do k1 = ka1,ke1
               v1 = sqrt(vvcan(k1,it))
               u1 = sqrt(1.d0-v1**2)
            do k2 = ka2,ke2
               v2 = sqrt(vvcan(k2,it))
               u2 = sqrt(1.d0-v2**2)
               uv = v1*v2*(v1*v2 + u1*u2)
                if(lpar) then
                  s = zero
                  do ifg = 1,2
                  do mss  = 1,2
                     if (mss.eq.1) then
                        ml2 = kk2
                        ms2 = 2
                     else
                        ml2 = kk2-1
                        ms2 = 1
                     endif
                     ipos1 = (mss -1 + 2*(ifg-1))*MG
                     ipos2 = (ms2-1 + 2*(ifg-1))*MG
		     s1 = zero
                     do il = 1,NGL
                     do ih = 1,NGH
                        i1 = ipos1 + ihl(ih,il)
                        i2 = ipos2 + ihl(ih,il)
                        s1 = s1 + wc(i1,k1)*dr(i2,k2) +
     &                            wc(i1,k1)*wc(i2,k2)*ml2/rb(il)
                     enddo
                     enddo
                     s = s - (3-2*mss)*(3-2*ifg)*s1
                  enddo   ! mss
                  enddo   ! ifg
                  sr = sr + half*s**2 * uv
               endif
            enddo   ! k2
            enddo   ! k1
         endif
      enddo   ! ib2
      enddo   ! ib1
      sez = two * sz * fac
      ser = two * sr * fac
      se  = ser + sez
      cmas(it) = sd + se
c
      if (lpr) then
         write(l6,100) ' centmas z',it,sdz,sez,sdz+sez
         write(l6,100) ' centmas r',it,sdr,ser,sdr+ser
         write(l6,100) ' centmas  ',it,sd,se,cmas(it),0.375d0*hom
      endif
c
      enddo   ! it
      cmas(3) = cmas(1) + cmas(2)
c      write(l6,100) ' centmas T',3,cmas,0.75d0*hom
c
      if (lpr) then
      write(l6,*) '****** END CENTMAS *********************************'
      endif
c
  100 format(a,i3,3f12.7,f15.7)
c
c      CALL zeit(2,'CENTMAS_',6,.false.)

      return
c-end-CENTMAS
      end

c======================================================================c

      subroutine coulom(lpr)

c======================================================================c
c
c     calculation of the Coulomb field
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use coulmb;
      use procou;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
!ab   common /coulmb/ cou(MG),drvp(MG)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
!ab   common /procou/ ggc(MG,MG)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      if (icou.eq.0) return
c
c      CALL zeit(1,'COULOM__',6,.false.)
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN COULOM *******************************'
c
      do i = 1,MG
         s = zero
	     do k = 1,MG
	        s = s + ggc(i,k) * drvp(k)
	     enddo   ! k
	     cou(i) = s
      enddo   ! i
c
c      if (lpr) call prigh(1,cou,one,'Coulom')
c
      if (lpr)
     &write(l6,*) ' ****** END COULOM *********************************'
c
c      CALL zeit(2,'COULOM__',6,.false.)
c
      return
C-end-COULOM
      end
c======================================================================c

      subroutine greecou(lpr)

c======================================================================c
c
C     calculation of the Coulomb-propagator
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use procou;
      use gaussh;
      use gaussl;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
!ab   common /procou/ ggc(MG,MG)
!ab   common /gaussh/ xh(0:NGH),wh(0:NGH),zb(0:NGH)
!ab   common /gaussl/ xl(0:NGL),wl(0:NGL),sxl(0:NGL),rb(0:NGL)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
c
      ihl(ih,il) = 1+ih + il*(NGH+1)
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN GREECOU ******************************'
c
c
      if (icou.gt.0) then
         f = one/(alphi*two*pi)
         do kl = 0,NGL
            r   = rb(kl)
            r4r = -4*r
         do kh = 0,NGH
            z   = zb(kh)
c
            do il = 0,NGL
               r1  = rb(il)
               r2r = (r+r1)**2
               r4s = r4r*r1
            do ih = 0,NGH
               z1  = zb(ih)
c
               s   = zero
               ivm = 1
               do m = 1,2
                  y  = (z+ivm*z1)**2+r2r
                  x  = one + r4s/y
                  xx = x*x
c
c                 expansion of the eliptic integral
                  e = one + 0.44325141463*x + 0.06260601220*xx
     &                + 0.04757383546*xx*x + 0.01736506451*xx*xx
                  if (x.gt.1d-9) then
                     e = e - dlog(x) *
     &                  ( 0.24998368310*x     + 0.09200180037*xx
     &                   +0.04069697526*xx*x + 0.00526449639*xx*xx)
                  endif
                  s = s + e*sqrt(y)
                  ivm = -ivm
               enddo   ! m
c
               ggc(ihl(ih,il),ihl(kh,kl)) =  f*s
c
            enddo   ! ih
            enddo   ! il
         enddo   ! kh
         enddo   ! kl
c
c        if (lpr)
c    &   call aprint(1,1,6,MG,MG,MG,ggc,' ',' ','VC')
      endif   ! icou
c
      if (lpr)
     &write(l6,*) ' ****** END GREECOU ********************************'
      return
c-end-GREECOU
      end

c======================================================================c

      subroutine cstrpot(lpr)
c======================================================================c
c
c     CALCULATION OF THE CONSTRAINING POTENTIAL for axial case
c
c we consider an constraining operator Q
c and require a fixed expectation value      <Q>
c actual value of <Q> = Tr(Q*ro) (density)   <Q>  = d_Q   = d_r2
c required (wanted) value of <Q>                    w_Q   = w_r2
c strength parameter of the constraint              x_Q   = x_r2
c
c constraining field:                        V    -> V - x_Q * Q
c
c augmented constraint (ic_Q = 3):
c   final    wanted value               w_Q (from common /con_Q/
c   variable wanted value               v_Q
c   initial  wanted value               v_Q = w_Q
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use constr;
      use gaussh;
      use gaussl;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
      common /baspar/ hom,hb0,b0
!ab   common /constr/ vc(MG,2)
      common /con_b2/ betac,q0c,cquad,c0,alaq,calcq0,icstr
      common /con_b2a/ fac0
!ab   common /gaussh/ xh(0:NGH),wh(0:NGH),zb(0:NGH)
!ab   common /gaussl/ xl(0:NGL),wl(0:NGL),sxl(0:NGL),rb(0:NGL)
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      fac0=zero
      if (icstr.eq.0) return
c
      if (lpr) then
      write(l6,*) '****** BEGIN CSTRPOT *******************************'
      endif
c
      alaq = alaq+c0*(q0c-calcq0)
      q0cst= q0c+alaq/c0
      fac0 = -c0*(q0cst-calcq0)
c
c---- calculation of new constraining fields VC
      do il = 0,NGL
         r = rb(il)
         rr = r*r
         do ih = 0,NGH
            i = ih+1 + il*(NGH+1)
            z = zb(ih)
            zz = z*z
            v = fac0*(2*zz - rr)
            vc(i,1) = v
            vc(i,2) = v
         enddo   ! ih
      enddo   ! il

c
c      if (lpr) then
c         call prigh(0,sig,one,'V_constr')
c         call prigh(1,vc,hbc,'V_c')
c      endif
c
      if (lpr) then
      write(l6,*) '****** END CSTRPOT *********************************'
      read*
!ab   l6 = lx
      endif
c
      return
c-end-CSTRPOT
      end
c======================================================================c

      subroutine cstrx(ic,x,c,w,d,slam,t)

c======================================================================c
      implicit real*8 (a-h,o-z)
c
      character t*(*)
c
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      data iit/800/
c
      x = zero
      v = zero
      if (ic.eq.0) return
      if (ic.eq.1) then
         x = c * (v - d)
      elseif (ic.eq.2) then
         x = c * (w - d)
      elseif (ic.eq.3) then
         if (ii.eq.1) then
            slam = zero
         else
            slam = slam + c*(w - d)
         endif   ! ii
         v = w + slam/c
         if (ii.gt.iit) ic = 1
         x = c * (v - d)
      else
         stop 'in CSTRC: ic not properly defined'
      endif
c
      write(6   ,100) ii,t,ic,w,d,v,x,si

!ab   write(lstr,100) ii,t,ic,w,d,v,x,si

c     read*
  100 format(i4,2x,a,i3,5f10.6,f16.6)
c
      x = x/hbc
c
      return
c-end-CSTRX
      end

c======================================================================c

      subroutine default(lpr)

c======================================================================c
c
c     Default for Relativistic Mean Field spherical
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use constr;
      use couplf;
      use broyde;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
      character tp*1,tis*1,tit*8,tl*1                   ! textex

      common /baspar/ hom,hb0,b0
!ab   common /constr/ vc(MG,2)
!ab   common /couplf/ ff(MG,4,2)
      common /couplg/ ggmes(4),lmes(4)
      common /dforce/ a_m(4),b_m(4),c_m(4),d_m(4),dsat
      common /fermi / ala(2),tz(2)
      common /initia/ inin,inink
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
!ab   common /broyde2/ ibroyd
c


c======================================================================c
c---- signs and factorials
c-----------------------------------------------------------------------
      call gfv(n)
      if (n.ne.IGFV) stop 'in DEFAULT: IGFV wrong for GFV'
      call binom(n)
      if (n.ne.IGFV) stop 'in DEFAULT: IGFV wrong for BINOM'
c======================================================================c


c======================================================================c
c     center off mass correction:
c-----------------------------------------------------------------------
c
c     a) the kinetic energy in the variation can be
c           <T>  = hbc**2/(2*amu)*Delta   (no correction)
c        or <T>' = hbc**2/(2*amu)*(1-1/A)*Delta
c                  the diagonal part of ecm=<P**2>/2Am is included
c     b) the total energy is
c        <T> + <V> - <P**2>/2Am
c
c     icm: 0   variation of <T>
c              hb0 = hbc**2/(2*amu)       ecm = 3/4*hom
c          1   variation of <T>'
c              hb0 = hb0*(1-1/A)          ecm = 3/4*hom
c          2   variation of <T>
c              hb0 = hbc**2/(2*amu)       ecm = <P**2>/2M
c
      icm = 0

c
c======================================================================c






      !ab
c=======================================================================c
c=======================================================================c
c=======================================================================c
      intape = 123;
      open( intape , file = 'dirqfam.dat' , status = 'old' );
      read( intape , * );read( intape , * );read( intape , * );
      read( intape , * );read( intape , * );read( intape , * );
      read( intape , * );read( intape , * );read( intape , * );
      read( intape , * );read( intape , * );read( intape , * );
      read( intape , * );read( intape , * );read( intape , * );
      read( intape , * );read( intape , * );read( intape , * );
      read( intape , * );read( intape , * );read( intape , * );
      read( intape , '(18x,i9)' ) icou;
      close( intape );
c=======================================================================c
c=======================================================================c
c=======================================================================c
      !ab


c======================================================================c
c     Coulomb-Field
c-----------------------------------------------------------------------
c     icou: Coulomb-field:  0   not at all
c                           1   only direct term
c                           2   plus exchange


c======================================================================c

c======================================================================c
c     pairing
c----------------------------------------------------------------------c
      do it = 1,2
         del(it)  = zero
         spk(it)  = zero
         spk0(it) = zero
         ala(it)  = -7.0d0     ! chemical potential
      enddo   ! it
c======================================================================c

c======================================================================c
c     iteration
c----------------------------------------------------------------------c
c
      maxi = 500             ! maximal number of iteration
c     maxi = 1               ! remove
      si   = one             ! actual error in the main iteration
      epsi = 1.d-6           ! accuracy for the main iteration
      iaut  = 1              ! automatic change of xmix: 0 (no) 1 (yes)
      inxt  = 3              ! number of iterations until next question
      xmix  = 0.5d0          ! starting xmix
      xmax  = 0.7d0          ! maximal value for xmix
      xmix0 = xmix           ! basic xmix, where it returns to
      ibroyd= 1
c======================================================================c

c======================================================================c
c---- parameters of the initial potentials
c----------------------------------------------------------------------c
c     inin = 0: fields read, 1: saxon-wood,
      inin  = 1
c     inink = 0: pairing potential read, 1: pairing potential monopol
      inink = 1
c
c     oscillator length b0 (is calcuated for b0 <= 0)
      b0 = -2.320
c======================================================================c

c======================================================================c
c     constraining fields
c----------------------------------------------------------------------c
      icstr=0
      alaq=0.d0
      call mzero(MG,MG,2,vc)
c======================================================================c

c======================================================================c
c---- tapes
c----------------------------------------------------------------------c
      l6   = 10
      lin  = 3
      lou  = 6
      lwin = 1
      lwou = 2
      lplo = 11
      laka = 12
      lvpp = 13
c======================================================================c

c======================================================================c
c---- preparation of density dependence
c----------------------------------------------------------------------c
      do m = 1,4
         lmes(m) = 0
         ggmes(m) = zero
         a_m(m) = zero
         b_m(m) = zero
         c_m(m) = zero
         d_m(m) = zero
         do i = 1,MG
            ff(i,m,1) = one
	    ff(i,m,2) = zero
         enddo   ! i
      enddo   ! m
c======================================================================c
c
      return
c-end-DEFAULT
      end
c======================================================================c

      blockdata block1

c======================================================================c
      use dirqfampar;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'

      character tp*1,tis*1,tit*8,tl*1                   ! textex
      common /physco/ hbc,alphi,r0
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)

c---- fixed texts
      data tp/'+','-'/,tis/'n','p'/,tit/'Neutron:','Proton: '/
      data tl/'s','p','d','f','g','h','i','j','k','l','m',
     &            'n','o','P','q','r','S','t','u','v','w',
     &            'x','y','z','0','0','0','0','0','0','0'/
c
c---- physical constants
      data hbc/197.328284d0/,r0/1.2d0/,alphi/137.03602/

c-end-BLOCK1
      end
C======================================================================c

      subroutine delta(it,lpr)

c======================================================================c
c
c     calculats the pairing field
c     for separable pairing: Tian,Ma,Ring, PRB 676, 44 (2009)
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use deltamod;
      use blokap;
      use bloosc;
      use quaosc;
      use deldel;
      use rokaos;
      use vvvikf;
      use tmrwnn;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical  lpr
c
!ab   character tb*6                                            ! blokap
!ab   character tt*8                                            ! quaosc
      character tp*1,tis*1,tit*8,tl*1                           ! textex
c
!ab   dimension pnn(NNNX)
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
!ab   common /quaosc/ nt,nnn(NTX,5),tt(NTX)
!ab   common /deldel/ de(NHHX,NB2X)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
!ab   common /rokaos/ rosh(NHHX,NB2X),aka(MVX,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
!ab   common /vvvikf/ mv,ipos(NBX),nib(MVX),nni(2,MVX)
!ab   common /tmrwnn/ wnn(MVX,NNNX),nnmax
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
      common /tmrpar/ gl(2),ga

      if (lpr) then
      write(l6,*) ' ****** BEGIN DELTA ********************************'
      endif
c
c----------------------------------------------------------------------c
c     calculation of PNN
c----------------------------------------------------------------------c
      do nn = 1,nnmax
c        nx = 5
c        call aprint(3,1,1,nh,nx,nx,aka(1,it),' ',' ','KA++')
c        call aprint(3,1,1,nh,nx,nx,wnn(1,nn),' ',' ','WNN++')
         s = zero
         do i = 1,mv
            s = s + wnn(i,nn)*aka(i,it)
         enddo   ! i
         pnn(nn) = s
      enddo  ! nn
c     if (lpr) then
c        write(6,*) 'PNN ',mv,nnmax
c        do nn = 1,nnmax
c           write(6,101) ' nn =',nn,' PNN =',pnn(nn)
c        enddo   ! nn
c        write(l6,100) 'pnn',nnmax,(pnn(nn),nn=1,8)
c     endif
c
      g = half*gl(it)
      i12 = 0
      do ib = 1,nb
         i0 = ia(ib,1)
         nf = id(ib,1)
         ng = id(ib,2)
         nh = nf + ng
         m  = ib + (it-1)*NBX
         do n2 =  1,nf
         do n1 = n2,nf
            i12 = i12 + 1
            s = zero
            do nn = 1,nnmax
               s = s + wnn(i12,nn)*pnn(nn)
            enddo   ! nn
            de(n1+(n2-1)*nh,m) = -g*s
            de(n2+(n1-1)*nh,m) = -g*s
         enddo  ! n1
         enddo  ! n2
c
      if (lpr.and.ib.eq.1) then
         k0 = ia(ib,1)+1
         nx = nf
         nx = 5
         call aprint(1,3,1,nh,nx,nx,de(1,m),tt(k0),tt(k0),'DE++')
      endif
      enddo   ! ib
c
      if (lpr) then
      write(l6,*) ' ****** END DELTA **********************************'
      endif
c
  100 format(a,i6,8f10.6)
  101 format(a,i5,a,f15.6)
c
      return
C-end-DELTA
      end
c======================================================================c

      subroutine densit(lpr)

c=====================================================================c
C
c     calculates the densities in r-space at Gauss-meshpoints
C
c---------1---------2---------3---------4---------5---------6---------7-
      use dirqfampar;
      use densitmod;
      use blokap;
      use bloosc;
      use coulmb;
      use dens;
      use gaucor;
      use gaussh;
      use gaussl;
      use gfvmod;
      use herpol;
      use lagpol;
      use quaosc;
      use rhorho;
      use rokaos;

      implicit real*8 (a-h,o-z)
!ab   include 'dirqfam.par'
c
      logical lpr
c
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character nucnam*2                                        ! nucnuc
!ab   character tb*6                                            ! blokap
!ab   character tt*8                                            ! quaosc
c
      dimension rsh(2)
!ab   dimension drs(MG,2),drv(MG,2)
c
      common /baspar/ hom,hb0,b0
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
!ab   common /coulmb/ cou(MG),drvp(MG)
      common /defbas/ beta0,q,bp,bz
!ab   common /dens  / ro(MG,4),dro(MG,4)
!ab   common /gaucor/ ww(MG)
!ab   common /gaussh/ xh(0:NGH),wh(0:NGH),zb(0:NGH)
!ab   common /gaussl/ xl(0:NGL),wl(0:NGL),sxl(0:NGL),rb(0:NGL)
!ab   common /gfviv / iv(-IGFV:IGFV)
!ab   common /herpol/ qh(0:NZX,0:NGH),qh1(0:NZX,0:NGH)
!ab   common /lagpol/ ql(0:2*NRX,0:MLX,0:NGL),ql1(0:2*NRX,0:MLX,0:NGL)
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,npr(3),nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
!ab   common /quaosc/ nt,nz(NTX),nr(NTX),ml(NTX),ms(NTX),np(NTX),tt(NTX)
!ab   common /rhorho/ rs(MG,2),rv(MG,2)
!ab   common /rokaos/ rosh(NHHX,NB2X),aka(MVX,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
c
      ihl(ih,il) = 1+ih +il*(NGH+1)
c
      if (lpr) then
      write(l6,*) ' ****** BEGIN DENSIT *******************************'
      endif
c
      ap2  = one/(b0*bp)**2
      az2  = one/(b0*bz)**2
c
      do it = 1,2
      do i = 1,MG
         rs(i,it)  = zero
         rv(i,it)  = zero
         drs(i,it) = zero
         drv(i,it) = zero
      enddo   ! i
      enddo   ! it
c
c---- loop over K-parity-blocks
      il = 0
      do ib = 1,nb
         nf = id(ib,1)
         ng = id(ib,2)
         nh = nf + ng
c
c------- loop over contributions from large and small components
         do ifg = 1,2
            n0  = (ifg-1)*nf
            nd  = id(ib,ifg)
            i0  = ia(ib,ifg)
            ivv = iv(ifg)
c
c------- loop of oscillator basis states n2 and n1
         do n2 =  1,nd
            nz2 = nz(n2+i0)
            nr2 = nr(n2+i0)
            ml2 = ml(n2+i0)
         do 20 n1 =  n2,nd
            nz1 = nz(n1+i0)
            nr1 = nr(n1+i0)
            ml1 = ml(n1+i0)
            if (ml1.ne.ml2) goto 20
            ll = ml1**2
            i12 = 2 - n2/n1
            rsh(1) = rosh(n0+n1+(n0+n2-1)*nh,ib)*i12
            rsh(2) = rosh(n0+n1+(n0+n2-1)*nh,ib+NBX)*i12

c-------    loop over the mesh-points
            ben  = az2*(nz1+nz2+1) + ap2*2*(nr1+nr2+ml1+1)
            do il = 0,NGL
               qlab  = ql(nr1,ml1,il)*ql(nr2,ml1,il)
               qltab = ap2*(ql1(nr1,ml1,il)*ql1(nr2,ml1,il) +
     &                      ll*qlab/xl(il))
               bfl   = ap2*xl(il) - ben
c
            do ih = 0,NGH
               qhab  = qh(nz1,ih)*qh(nz2,ih)
               qh1ab = az2*qh1(nz1,ih)*qh1(nz2,ih)
               sro   = qlab*qhab
               stau  = qh1ab*qlab + qhab*qltab
               sdro  = sro*(az2*xh(ih)**2 + bfl)
c
c-------       scalar and vector and pairing density
               do it = 1,2
                  fgr = rsh(it)*sro
                  rs(ihl(ih,il),it) = rs(ihl(ih,il),it) - ivv*fgr
                  rv(ihl(ih,il),it) = rv(ihl(ih,il),it) + fgr
c
c-------          delta rho
                  sdt = 2*rsh(it)*(sdro+stau)
                  drs(ihl(ih,il),it) = drs(ihl(ih,il),it) - ivv*sdt
                  drv(ihl(ih,il),it) = drv(ihl(ih,il),it) + sdt
               enddo   ! it
            enddo   ! ih
            enddo   ! il

   20    continue   ! n1
         enddo   ! n2
         enddo   ! ifg
      enddo   ! ib
c
c---- check, whether integral over dro vanishes
      s1 = zero
      s2 = zero
      do i = 1,MG
         s1 = s1 + drv(i,1)
         s2 = s2 + drv(i,2)
      enddo
      if (lpr) write(l6,*) ' integral over dro',s1,s2
c
c
c---- normalization and renormalization to particle number
      do it = 1,2
         s  = zero
         do i = 1,MG
            s  =  s + rv(i,it)
         enddo
         if (lpr) write(l6,'(a,i3,2f15.8)') '  Integral over rv:',it,s
         s  = npr(it)/s
         do i = 1,MG
            rv(i,it)  = s*rv(i,it)
            rs(i,it)  = s*rs(i,it)
            drs(i,it) = s*drs(i,it)
            drv(i,it) = s*drv(i,it)
         enddo
c
c------- printout of the density in configuration space
c         if (lpr) then
c            call prigh(2,rs(1,it),one,'RS  '//tis(it))
c            call prigh(2,rv(1,it),one,'RV  '//tis(it))
c            call prigh(2,drs(1,it),one,'DROS'//tis(it))
c            call prigh(2,drv(1,it),one,'DROV'//tis(it))
c         endif
      enddo  ! it
c
      do i = 1,MG
         f        = one/ww(i)
         ro(i,1)  = f*( + rs(i,1) + rs(i,2) )
         ro(i,2)  = f*( + rv(i,1) + rv(i,2) )
         ro(i,3)  = f*( - rs(i,1) + rs(i,2) )
         ro(i,4)  = f*( - rv(i,1) + rv(i,2) )
         dro(i,1) = f*(   drs(i,1) + drs(i,2) )
         dro(i,2) = f*(   drv(i,1) + drv(i,2) )
         dro(i,3) = f*( - drs(i,1) + drs(i,2) )
         dro(i,4) = f*( - drv(i,1) + drv(i,2) )
         drvp(i)  = drv(i,2)
      enddo   ! i
c      if (lpr) then
c         call prigh(1,ro(1,1),one,'RO-sig')
c         call prigh(1,ro(1,2),one,'RO-ome')
c         call prigh(1,ro(1,3),one,'RO-del')
c         call prigh(1,ro(1,4),one,'RO-rho')
c         ix = 5
c         do i = 1,ix
c            write(6,100) i,(ro(i,m),m=1,4)
c         enddo   ! i
c      endif
c
  100 format(i5,4f15.10)
c
      if (lpr) then
      write(l6,*) ' ****** END DENSIT *********************************'
      endif
c
      return
C-end-DENSIT
      end
c======================================================================c

      subroutine denssh(it,lpr)

c=====================================================================c
C
c     calculates the densities in oscillator basis
C
c---------1---------2---------3---------4---------5---------6---------7-
      use dirqfampar;
      use blodir;
      use blokap;
      use bloosc;
      use quaosc;
      use rokaos;
      use waveuv;

      implicit real*8 (a-h,o-z)
!ab   include 'dirqfam.par'
c
      logical lpr
c
!ab   character tb*6                                            ! blokap
!ab   character tt*8                                            ! quaosc
      character tp*1,tis*1,tit*8,tl*1                           ! textex
c
!ab   common /blodir/ ka(NBX,4),kd(NBX,4)
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
!ab   common /quaosc/ nt,nz(NTX),nr(NTX),ml(NTX),ms(NTX),np(NTX),tt(NTX)
!ab   common /rokaos/ rosh(NHHX,NB2X),aka(MVX,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
!ab   common /waveuv/ fguv(NHBX,KX,4),equ(KX,4)
c
      if (lpr) then
      write(l6,*) ' ****** BEGIN DENSSH *******************************'
      endif
c
c
      sp  = zero
      il = 0
c=====================================================================c
      do ib = 1,nb            ! loop over the blocks
c=====================================================================c
         nf  = id(ib,1)
         ng  = id(ib,2)
         nh  = nf + ng
         k1  = ka(ib,it) + 1
         ke  = ka(ib,it) + kd(ib,it)
         k1a = ka(ib,it+2) + 1
         kea = ka(ib,it+2) + kd(ib,it+2)
         mul = mb(ib)
         m   = ib + (it-1)*NBX
         if (lpr.and.ib.eq.1) write(l6,'(/,a,1x,a)') tb(ib),tis(it)
c
c------- calculation of ro
         do n2 =   1,nh
         do n1 =  n2,nh
            sr = zero
            do k = k1,ke
               sr = sr + fguv(nh+n1,k,it)*fguv(nh+n2,k,it)
            enddo    ! k
            do k = k1a,kea                 ! no-sea approximation
               sr = sr + fguv(nh+n1,k,it+2)*fguv(nh+n2,k,it+2)
            enddo                          ! no-sea approximation
            sr = mul*sr
            rosh(n1+(n2-1)*nh,m) = sr
            rosh(n2+(n1-1)*nh,m) = sr
         enddo   ! n1
         enddo   ! n2
c
c------- contributions of large components f*f to kappa
         i0  = ia(ib,1)
         do n2 =   1,nf
         do n1 =  n2,nf
            ml1=ml(i0+n1)
            ml2=ml(i0+n2)
            i12 = 2 - n2/n1
            il  = il + 1
            sk = zero
            do k = k1,ke
               sk = sk + fguv(nh+n1,k,it)*fguv(n2,k,it)
            enddo    ! k
            do k = k1a,kea                   ! no-sea approximation
               sk = sk + fguv(nh+n1,k,it+2)*fguv(n2,k,it+2)
            enddo                            ! no-sea approximation
            sk = mul*sk
            if (ml1.ne.ml2) sk = zero      ! remove question ???
            aka(il,it) = i12*sk
            if (n1.eq.n2) sp = sp + aka(i12,it)
         enddo   ! n1
         enddo   ! n2
         spk(it)=sp
c
      if (lpr.and.ib.eq.1) then
         i0 = ia(ib,1) + 1
         write(l6,100) 'block ',ib,': ',tb(ib),' ',tit(it)
         nx = nh
         nx = 5
         call aprint(1,3,1,nh,nx,nx,rosh(1,m),tt(i0),tt(i0),'RO')
         nx = nf
         nx = 5
         call aprint(3,3,1,nh,nx,nx,aka(1,it),tt(i0),tt(i0),'AK')
      endif
c
      enddo   ! ib
c
      if (lpr) then
      write(l6,*) ' ****** END DENSSH *********************************'
      endif
c
  100 format(//,a,i2,4a)
c
      return
C-end-DENSSH
      end

c======================================================================c

      subroutine dinout(is,lpr)

c======================================================================c
c
c     IS = 1 : reads pairing field Delta from tape
c     IS = 2 : writes pairing field Delta to tape
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use blokap;
      use bloosc;
      use quaosc;
      use deldel;
      use ptenso;
      use vvvikf;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
      logical lpr
c
!ab   character tt*8                                            ! quaosc
!ab   character tb*6                                            ! blokap
c
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
!ab   common /quaosc/ nt,nz(NTX),nr(NTX),ml(NTX),ms(NTX),np(NTX),tt(NTX)
!ab   common /deldel/ de(NHHX,NB2X)
      common /initia/ inin,inink
!ab   common /ptenso/ aka(MVX,2),dkd(MVX,2)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
!ab   common /vvvikf/ mv,ipos(NBX),nib(MVX),nni(2,MVX)
c
      if (is.eq.1.and.inink.ne.0) return

      if (lpr) then
      write(l6,*) '****** BEGIN DINOUT ********************************'
      endif
c
c
c======================================================================c
c---- reads pairing potential
c======================================================================c
      if (is.eq.1) then
         call mzero(NHHX,NHHX,NB2X,de)


      !ab
c=======================================================================c
         open( laka , file = './output/GS_output/dirhb.del' ,
     &         status = 'unknown' , form = 'unformatted' );
c=======================================================================c
      !ab


         read(laka) mv0
         if (mv0.ne.mv) stop 'in DINOUT: mv wrong'
         read(laka) dkd
         do it = 1,2
            il = 0
            do ib = 1,nb
   	       nf = id(ib,1)
   	       ng = id(ib,2)
               nh = nf + ng
               m  = ib + (it-1)*NBX
               do n2 =  1,nf
               do n1 = n2,nf
                  il = il + 1
                  de(n1+(n2-1)*nh,m) = dkd(il,it)
                  de(n2+(n1-1)*nh,m) = dkd(il,it)
               enddo   ! n1
               enddo   ! n2
            enddo   ! ib
         enddo   ! it
         close(laka)
c======================================================================c
c==== writing of the pairing potential
c======================================================================c
      elseif (is.eq.2) then



      !ab
c=======================================================================c
         open( laka , file = './output/GS_output/dirhb.del' ,
     &         status = 'unknown' , form = 'unformatted' );
c=======================================================================c
      !ab


         rewind(laka)
         write(laka) mv

         do it = 1,2             ! loop over neutron, proton
c            write(laka) mv
            il = 0
            do ib = 1,nb
   	       nf = id(ib,1)
   	       ng = id(ib,2)
               nh = nf + ng
               m  = ib + (it-1)*NBX
               do n2 =  1,nf
               do n1 = n2,nf
                  il = il + 1
                  dkd(il,it) = de(n1+(n2-1)*nh,m)
               enddo   ! n1
               enddo   ! n2
            enddo   ! ib
c            call writms(laka,mv,dkd(1,it))
         enddo   ! it
         write(laka) dkd
      endif   ! is
c
      close(laka)
      if(lpr) write(l6,*) 'pairing fields written to dirhb.del'
c
      if (lpr) then
      write(l6,*) ' ****** END DINOUT *********************************'
      endif
c
      return
c-end-DINOUT
      end

c======================================================================c

      subroutine dirhb(it,lpr)

c======================================================================c
c
c     solves the RHB-Equation
c     IT    = 1 for neutrons
c     IT    = 2 for protons
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use dirhbmod;
      use blodir;
      use blokap;
      use bloosc;
      use deldel;
      use gamgam;
      use quaosc;
      use waveuv;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr,lprl
c
      character*1 bbb
!ab   character*8 tbb(NHBX)
      character tp*1,tis*1,tit*8,tl*1                           ! textex
!ab   character tb*6                                            ! blokap
!ab   character tt*8                                            ! quaosc
      character nucnam*2                                        ! nucnuc
c
!ab   dimension hb(NHBQX),e(NHBX),ez(NHBX)
c
!ab   common /blodir/ ka(NBX,4),kd(NBX,4)
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
!ab   common /deldel/ de(NHHX,NB2X)
      common /fermi / ala(2),tz(2)
!ab   common /gamgam/ hh(NHHX,NB2X)
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,npr(3),nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
!ab   common /quaosc/ nt,nnn(NTX,5),tt(NTX)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
!ab   common /waveuv/ fguv(NHBX,KX,4),equ(KX,4)
c
      data maxl/200/,epsl/1.d-8/,bbb/'-'/,lprl/.false./
c
      if (lpr) then
      write(l6,*) ' ****** BEGIN DIRHB ********************************'
      endif

      !ab
      dd    = 0.d0;
      sn    = 0.d0;
      !ab

      dl    = 100.d0
      xh    = ala(it) + dl
      xl    = ala(it) - dl
      al    = ala(it)
c
c======================================================================c
      do lit = 1,maxl         ! loop over lambda-iteration
c======================================================================c
         snold = sn
         sn  = zero
         klp = 0
         kla = 0
c======================================================================c
      do ib = 1,nb            ! loop over differnt blocks
c======================================================================c
         mul  = mb(ib)
         nf   = id(ib,1)
         ng   = id(ib,2)
         nh   = nf + ng
         nhb  = nh + nh
         m    = ib + (it-1)*NBX
c
c------- calculation of the RHB-Matrix:
         do n2 = 1,nh
         do n1 = n2,nh
            hb(   n1+(   n2-1)*nhb) =  hh(n1+(n2-1)*nh,m)
            hb(nh+n1+(nh+n2-1)*nhb) = -hh(n1+(n2-1)*nh,m)
            hb(nh+n1+(   n2-1)*nhb) =  de(n1+(n2-1)*nh,m)
            hb(nh+n2+(   n1-1)*nhb) =  de(n2+(n1-1)*nh,m)
         enddo
            hb(   n2+(   n2-1)*nhb) =  hb(n2+(n2-1)*nhb) - al
            hb(nh+n2+(nh+n2-1)*nhb) = -hb(n2+(n2-1)*nhb)
         enddo
c
c------- Diagonalization:
         if (lpr) then
            i0f = ia(ib,1)
            do n = 1,nh
               tbb(n)    = tt(i0f+n)
               tbb(nh+n) = tbb(n)
            enddo
            write(l6,'(/i3,a,1x,a)') ib,'. Block ',tb(ib)
            nx = nhb
            nx = 5
            call aprint(2,3,6,nhb,nx,nx,hb,tbb,tbb,'HB')
            read*
         endif
         call sdiag(nhb,nhb,hb,e,hb,ez,+1)
c
c------- store eigenvalues and wave functions
c------- particles
         ka(ib,it) = klp
         do k = 1,nf
	    klp = klp + 1
            equ(klp,it) = e(nh+k)
            do n = 1,nhb
               fguv(n,klp,it) = hb(n+(nh+k-1)*nhb)
            enddo
            v2 = zero
            do n = 1,nh
               v2 = v2 + fguv(nh+n,klp,it)**2
            enddo
            if (v2.lt.zero) v2 = zero
            if (v2.gt.one)  v2 = one
            sn = sn + v2*mul
         enddo
	 kd(ib,it) = klp - ka(ib,it)
c
c------- anti-particles
	 ka(ib,it+2) = kla
         do k = 1,ng
	    kla = kla + 1
            equ(kla,it+2) = e(ng-k+1)
            do n = 1,nhb
               fguv(n,kla,it+2) = hb(n+(ng-k)*nhb)
            enddo
            v2 = zero
            do n = 1,nh                   ! no-sea approximation
               v2 = v2 + fguv(nh+n,kla,it+2)**2
            enddo                         ! no-sea approximation
            sn = sn + v2*mul
         enddo
	 kd(ib,it+2) = kla - ka(ib,it+2)
c
c======================================================================c
      enddo   ! ib
c======================================================================c
      if (lit.gt.1) dd = (sn - snold)/(al - alold)
c------- calculation of a new lambda-value
      alold = al
      dn    = sn - tz(it)
      if (dn.lt.zero) then
          xl = al
      else
          xh = al
      endif
      if (lit.eq.1) then
         if(dabs(dn).le.0.1d0) then
            al = al - dn
         else
             al = al - 0.1d0*sign(one,dn)
         endif
      else
c           secant method
         if (dd.eq.zero) dd = 1.d-20
         al    = al - dn/dd
         if (al.lt.xl.or.al.gt.xh) then
c              bisection
            al = half*(xl+xh)
            bbb = 'B'
         endif
      endif
      if (abs(al-alold).lt.epsl) goto 30
c
      if (lprl.or.lit.gt.10) then
          write(l6,113) lit,'. L-Iteration: ',bbb,alold,dn,al
          write(6,113)  lit,'. L-Iteration: ',bbb,alold,dn,al
          bbb = ' '
      endif
c
c---- end of lambda-loop
      enddo
      write(l6,'(a,i4,a)')
     &     ' Lambda-Iteration interupted after',lit-1,' steps'
      stop
   30 if (lprl) then
         write(l6,101) lit,'. Lambda-Iteration successful:',it,al,dn,sn
c        write(6,101) lit,'. Lambda-Iteration successful:',it,al,dn,sn
      endif
      ala(it) = al

      if (lpr) then
      write(l6,*) ' ****** END DIRHB **********************************'
      endif
c
  101 format(i4,a,i4,3f13.8)
  113 format(i4,a,a1,3f13.8)
c
      return
C-end-DIRHB
      end

C=======================================================================

      subroutine nucleus(is,npro,te)

C=======================================================================
C
C     is = 1 determines the symbol for a given proton number npro
c          2 determines the proton number for a given symbol te
c
C-----------------------------------------------------------------------
C
      PARAMETER (MAXZ=140)
C
      CHARACTER TE*2,T*(2*MAXZ+2)
C
      T(  1: 40) = '  _HHeLiBe_B_C_N_O_FNeNaMgAlSi_P_SClAr_K'
      T( 41: 80) = 'CaSsTi_VCrMnFeCoNiCuZnGaGeAsSeBrKrRbSr_Y'
      T( 81:120) = 'ZrNbMoTcRuRhPdAgCdInSnSbTe_IXeCsBaLaCePr'
      T(121:160) = 'NdPmSmEuGdTbDyHoErTmYbLuHfTa_WReOsIrPtAu'
      T(161:200) = 'HgTlPbBiPoAtRnFrRaAcThPa_UNpPuAmCmBkCfEs'
      T(201:240) = 'FmMdNoLrRfHaSgNsHsMrDsRt1213141516171819'
      T(241:280) = '2021222324252627282930313233343536373839'
      T(281:282) = '40'
c
c ... Rf is called also as Ku (kurchatovium)
c ... Ha: IUPAC calls it as dubnium (Db). J.Chem.Educ. 1997, 74, 1258
c ... Ha is called also as Db (Dubnium)
c
      if (is.eq.1) then
         if (npro.lt.0.or.npro.gt.maxz) stop 'in NUCLEUS: npro wrong'
         te = t(2*npro+1:2*npro+2)
         return
      else
c
         do np = 0,maxz
            if (te.eq.t(2*np+1:2*np+2)) then
               npro = np
	       if (npro.gt.maxz) write(6,100) TE
               return
            endif
         enddo
c
         write(6,100) TE
  100    format(//,' NUCLEUS ',A2,'  UNKNOWN')
      endif
c
      stop
C-END-NUCLEUS
      END
C=======================================================================

      subroutine gfv(ngfv)

C=======================================================================
C
C     Calculates sign, sqrt, factorials, etc. of integers and half int.
C
c     iv(n)  =  (-1)**n
c     sq(n)  =  sqrt(n)
c     sqi(n) =  1/sqrt(n)
c     sqh(n) =  sqrt(n+1/2)
c     shi(n) =  1/sqrt(n+1/2)
c     fak(n) =  n!
c     fad(n) =  (2*n+1)!!
c     fdi(n) =  1/(2*n+1)!!
c     fi(n)  =  1/n!
c     wf(n)  =  sqrt(n!)
c     wfi(n) =  1/sqrt(n!)
c     wfd(n) =  sqrt((2*n+1)!!)
c     gm2(n) =  gamma(n+1/2)
c     gmi(n) =  1/gamma(n+1/2)
c     wg(n)  =  sqrt(gamma(n+1/2))
c     wgi(n) =  1/sqrt(gamma(n+1/2))
C
C-----------------------------------------------------------------------
      use dirqfampar;
      use gfvmod;

      implicit double precision (a-h,o-z)
c
!ab   parameter (IGFV = 100)
c
!ab   common /gfviv / iv(-IGFV:IGFV)
!ab   common /gfvsq / sq(0:IGFV)
!ab   common /gfvsqi/ sqi(0:IGFV)
!ab   common /gfvsqh/ sqh(0:IGFV)
!ab   common /gfvshi/ shi(0:IGFV)
!ab   common /gfvfak/ fak(0:IGFV)
!ab   common /gfvfad/ fad(0:IGFV)
!ab   common /gfvfi / fi(0:IGFV)
!ab   common /gfvfdi/ fdi(0:IGFV)
!ab   common /gfvwf / wf(0:IGFV)
!ab   common /gfvwfi/ wfi(0:IGFV)
!ab   common /gfvwfd/ wfd(0:IGFV)
!ab   common /gfvgm2/ gm2(0:IGFV)
!ab   common /gfvgmi/ gmi(0:IGFV)
!ab   common /gfvwg / wg(0:IGFV)
!ab   common /gfvwgi/ wgi(0:IGFV)
      common /mathco/ zero,one,two,half,third,pi
c
c---- mathemathical constants
c     data zero/0.0d0/,one/1.d0/,two/2.d0/
c     data half/0.5d0/,third/0.333333333333333333d0/
c     data pi/3.141592653589793d0/
c
      zero  = 0.d0
      one   = 1.d0
      two   = 2.d0
      half  = one/two
      third = one/3.d0
      pi    = 4*atan(one)
c
      iv(0)  = +1
      sq(0)  =  zero
      sqi(0) =  1.d30
      sqh(0) =  sqrt(half)
      shi(0) =  1/sqh(0)
      fak(0) =  one
      fad(0) =  one
      fi(0)  =  one
      fdi(0) =  one
      wf(0)  =  one
      wfi(0) =  one
      wfd(0)=  one
c     gm2(0) = Gamma(1/2) = sqrt(pi)
      gm2(0) =  sqrt(pi)
      gmi(0) =  1/gm2(0)
      wg(0)  =  sqrt(gm2(0))
      wgi(0) =  1/wg(0)
      do i = 1,IGFV
         iv(i)  = -iv(i-1)
         iv(-i) = iv(i)
         sq(i)  = dsqrt(dfloat(i))
         sqi(i) = one/sq(i)
         sqh(i) = sqrt(i+half)
         shi(i) = one/sqh(i)
         fak(i) = i*fak(i-1)
	 fad(i) = (2*i+1)*fad(i-1)
         fi(i)  = one/fak(i)
         fdi(i) = one/fad(i)
         wf(i)  = sq(i)*wf(i-1)
         wfi(i) = one/wf(i)
	 wfd(i) = sqrt(fad(i))
         gm2(i) = (i-half)*gm2(i-1)
         gmi(i) = one/gm2(i)
         wg(i)  = sqh(i-1)*wg(i-1)
         wgi(i) = one/wg(i)
      enddo
c
      ngfv = IGFV
c
c     write(6,*) ' ****** END GFV *************************************'
      return
c-end-GFV
      end
c======================================================================c
c
      subroutine binom(ngfv)
c
c======================================================================c
C     THE ARRAY OF BINOMIAL COEFFICIENTS
C     BIN(I,J)= = I!/J!/(I-J)!
c----------------------------------------------------------------------c
      use dirqfampar;
      use gfvmod;

      implicit double precision (a-h,o-z)
c
!ab   parameter (IGFV = 100)
c
!ab   common /bin0/ bin(0:IGFV,0:IGFV)
c
      do i = 0,IGFV
         do k = 0,IGFV
            bin(i,k) = 0.d0
         enddo   ! k
         bin(i,0)=1.d0
         bin(i,i)=1.d0
         do k = 1,i/2
            bin(i,k)   = dnint(bin(i,k-1)/dfloat(k)*dfloat(i-k+1))
            bin(i,i-k) = bin(i,k)
         enddo
      enddo
      ngfv = IGFV
      return
c-end-BINOM
      end
c======================================================================c

      function itestc()

c======================================================================c
c
C    yields 1, if interrupted
c           2, if convergence
c
c    the iteration is determined by the parameter XMIX
c    it can be fixed, or automatically adjusted according to
c
c     IAUT = 0 fixed value for XMIX
c          = 1 automatic adjustment of XMIX
c
c     INXT = 0 with out inputs from the console
c          > 0 after INXT iterations INX and XMIX are read
c
c     INX  = 0 immediate interruption of the iteration
c          > 0 further INX steps with fixed XMIX, which is read
c          < 0 further ABS(INX) steps with automatic change of XMIX
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
      integer itestc
c
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
c
      if (ii.ge.2.and.si.lt.epsi) then
         itestc = 2
         return
      endif
c
c     if you want to change xmix at the console, please remove the
c     next line !
c     inxt = maxi
c     change of XMIX by reading from the console:
c      if (inxt.eq.ii) then
c         write(6,*)
c     &   'next stop? (0 right now, >0 fixed xmix, <0 autom. xmix)'
c         read(*,*)  inx
c        for running the program without console
        inx=-30
c
c         if (inx.eq.0) then
c            itestc = 1
c            return
c         endif
c         if (inx.lt.0) then
c            iaut = 1
c         endif
c         if (inx.gt.0) then
c            iaut = 0
c         endif
         inxt = ii+iabs(inx)
c         write(6,*) 'new value for xmix?'
c         read(*,*) xmix
c        write(6,*) inxt,xmix
c         xmix0 = xmix
c      endif
c
c     automatic change of XMIX:
c      if ((si.lt.siold).and.iaut.eq.1) THEN
c         xmix = xmix * 1.04
c         if (xmix.gt.xmax) xmix = xmax
c      else
c         xmix = xmix0
c      endif
      siold  = si
      itestc = 0
c
      return
c-end-ITESTC
      end
c======================================================================c

       subroutine gaush(fak,lpr)

c======================================================================c
c
c     Gauss-Hermite integration data
c     ------------------------------
c     for integration from minus to plus infinity
c                     or from 0 to infinity
c
c     ph  =  wh * exp(-xh**2)
c
c     whh=ph
c
c     \int_-\infty^+\infty  f(z) exp(-z**2) dz  =   \sum_i f(xh(i)) ph(i)
c     possible alternative
c     \int_-\infty^+\infty  f(z) dz             =   \sum_i f(xh(i)) wh(i)
c
c----------------------------------------------------------------------
      use dirqfampar;
      use gaushmod;
      use gaussh;

      implicit double precision (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
!ab   dimension x(2*NGH),w(2*NGH)
c
!ab   common /gaussh/ xh(0:NGH),wh(0:NGH),xb(0:NGH)
      common /basnnn/ n0f,n0b
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      if (lpr) then
      write(l6,*) '****** BEGIN GAUSH ********************************'
      endif
c
c
      CALL gauher(x,w,NGH2)
c
      if (lpr) write(l6,100) 'Gauss Hermit integration',NGH
      do i = 1,NGH
         xh(i) = x(NGH+1-i)
         ph    = fak*w(NGH+1-i)
         wh(i) = ph*dexp(xh(i)*xh(i))
         if (lpr) write(6,102) i,xh(i),ph,wh(i)
      enddo
      xh(0)=1.d-10
      wh(0)=1.d-10
c
  100 format('  GAUSH:  ',a,' NGH =',i3,/,2x,52(1h-))
  102 format(i3,3d20.11)
c
      if (lpr) then
      write(l6,*) '****** END GAUSH **********************************'
      endif
c
      return
c-end-GAUSH
      end
c======================================================================c
c
      SUBROUTINE gauher(x,w,n)
c
c======================================================================c
      INTEGER n,MAXIT
      DOUBLE PRECISION w(n),x(n)
      DOUBLE PRECISION EPS,PIM4
      PARAMETER (EPS=3.D-14,PIM4=.7511255444649425D0,MAXIT=10)
      INTEGER i,its,j,m
      DOUBLE PRECISION p1,p2,p3,pp,z,z1
      m=(n+1)/2
      !ab
      z = 0.d0;
      !ab
      do 13 i=1,m
        if(i.eq.1)then
          z=sqrt(float(2*n+1))-1.85575*(2*n+1)**(-.16667)
        else if(i.eq.2)then
          z=z-1.14*n**.426/z
        else if (i.eq.3)then
          z=1.86*z-.86*x(1)
        else if (i.eq.4)then
          z=1.91*z-.91*x(2)
        else
          z=2.*z-x(i-2)
        endif
        do 12 its=1,MAXIT
          p1=PIM4
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=z*sqrt(2.d0/j)*p2-sqrt(dble(j-1)/dble(j))*p3
11        continue
          pp=sqrt(2.d0*n)*p2
          z1=z
          z=z1-p1/pp
          if(abs(z-z1).le.EPS)goto 1
12      continue
        stop 'too many iterations in gauher'
1       x(i)=z
        x(n+1-i)=-z
        w(i)=2.d0/(pp*pp)
        w(n+1-i)=w(i)
13    continue
c
      return
c-end-GAUHER
      END
c======================================================================c

       subroutine gausl(lpr)

c======================================================================c
c
c     Gauss-Laguerre integration data
c     -------------------------------
c     ph  =  wh * exp(-xh**2)
c     pl  =  wl * exp(-xl)
c
c     \int_0^\infty f(x) exp(-x) de     =   \sum_i f(xl(i)) pl(i)
c     \int_0^\infty f(x) dx             =   \sum_i f(xl(i)) wl(i)
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use gaussl;

      implicit double precision (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
!ab   common /gaussl/ xl(0:ngl),wl(0:ngl),sxl(0:ngl),rb(0:ngl)
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
c     dimension xll12(12),wll12(12)
c
c     gauss-laguerre -meshpoints for n = 12
c     DATA xll12/.11572211736D+00,  .61175748452D+00,  .15126102698D+01,
c    2           .28337513377D+01,  .45992276394D+01,  .68445254531D+01,
c    3           .96213168425D+01,  .13006054993D+02,  .17116855187D+02,
c    4           .22151090379D+02,  .28487967251D+02,  .37099121044D+02/

c
c     gauss-laguerre -weights    for n = 12
c     DATA wll12/.26473137106D+00,  .37775927587D+00,  .24408201132D+00,
c    2           .90449222212D-01,  .20102381155D-01,  .26639735419D-02,
c    3           .20323159266D-03,  .83650558568D-05,  .16684938765D-06,
c    4           .13423910305D-08,  .30616016350D-11,  .81480774674D-15/
c
c
      lx = l6
      if (lpr) then
      l6 = 6
      write(l6,*) '****** BEGIN GAUSL ********************************'
      endif
c
c     ngl0=ngl
c
c     if (ngl0.eq.12) then
c        do i = 1,ngl
c           xl(i)  = xll12(i)
c           wl(i)  = wll12(i)*exp(xl(i))
c           sxl(i) = sqrt(xl(i))
c        enddo
c     else
c        stop 'in GAUSL: # of Laguerre mesh points not implemented'
c     endif
c
c
c
      CALL gaulag(xl(1),wl(1),ngl,zero)
c
      if (lpr) then
         write(6,100) ngl
         write(l6,101)
         do i = 1,ngl
            write(6,102) i,xl(i),wl(i),wl(i)*exp(xl(i))
         enddo   ! i
      endif
c
      do i = 1,ngl
         wl(i)  = wl(i)*exp(xl(i))
         sxl(i) = sqrt(xl(i))
      enddo   ! i
      xl(0)=1.d-10
      wl(0)=1.d-10
      sxl(0)=sqrt(xl(0))
c
      if (lpr) then
      write(l6,*) '****** END GAUSL **********************************'
      read*
      l6 = lx
      endif
c
  100 format('  GAUSL:  G-L-Integration  ngl =',i3)
  101 format(2x,33(1h-))
  102 format(i3,3d20.11)
c
      return
c-end-GAUSL
      end
c======================================================================c
      DOUBLE PRECISION FUNCTION gammln(xx)
c======================================================================c
c     implicit real*8(a-h,o-z)
      DOUBLE PRECISION xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
c======================================================================c
      SUBROUTINE gaulag(x,w,n,alf)
c======================================================================c
c     implicit real*8(a-h,o-z)
      INTEGER n,MAXIT
      DOUBLE PRECISION alf,w(n),x(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.D-14,MAXIT=10)
CU    USES gammln
      INTEGER i,its,j
      DOUBLE PRECISION ai,gammln
      DOUBLE PRECISION p1,p2,p3,pp,z,z1
      do 13 i=1,n
        if(i.eq.1)then
          z=(1.+alf)*(3.+.92*alf)/(1.+2.4*n+1.8*alf)
        else if(i.eq.2)then
          z=z+(15.+6.25*alf)/(1.+.9*alf+2.5*n)
        else
          ai=i-2
          z=z+((1.+2.55*ai)/(1.9*ai)+1.26*ai*alf/(1.+3.5*ai))*
     *(z-x(i-2))/(1.+.3*alf)
        endif
        do 12 its=1,MAXIT
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
            p3=p2
            p2=p1
            p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j
11        continue
          pp=(n*p1-(n+alf)*p2)/z
          z1=z
          z=z1-p1/pp
          if(abs(z-z1).le.EPS)goto 1
12      continue
        stop 'too many iterations in gaulag'
1       x(i)=z
        w(i)=-exp(gammln(alf+n)-gammln(dfloat(n)))/(pp*n*p2)
13    continue
      return
      END
c======================================================================c

      subroutine aprint(is,it,ns,ma,n1,n2,a,t1,t2,text)

c======================================================================c
C
C     IS = 1    Full matrix
C          2    Lower diagonal matrix
c          3    specially stored symmetric matrix
C
C     IT = 1    numbers for rows and columns
C          2    text for rows and numbers for columns
C          3    text for rows and columns
C
C     NS = 1     FORMAT  8F8.4      80 Coulums
C     NS = 2     FORMAT  8f8.2      80 Coulums
C     NS = 3     FORMAT 17F4.1      80 Coulums
C     NS = 4     FORMAT 30F4.1     120 Coulums
C     NS = 5     FORMAT  5F12.8     80 Coulums
C     NS = 6     FORMAT  5F12.4     80 Coulums
C     NS = 7     FORMAT  4E13.6     80 Coulums
C     NS = 8     FORMAT  8E15.8    130 Coulums
c
c----------------------------------------------------------------------c
      implicit double precision (a-h,o-z)
C
      character*8 t1(n1),t2(n2)
      character text*(*)
C
      dimension a(ma*n2)
C
      character*30 fmt1,fmt2
      character*20 fti,ftt,fmt(8),fmti(8),fmtt(8)
      dimension nsp(8)
c
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      data nsp/8,8,17,30,5,5,4,8/
      data fmt /'8f8.4)',            '8F8.2)',
     &          '17f4.1)',           '30f4.1)',
     &          '5f12.8)',           '5f12.4)',
     &          '4e13.6)',           '8e15.8)'/
      data fmti/'(11x,8(i4,4x))',    '(11x,8(i4,4x))',
     &          '(11x,17(1x,i2,1x))','(11x,30(1x,i2,1x))',
     &          '(11x,6(i4,8x))',    '(11x,10(i4,8x))',
     &          '(11x,5(i4,9x))',    '(11x,8(i4,11x))'/
      data fmtt/'(11x,8a8)',         '(11x,8a8)',
     &          '(11x,17a4)',        '(11x,30a4)',
     &          '(11x,5(a8,2x))',    '(11x,5(a8,4x))',
     &          '(11x,4(a8,5x))',    '(11x,8(a8,7x))'/
C
      fmt1   = '(4x,i3,4x,' // fmt(ns)
      fmt2   = '(1x,a8,2x' // fmt(ns)
      fti    = fmti(ns)
      ftt    = fmtt(ns)
      nspalt = nsp(ns)

C
      if (ma.eq.0.or.n1.eq.0.or.n2.eq.0) return
      write(l6,'(//,3x,a)') text
C
      ka = 1
      ke = nspalt
      nteil = n2/nspalt
      if (nteil*nspalt.ne.n2) nteil = nteil + 1
C
      do  10  nt = 1,nteil
      if (n2.gt.nspalt)  write(L6,100)  nt
  100 format(//, 10x,'Part',i5,' of the Matrix',/)
      if(nt.eq.nteil) ke = n2
      if (it.lt.3) then
        write(L6,fti) (k,k=ka,ke)
      else
        write(L6,ftt) (t2(k),k=ka,ke)
      endif
C
      do 20  i=1,n1
         kee=ke
         if (is.ge.2.and.ke.gt.i) kee=i
         if (ka.gt.kee) goto 20
         if (is.eq.3) then
            if (it.eq.1) then
              write(l6,fmt1) i,(a(i+(k-1)*(n1+n1-k)/2),k=ka,kee)
            else
              write(l6,fmt2) t1(i),(a(i+(k-1)*(n1+n1-k)/2),k=ka,kee)
            endif
         else
            if (it.eq.1) then
               write(l6,fmt1) i,(a(i+(k-1)*ma),k=ka,kee)
            else
               write(l6,fmt2) t1(i),(a(i+(k-1)*ma),k=ka,kee)
            endif
         endif
   20 continue
c
      ka=ka+nspalt
      ke=ke+nspalt
   10 continue
C
      return
C-end-APRINT
      end
C=======================================================================
c
      subroutine mzero(ma,m1,m2,aa)
c
C=======================================================================
c
c     AA -> 0
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
C
      dimension aa(ma,m2)
c
      do k = 1,m2
      do i = 1,m1
         aa(i,k) = 0.d0
      enddo   ! i
      enddo   ! k
c
      return
c-end-MZERO
      END
C=======================================================================

      subroutine lingd(ma,mx,n,m,a,x,d,ifl)

C=======================================================================
C
C     solves the system of linear equations A*X = B
C     at the beginning the matrix B is stored in X
C     during the calculation it will be overwritten
C     D is the determinant of A
C
C-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
C
c     solves the system A*X=B, where B is at the beginning on X
c     it will be overwritten lateron, d is the determinant
C
      dimension  a(ma,m),x(mx,m)
C
      data tollim/1.d-10/,one/1.d0/,zero/0.d0/
C
      !ab
      cp = -1;
      !ab

      ifl=1
      p=zero
      do 10 i=1,n
         q=zero
         do 20 j=1,n
 20      q=q+ abs(a(i,j))
         if (q.gt.p)   p=q
 10   continue
      tol=tollim*p
      d=one
      do 30 k=1,n
         p=zero
         do 40 j=k,n
            q = abs(a(j,k))
            if (q.lt.p) goto 40
            p=q
            i=j
 40      continue
         if (p.gt.tol) goto 70
         write (6,200) ('-',j=1,80),tol,i,k,a(i,k),('-',j=1,80)
  200    format (/1x,80a1/' *****  ERROR IN LINGD , TOLERANZ =',e10.4,
     1 ' VALUE OF A(',i3,',',i3,') IS ',e10.4/1x,80a1)
         ifl=-1
         return
   70    cp=one/a(i,k)
         if (i.eq.k) goto 90
         d=-d
         do 81 l=1,m
            cq=x(i,l)
            x(i,l)=x(k,l)
   81       x(k,l)=cq
         do 80 l=k,n
            cq=a(i,l)
            a(i,l)=a(k,l)
   80       a(k,l)=cq
   90       d=d*a(k,k)
            if (k.eq.n) goto 1
            k1=k+1
            do 120 i=k1,n
               cq=a(i,k)*cp
               do 106 l=1,m
  106             x(i,l)=x(i,l)-cq*x(k,l)
               do 120 l=k1,n
  120             a(i,l)=a(i,l)-cq*a(k,l)
   30 continue
    1 do 126 l=1,m
  126    x(n,l)=x(n,l)*cp
         if (n.eq.1) return
         n1=n-1
         do 140 k=1,n1
            cp=one/a(n-k,n-k)
            do 140 l=1,m
               cq=x(n-k,l)
               do 141 i=1,k
  141             cq=cq-a(n-k,n+1-i)*x(n+1-i,l)
  140          x(n-k,l)=cq*cp
c
c
      return
c-end-LINGD
      end
C=======================================================================

      subroutine sdiag(nmax,n,a,d,x,e,is)

C=======================================================================
C
C     A   matrix to be diagonalized
C     D   eigenvalues
C     X   eigenvectors
C     E   auxiliary field
C     IS = 1  eigenvalues are ordered and major component of X is positiv
C          0  eigenvalues are not ordered
C-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
      dimension a(nmax,nmax),x(nmax,nmax),e(n),d(n)
C
      data tol,eps/1.e-32,1.e-10/
C
      if (n.eq.1) then
         d(1)=a(1,1)
         x(1,1)=1.
         return
      endif
c
      do 10 i=1,n
      do 10 j=1,i
   10    x(i,j)=a(i,j)
c
ccc   householder-reduktion
      i=n
   15 if (i-2) 200,20,20
   20 l=i-2
      f=x(i,i-1)
      g=f
      h=0
      if (l) 31,31,32
   32 do 30 k=1,l
   30 h=h+x(i,k)*x(i,k)
   31 s=h+f*f
      if (s-tol) 33,34,34
   33 h=0
      goto 100
   34 if (h) 100,100,40
   40 l=l+1
      g= dsqrt(s)
      if (f.ge.0.) g=-g
      h=s-f*g
      hi=1.d0/h
      x(i,i-1)=f-g
      f=0.0
      if (l) 51,51,52
   52 do 50 j=1,l
      x(j,i)=x(i,j)*hi
      s=0.0
      do 55 k=1,j
   55 s=s+x(j,k)*x(i,k)
      j1=j+1
      if (l-j1) 57,58,58
   58 do 59 k=j1,l
   59 s=s+x(k,j)*x(i,k)
   57 e(j)=s*hi
   50 f=f+s*x(j,i)
   51 f=f*hi*.5d0
c
      if (l) 100,100,62
   62 do 60 j=1,l
      s=x(i,j)
      e(j)=e(j)-f*s
      p=e(j)
      do 65 k=1,j
   65 x(j,k)=x(j,k)-s*e(k)-x(i,k)*p
   60 continue
  100 continue
      d(i)=h
      e(i-1)=g
      i=i-1
      goto 15
c
ccc   Bereitstellen der Transformationmatrix
  200 d(1)=0.0
      e(n)=0.0
      b=0.0
      f=0.0
      do 210 i=1,n
      l=i-1
      if (d(i).eq.0.) goto 221
      if (l) 221,221,222
  222 do 220 j=1,l
      s=0.0
      do 225 k=1,l
  225 s=s+x(i,k)*x(k,j)
      do 226 k=1,l
  226 x(k,j)=x(k,j)-s*x(k,i)
  220 continue
  221 d(i)=x(i,i)
      x(i,i)=1
      if (l) 210,210,232
  232 do 230 j=1,l
      x(i,j)=0.0
  230 x(j,i)=0.0
  210 continue
c
ccc   Diagonalisieren der Tri-Diagonal-Matrix
      DO 300 L=1,N
      h=eps*( abs(d(l))+ abs(e(l)))
      if (h.gt.b) b=h
c
ccc   Test fuer Splitting
      do 310 j=l,n
      if ( abs(e(j)).le.b) goto 320
  310 continue
c
ccc   test fuer konvergenz
  320 if (j.eq.l) goto 300
  340 p=(d(l+1)-d(l))/(2*e(l))
      r= dsqrt(p*p+1.d0)
      pr=p+r
      if (p.lt.0.) pr=p-r
      h=d(l)-e(l)/pr
      do 350 i=l,n
  350 d(i)=d(i)-h
      f=f+h
c
ccc   QR-transformation
      p=d(j)
      c=1.d0
      s=0.0
      i=j
  360 i=i-1
      if (i.lt.l) goto 362
      g=c*e(i)
      h=c*p
      if ( abs(p)- abs(e(i))) 363,364,364
  364 c=e(i)/p
      r= dsqrt(c*c+1.d0)
      e(i+1)=s*p*r
      s=c/r
      c=1.d0/r
      goto 365
  363 c=p/e(i)
      r= dsqrt(c*c+1.d0)
      e(i+1)=s*e(i)*r
      s=1.d0/r
      c=c/r
  365 p=c*d(i)-s*g
      d(i+1)=h+s*(c*g+s*d(i))
      do 368 k=1,n
         h=x(k,i+1)
         x(k,i+1)=x(k,i)*s+h*c
  368    x(k,i)=x(k,i)*c-h*s
      goto 360
  362 e(l)=s*p
      d(l)=c*p
      if ( abs(e(l)).gt.b) goto 340
c
ccc   konvergenz
  300 d(l)=d(l)+f
c
      if (is.eq.0) return
ccc   ordnen der eigenwerte
      do 400 i=1,n
      k=i
      p=d(i)
      j1=i+1
      if (j1-n) 401,401,400
  401 do 410 j=j1,n
      if (d(j).ge.p) goto 410
      k=j
      p=d(j)
  410 continue
      if (k.eq.i) goto 400
!ab 420 if (k.eq.i) goto 400
      d(k)=d(i)
      d(i)=p
      do 425 j=1,n
      p=x(j,i)
      x(j,i)=x(j,k)
  425 x(j,k)=p
  400 continue
c
c     signum
      do k = 1,n
         s = 0.0d0
         !ab
         im = 1;
         !ab
         do i = 1,n
            h = abs(x(i,k))
            if (h.gt.s) then
               s  = h
               im = i
            endif
         enddo   ! i
         if (x(im,k).lt.0.0d0) then
            do i = 1,n
               x(i,k) = - x(i,k)
       	    enddo
         endif
      enddo   ! k
c
      return
c-end-SDIAG
      end
c=====================================================================c

      double precision function trace_abt(ma,mb,n1,n2,aa,bb)

c======================================================================c
c
c     calculaties the trace( A * BT )
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
c
      dimension aa(ma,n2),bb(mb,n2)
c
      s = 0.d0
      do i = 1,n1
      do k = 1,n2
         s = s + aa(i,k)*bb(i,k)
      enddo
      enddo
      trace_abt = s
c
      return
c-end-TRACE_ABT
      end
c=====================================================================c

      double precision function trace(ma,n,aa)

c======================================================================c
c
c     calculaties the trace( A )
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
c
      dimension aa(ma,n)
c
      s = 0.d0
      do i = 1,n
         s = s + aa(i,i)
      enddo
      trace = s
c
      return
c-end-TRACE
      end
c=====================================================================c

      double precision function trabt(ma,mb,n1,n2,aa,bb)

c======================================================================c
c
c     calculaties the trace( A * BT )
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
c
      dimension aa(ma,n2),bb(mb,n2)
c
      s = 0.d0
      do i = 1,n1
      do k = 1,n2
         s = s + aa(i,k)*bb(i,k)
      enddo
      enddo
      trabt = s
c
      return
c-end-TRABT
      end
c=====================================================================c

      subroutine matransa(ma,n,aa,dd,dad,zz)

c======================================================================c
c
c     transforms the matrix AA to DAD =  DDT * AA * AD
c
c     ZZ is a two-dimensional auxiliary field
c     the matrix DAD my be identical with the matrix AA
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
c
      dimension aa(ma,n),dd(ma,n),dad(ma,n),zz(ma,n)
c
      do k = 1,n
      do i = 1,n
         s = 0.d0
         do l = 1,n
            s = s + aa(i,l) * dd(l,k)
         enddo
         zz(i,k) = s
      enddo   ! i
      enddo   ! k
      do k = 1,n
      do i = 1,n
         s = 0.d0
         do l = 1,n
            s = s + dd(l,i)*zz(l,k)
         enddo
         dad(i,k) = s
      enddo   ! i
      enddo   ! k
c
      return
c-end-MATRANSA
      end
c=====================================================================c
      subroutine matransb(ma,n,aa,dd,dad,z)

c======================================================================c
c
c     transforms the matrix AA to DAD =  DDT * AA * DD
c
c     Z is a one-dimensional auxiliary field
c     DAD must be different from AA  !
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
c
      dimension aa(ma,n),dd(ma,n),dad(ma,n),z(n)
c
      do k = 1,n
         do l = 1,n
            s = 0.d0
            do m = 1,n
               s = s + aa(l,m) * dd(m,k)
            enddo   ! m
            z(l) = s
         enddo   ! l
         do i = 1,n
            s = 0.d0
            do l = 1,n
               s = s + dd(l,i)*z(l)
            enddo   ! l
            dad(i,k) = s
         enddo   ! i
      enddo   ! k
c
      return
c-end-MATRANSB
      end

c======================================================================c
      subroutine intpol(n,v,vh)
c
c======================================================================c
c
C     Interpolates potentials between mesh points  VH(i) = V(i+1/2)
c     Lagrange 4-point interpolation
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
c
      dimension v(n),vh(n-1)
c
      data c1/0.0625d0/
c
      ie = n-2
      do i = 2,ie
         vh(i)  = (9*(v(i)+v(i+1)) - v(i-1) - v(i+2) )*c1
      enddo
      vh(1)   = (5*(v(1)-v(3)  ) + 15*v(2)   + v(4)  )*c1
      vh(n-1) = (5*(v(n)-v(n-2)) + 15*v(n-1) + v(n-3))*c1
      vh(n)   = 0.d0
c
      return
c-end-INTPOL
      end
c=======================================================================

      subroutine ordi(n,e,mu)

c=======================================================================
c
C     orders a set of numbers according to their size
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
      dimension e(n),mu(n)
c
      do 10 i = 1,n
         k  = i
         p  = e(i)
         if (i.lt.n) then
            do 20 j = i+1,n
               if (e(j).lt.p) then
                  k = j
                  p = e(j)
               endif
   20       continue
            if (k.ne.i) then
               e(k)  = e(i)
               e(i)  = p
               mk    = mu(k)
               mu(k) = mu(i)
               mu(i) = mk
            endif
         endif
   10 continue
c
      return
c-end-ORDI
      end
c=======================================================================

      subroutine ordx(n,e,a1,a2,bb)

c=======================================================================
c
C     orders a set of numbers according to their size
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
      dimension e(n),a1(n),a2(n),bb(n,n)
c
      do 10 i = 1,n
         k  = i
         p  = e(i)
         if (i.lt.n) then
            do 20 j = i+1,n
               if (e(j).lt.p) then
                  k = j
                  p = e(j)
               endif
   20       continue
            if (k.ne.i) then
               e(k)  = e(i)
               e(i)  = p
               x     = a1(k)
               a1(k) = a1(i)
               a1(i) = x
               x     = a2(k)
               a2(k) = a2(i)
               a2(i) = x
               do j = 1,n
                  x       = bb(j,k)
                  bb(j,k) = bb(j,i)
                  bb(j,i) = x
               enddo
            endif
         endif
   10 continue
c
      return
c-end-ORD3
      end
c======================================================================c
c
      subroutine simps(f,n,h,result)
c
c======================================================================c
c
C     This routine performs Simpson's rule integration of a
C     function defined by a table of equispaced values.
C     Parameters are:
C       F  ---   Array of values of the function
C       N  ---   Number of points
C       H  ---   The uniform spacing between X values
C       RESULT  ---   Estimate of the integral that is returned to caller.
c
c----------------------------------------------------------------------c
c
      implicit real*8 (a-h,o-z)
C
      dimension f(n)
c
      common /mathco/ zero,one,two,half,third,pi
c
      data c3d8/0.375d0/
c
c
C     Check to see if number of panels is even.  Number of panels is N-1.
      npanel = n-1
      nhalf  = npanel/2
      nbegin = 1
      result = zero
      if ((npanel-2*nhalf).ne.0) then
C
C     Number of panels is odd.  Use 3/8 rule on first three
C     panels, 1/3 rule on rest of them.
C
         result = h*c3d8*( f(1) + 3*(f(2)+f(3)) + f(4) )
         if (n.eq.4) then
            return
         endif
         nbegin=4
      endif
C
C Apply 1/3 rule - add in first, second, last values
C
      result = result + h*third*( f(nbegin) + 4*f(nbegin+1) + f(n) )
      nbegin = nbegin+2
      if (nbegin.eq.n) then
            return
         else
            x = zero
            nend = n-2
            do i = nbegin,nend,2
               x = x + f(i) + 2*f(i+1)
            enddo
            result = result + h*two*third*x
c
            return
      endif
c
c-end-SIMPS
      end
c=======================================================================
      integer function kap(l,j)
c======================================================================c
c
c     kappa = (2j+1)(l-j)
      kap = j * (2*(l-j) + 1)
c
      return
c-end-KAP
      end
c======================================================================c
      integer function lfkap(kappa)
c======================================================================c
      if (kappa.gt.0) then
         lfkap = kappa
      else
         lfkap = - kappa - 1
      endif
      return
c-end-LFKAP
      end
c======================================================================c
      integer function lgkap(kappa)
c======================================================================c
      if (kappa.gt.0) then
         lgkap = kappa - 1
      else
         lgkap = - kappa
      endif
      return
c-end-LGKAP
      end
c======================================================================c
      integer function lfgkap(kappa,is)
c======================================================================c
      if (is.eq.1) then
         if (kappa.gt.0) then
            lfgkap = kappa
         else
            lfgkap = - kappa - 1
         endif
      else
         if (kappa.gt.0) then
            lfgkap = kappa - 1
         else
            lfgkap = - kappa
         endif
      endif
      return
c-end-LFGKAP
      end
c=======================================================================

      subroutine readms(ltape,n,aa)

c=======================================================================
c
C     reads from tape
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
      dimension aa(n)
c
      read(ltape) aa
c
      return
c-end-READMS
      end
c=======================================================================

      subroutine writms(ltape,n,aa)

c=======================================================================
c
C     writes on tape
c
c-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
C
      dimension aa(n)
c
      write(ltape) aa
c
      return
c-end-WRITMS
      end
c===========================================================================
      SUBROUTINE check_abs_le(x1,x2,text)
c===========================================================================
      REAL*8   x1,x2
      CHARACTER TEXT*(*)
      IF (abs(x1).gt.x2) THEN
         WRITE(6,*) 'STOP in ',text,x1,x2
         STOP
      ENDIF
      RETURN
      END   ! SUBROUTINE check_abs_le
c===========================================================================
      SUBROUTINE icheck_n1len2(n1,n2,text)
c===========================================================================
      INTEGER   n1,n2
      CHARACTER TEXT*(*)
      IF (n1.GT.n2) THEN
         WRITE(6,*) 'STOP in ',text,n1,n2
         STOP
      ENDIF
      RETURN
      END   ! SUBROUTINE icheck_n1len2
c===========================================================================
      SUBROUTINE icheck_n1nen2(n1,n2,text)
c===========================================================================
      INTEGER   n1,n2
      CHARACTER TEXT*(*)
      IF (n1.NE.n2) THEN
         WRITE(6,*) 'STOP in ',text,n1,n2
         STOP
      ENDIF
      RETURN
      END   ! SUBROUTINE icheck_n1len2
c===========================================================================
      subroutine maxa(n,a,imax,smax)
c===========================================================================
      implicit double precision (a-h,o-z)
c
      dimension a(n)
c
      smax = 0.d0
      do i = 1,n
         s = abs(a(i))
         if (s.gt.smax) then
            smax = s
            imax = i
         endif
      enddo   ! i
c     write(6,*) 'max',imax,s
c     read*
c
      return
c-end-MAXA
      end

c======================================================================c
c
      subroutine expect(lpr)
c
c======================================================================c
c
c     calculates expectation values
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use gaussh;
      use gaussl;
      use rhorho;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
      character nucnam*2                                        ! nucnuc
c
      dimension ekt(3),epart(3),ept(3)
      dimension xn(3),xs(3),r2(3),qq(3),hh(3)
      dimension bet2(3),bet4(3)
      dimension emes(4)
c
      common /baspar/ hom,hb0,b0
      common /centma/ cmas(3)
      common /erwar / ea,rms,betg,gamg
      common /con_b2/ betac,q0c,cquad,c0,alaq,calcq0,icstr
      common /con_b2a/ fac0
      common /fermi / ala(2),tz(2)
!ab   common /gaussh/ xh(0:NGH),wh(0:NGH),zb(0:NGH)
!ab   common /gaussl/ xl(0:NGL),wl(0:NGL),sxl(0:NGL),rb(0:NGL)
      common /masses/ amu,ames(4)
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,npr(3),nucnam
      common /pair  / del(2),spk(2),spk0(2)
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
!ab   common /rhorho/ rs(MG,2),rv(MG,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      ihl(ih,il) = 1+ih + il*(NGH+1)
c
      if (lpr) then
      write(l6,*) ' ****** BEGIN EXPECT *****************'
      endif
c
c
c
c======================================================================c
c---- particle number, radii, deformations
c======================================================================c
      do it = 1,3
         xn(it) = zero
         xs(it) = zero
         r2(it) = zero
         qq(it) = zero
         hh(it) = zero
      enddo
c
      do ih = 1,NGH
         z  = zb(ih)
         zz = z**2
      do il = 1,NGL
         rrp = rb(il)**2
c
c------- root mean square radius
         rr = zz + rrp
c
c------- quadrupole moment
         rq = (3*zz - rr)
c
c------- hexadecupole moment
         rh = (8*zz**2 - 24*zz*rrp + 3*rrp**2)
c
         i = ihl(ih,il)
         do it = 1,2
            xn(it) = xn(it) + rv(i,it)
            xs(it) = xs(it) + rs(i,it)
            r2(it) = r2(it) + rv(i,it)*rr
            qq(it) = qq(it) + rv(i,it)*rq
            hh(it) = hh(it) + rv(i,it)*rh
         enddo   ! it
      enddo   ! il
      enddo   ! ih
c
      do it = 1,2
         r2(it) = sqrt(r2(it)/xn(it))
      enddo
      r2(3) = sqrt((npr(1)*r2(1)**2+npr(2)*r2(2)**2)/amas)
      rc    = sqrt(r2(2)**2 + 0.64)
      rms   = r2(3)
c
      xn(3) = xn(1) + xn(2)
      xs(3) = xs(1) + xs(2)
      qq(3) = qq(1) + qq(2)
      calcq0= qq(3)
      hh(3) = hh(1) + hh(2)
      r00   = r0*amas**third
      fac2  = sqrt(5/(4*pi))/2.d0
      fac4  = sqrt(9/(4*pi))/8.d0
      do it = 1,3
         amr2 = 3/(4*pi)*xn(it)*r00**2
         amr4 = 3/(4*pi)*xn(it)*r00**4
         bet2(it) = fac2*qq(it)/amr2
         bet4(it) = fac4*hh(it)/amr4
      enddo
c
      betg  = bet2(3)
      gamg  = zero
      ecstr = fac0*qq(3)

c======================================================================c
c---- single particle energies, kinetic energies and pairing energies
c======================================================================c
      do it = 1,2
c
c
c------- kinetic energy
         ekt(it) = ekin(it)
c
c------- particle energy
         epart(it) = epar(it)
c
c------- pairing energy
         ept(it) = epair(it)

      enddo   ! it
      ekt(3)   = ekt(1) + ekt(2)
      ept(3)   = ept(1) + ept(2)
      epart(3) = epart(1) + epart(2)

c======================================================================c
c---- field energies
c======================================================================c
c
      call efield(emes,er,ecou)
      emest = emes(1) + emes(2) + emes(3) + emes(4)
c
c
c
c======================================================================c
c---- center off mass correction ecm = <P**2>/2Am
c======================================================================c
c
c     harmonic oscillator
      ecm0  = -0.75d0*hom
c
c     Reinhard
c     ecm0 = -17.2d0/(amas**(1.d0/5.d0))
c
      if (icm.lt.2) then
         cmas(1) = ecm0*npr(1)/amas
         cmas(2) = ecm0*npr(2)/amas
         cmas(3) = cmas(1) + cmas(2)
      endif
c
c======================================================================c
c---- Total energy
c======================================================================c
      etot0 = ekt(3) + emest + ecou  + ept(3)
      etot  = etot0 + cmas(3)
      etot1 = epart(3) - ecstr - emest - ecou - er + ept(3)
      etest = etot1 - etot0
      ea    = etot/amas
c
c======================================================================c
c---- printout
c======================================================================c
      if (lpr) then
c
      write(l6,'(/,28x,a,8x,a,9x,a)') 'neutron','proton','total'
c
c     particle number
      write(l6,200) ' particle number .....',xn
      write(l6,200) ' trace scalar density ',xs
c
c     Lambda
      write(l6,200) ' lambda ..............',ala
c
c     trace of kappa
      write(l6,200) ' spk .................',spk
c
c     rms-Radius
      write(l6,200) ' rms-Radius ..........',r2
c
c     charge-Radius
      write(l6,201) ' charge-Radius, R0....',rc
c
c     quadrupole-deformation
      write(l6,200) ' beta2................',bet2
c
c     hexadecupole-deformation
      write(l6,200) ' beta4 ...............',bet4

c     quadrupole-moment
      write(l6,200) ' quadrupole moment ...',qq
      write(l6,200) '    <Q20> ............',(qq(i)*fac2,i=1,3)
c
c     hexadecupole moment
      write(l6,200) ' hexadecupole moment .',hh
      write(l6,200) '    <Q40> ............',(hh(i)*fac4,i=1,3)
      write(l6,*) ' '
c
c     single-particle energy
      write(l6,200) ' Particle Energy .....',epart
      write(l6,202) ' Selfconsistency Test.',etest
      write(l6,*) ' '
c
c     kinetic energy
      write(l6,200) ' Kinetic Energy ......',ekt
c
c     sigma energy
      write(l6,202) ' E-scsc ..............',emes(1)
c
c     omega energy
      write(l6,202) ' E-scve ..............',emes(2)
c
c     delta-energy
      write(l6,202) ' E-vesc...............',emes(3)
c
c     rho-energy
      write(l6,202) ' E-veve...............',emes(4)
c
c     rearrangement energy
      write(l6,202) ' E-rearrangement .....',er
c
c     Coulomb energy
      write(l6,202) ' Coulomb direct ......',ecou
c
c     Constrained energy
      write(l6,202) ' Constrained energy ..',ecstr
c
c     pairing energy
      write(l6,200) ' Pairing Energy ......',ept
c
c     total energy without center of mass correction
      write(l6,202) ' Sum without E-cm ....',etot0
c
c
c     center of mass correction
      if (icm.lt.2) then
      write(l6,200) ' E-cm  3/4*hom .......',cmas
      !ab
      elseif( icm.eq.2 ) then
      !ab(ii is not even defined...)
!ab   elseif ( icm.eq.2 .and. ii.eq.0 ) then
      write(l6,200) ' E-cm <P**2>/2M ......',cmas
      endif
c
c     total energy
      write(l6,202) ' Total Energy ........',etot
c
c     energy per particle
      ea = etot/amas
      write(l6,202) ' E/A .................',ea
c
      endif   ! printout
c
      if (lpr) then
      write(l6,*) '****** END EXPECT *********************************'
      endif
c
  200 format(a,3f15.6)
  201 format(a,15x,2f15.6)
  202 format(a,30x,f15.6)
c
c---------------------------
  100 format(a,7x,3f15.6)
  101 format(a,37x,f15.6)
  102 format(a,35x,f15.6)
  103 format(a,i4,a)
c
      return
c-end-EXPECT
      end

c======================================================================c
c
      real*8 function ekin(it)
c
c======================================================================c
c
c     calculates the kinetic energy
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use ekinmod;
      use blodir;
      use blokap;
      use bloosc;
      use rokaos;
      use single;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
!ab   character tb*6                                            ! blokap
c
!ab   dimension h0(NHHX)
c
      common /baspar/ hom,hb0,b0
!ab   common /blodir/ ka(NBX,4),kd(NBX,4)
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
      common /physco/ hbc,alphi,r0
      common /masses/ amu,ames(4)
      common /mathco/ zero,one,two,half,third,pi
!ab   common /rokaos/ rosh(NHHX,NB2X),aka(MVX,2)
!ab   common /single/ sp(NFGX,NBX)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      ek = zero
      do ib = 1,nb
         nf  = id(ib,1)
         ng  = id(ib,2)
         nh  = nf + ng
         m   = ib + (it-1)*NBX
c
c------- construction of the free Dirac-operator H0
         emcc2 = 2*amu*hbc
         f = hbc/b0
         do n2 = 1,nf
            do n1 = 1,nf
               h0(n1+(n2-1)*nh) = zero
            enddo
            do n1 = 1,ng
               h0(nf+n1+(n2-1)*nh) = f*sp(n1+(n2-1)*ng,ib)
               h0(n2+(nf+n1-1)*nh) = h0(nf+n1+(n2-1)*nh)
            enddo
         enddo
         do n2 = nf+1,nh
            do n1 = nf+1,nh
               h0(n1+(n2-1)*nh) = zero
            enddo
            h0(n2+(n2-1)*nh) = - emcc2
         enddo
         ek = ek + trabt(nh,nh,nh,nh,h0,rosh(1,m))
      enddo   ! ib
      ekin = ek
c
      return
c-end-EKIN
      end
c======================================================================c
c
      real*8 function epar(it)
c
c======================================================================c
c
c     calculates the particle energy
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use eparmod;
      use blodir;
      use blokap;
      use bloosc;
      use gamgam;
      use rokaos;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
!ab   character tb*6                                            ! blokap
c
!ab   dimension ro(NHHX)
c
      common /baspar/ hom,hb0,b0
!ab   common /blodir/ ka(NBX,4),kd(NBX,4)
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
!ab   common /gamgam/ hh(NHHX,NB2X)
      common /masses/ amu,ames(4)
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ hbc,alphi,r0
!ab   common /rokaos/ rosh(NHHX,NB2X),aka(MVX,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      ep = zero
      do ib = 1,nb
         nf  = id(ib,1)
         ng  = id(ib,2)
         nh  = nf + ng
         m   = ib + (it-1)*NBX
         ep = ep + trabt(nh,nh,nh,nh,hh(1,m),rosh(1,m))
      enddo   ! ib
      epar = ep
c
      return
c-end-EPAR
      end
c======================================================================c
c
      real*8 function epair(it)
c
c======================================================================c
c
c     calculates the pairing energy for TMR pairing
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use blokap;
      use bloosc;
      use deldel;
      use rokaos;

      implicit real*8 (a-h,o-z)
!ab   include 'dirqfam.par'
c
!ab   character tb*6                                            ! blokap
c
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
!ab   common /deldel/ de(NHHX,NB2X)
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ hbc,alphi,r0
!ab   common /rokaos/ rosh(NHHX,NB2X),aka(MVX,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /tmrpar/ gl(2),ga
c
      s  = zero
      il = 0
      do ib = 1,nb
         nf = id(ib,1)
         ng = id(ib,2)
         nh = nf + ng
         m  = ib + (it-1)*NBX
         do n2 =  1,nf
         do n1 = n2,nf
            il = il + 1
            s  = s + de(n1+(n2-1)*nh,m)*aka(il,it)
         enddo
         enddo
      enddo   ! ib
      epair=-s/2
c
      return
c-end-EPAIR
      end
c======================================================================c
c
      subroutine efield(emes,er,ecou)
c
c======================================================================c
c
c     calculates  field energies
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use coulmb;
      use couplf;
      use dens;
      use fields;
      use gaucor;
      use rhorho;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      dimension emes(4)
c
!ab   common /coulmb/ cou(MG),drvp(MG)
      common /coupld/ ddmes(4)
!ab   common /couplf/ ff(MG,4,2)
      common /couplg/ ggmes(4),lmes(4)
!ab   common /dens  / ro(MG,4),dro(MG,4)
!ab   common /fields/ phi(MG,4)
!ab   common /gaucor/ ww(MG)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
!ab   common /rhorho/ rs(MG,2),rv(MG,2)
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
c
c======================================================================c
c---- field energies
c======================================================================c
c---- meson-fields
      if (ipc.eq.0) then
         er = zero
         do m = 1,4
            s = zero
            do i = 1,MG
               s  = s  + ggmes(m)*ff(i,m,1)*phi(i,m)*ro(i,m)*ww(i)
               er = er + ggmes(m)*ff(i,m,2)*phi(i,m)*ro(i,m)*ro(i,2)
     &                   *ww(i)
            enddo   ! i
            emes(m) = half*hbc*s
        enddo   ! m
        er = hbc*er        ! rearrangement term
c
c     point coupling
      elseif (ipc.eq.1) then
         er = zero
         do m = 1,4
            s = zero
            do i = 1,MG
               s  = s  + ggmes(m)*ff(i,m,1)*ro(i,m)**2*ww(i)
               er = er + ggmes(m)*ff(i,m,2)*ro(i,m)**2*ro(i,2)*ww(i)
           enddo   ! i
c
c          derivative terms
           do i = 1,MG
              s = s + ddmes(m)*ro(i,m)*dro(i,m)*ww(i)
           enddo   ! i
	       emes(m) = half*hbc*s
        enddo   ! m
        er = half*hbc*er    ! rearrangment term
c
      else
	    stop 'in EFIELD: ipc not properly defined'
      endif   ! ipc
c
c======================================================================c
c---- Coulomb energy
c======================================================================c
      ecou  = zero
      if (icou.ne.0) then
         do i = 1,MG
            ecou  = ecou + cou(i)*rv(i,2)
         enddo   ! i
      endif   ! icou
      ecou  = half*hbc*ecou
c
      return
c-end-EFIELD
      end

c======================================================================c

      subroutine field(lpr)

c======================================================================c
c
c     calculation of the meson-fields in the oscillator basis
c     the fields are given in (fm^-1)
c
c     meson fields:  sig(i) = phi(i,1)*ggsig/gsig
c                    ome(i) = phi(i,2)*ggome/gome
c                    del(i) = phi(i,3)*ggdel/gdel
c                    rho(i) = phi(i,4)*ggrho/grho
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use fieldmod;
      use couplf;
      use dens;
      use fields;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
!ab   logical lpr,lprs
      logical lpr;
c
!ab   dimension gm(4),so(MG),ph(MG)
c
!ab   common /couplf/ ff(MG,4,2)
      common /couplg/ ggmes(4),lmes(4)
      common /couplm/ gmes(4)
!ab   common /dens  / ro(MG,4),dro(MG,4)
!ab   common /fields/ phi(MG,4)
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      if (ipc.eq.1) return

      if (lpr) then
      write(l6,*) '****** BEGIN FIELD *********************************'
      endif

c---- loop over mesons
      do imes = 1,4
         do i = 1,MG
            so(i) = ff(i,imes,1)*ro(i,imes)
	     enddo
         call gordon(imes,so,phi(1,imes))
      enddo   ! imes

      if (lpr) then
         write(l6,102)
  102    format(/,6x,'sigma',10x,'omega',10x,'rho  ')
         do i = 1,MG
            write(6,103) i,phi(i,1)*gm(1),phi(i,2)*gm(2),
     &                     phi(i,3)*gm(3)
         enddo   ! i
      endif
C
      if (lpr) then
      write(l6,*) '****** END FIELD ***********************************'
      endif
cc
c     CALL zeit(2,'FIELD___',6,.false.)
c
  103 format(i3,4f15.10)
c
      return
C-end-FIELD
      end


c======================================================================c

      subroutine gordon(imes,so,phi)

c======================================================================c
C
C     SOLUTION OF THE KLEIN-GORDON-EQU. BY EXPANSION IN OSCILLATOR
c     imes: number of the meson
c     so:   source
c     phi:  meson field
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use gordonmod;
      use bosqua;
      use gaucor;
      use propag;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      dimension so(*), phi(*);
!ab   dimension so(MG),phi(MG)
!ab   dimension rn(NOX),sn(NOX)
c
!ab   common /bosqua/ nzb(NOX),nrb(NOX),NO
!ab   common /gaucor/ ww(MG)
      common /mathco/ zero,one,two,half,third,pi
!ab   common /propag/ gg(NOX,NOX,4),psi(NOX,MG)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
c
      do i = 1,MG
         so(i) = so(i) * ww(i)
      enddo   ! i
      do n = 1,no
         s = zero
         do i = 1,MG
            s = s + psi(n,i) * so(i)
         enddo   ! i
         rn(n) = s
      enddo   ! n
c
      do n1 = 1,no
         s = zero
         do n2 = 1,no
            s = s + gg(n1,n2,imes) * rn(n2)
         enddo   ! n2
         sn(n1) = s
      enddo   ! n1
c
      do i = 1,MG
         s = zero
         do n = 1,no
            s = s + psi(n,i) * sn(n)
         enddo   ! n
         phi(i) = s
      enddo   ! i
c
c
      return
c-end-GORDON
      end
c======================================================================c

      subroutine greemes(lpr)

c======================================================================c
c
C     calculation of the meson and Coulomb-propagator
c
c     calculates the meson-propagators GG
c     DD is the Laplace operator in oscillator space
c     GG = m**2/(-DD+m**2)
c
c
c     imes = 1:   sigma
c            2:   omega
c            3:   rho
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use greemesmod;
      use propag;
      use bosqua;
      use gaussh;
      use gaussl;
      use gfvmod;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
!ab   dimension dd(NOX,NOX),gi(NOX,NOX)
      common /defbas/ beta0,q,bp,bz
      common /baspar/ hom,hb0,b0
!ab   common /propag/ gg(NOX,NOX,4),psi(NOX,MG)
!ab   common /bosqua/ nzb(NOX),nrb(NOX),NO
      common /couplg/ ggmes(4),lmes(4)
!ab   common /gaussh/ xh(0:NGH),wh(0:NGH),zb(0:NGH)
!ab   common /gaussl/ xl(0:NGL),wl(0:NGL),sxl(0:NGL),rb(0:NGL)
!ab   common /gfvsq / sq(0:IGFV)
      common /masses/ amu,ames(4)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp

      if (ipc.eq.1) return
c
      if (lpr) then
      write(l6,*) ' ****** BEGIN GREEMES ******************************'
      endif
c
c
c     calculation of DD
      az2 = one/(b0*bz)**2
      ap2 = one/(b0*bp)**2
      fac = one/(b0**3*bz*bp*bp)
      do i2 = 1,no
         nz2 = nzb(i2)
         nr2 = nrb(i2)
         dd(i2,i2) = -(az2*(nz2+half)+ap2*(nr2+nr2+1))
c
         do i1 = 1,i2-1
            nz1 = nzb(i1)
            nr1 = nrb(i1)
            t   = zero
c
c           z*z
            if (nr1.eq.nr2) then
               if (nz2.eq.nz1+2) t = az2*half*sq(nz1+1)*sq(nz2)
               if (nz2.eq.nz1-2) t = az2*half*sq(nz1)*sq(nz2+1)
            endif
c
c           r*r
            if (nz1.eq.nz2) then
               if (nr2.eq.nr1+1) t = -ap2*nr2
               if (nr2.eq.nr1-1) t = -ap2*nr1
            endif
c
            dd(i1,i2) = t
            dd(i2,i1) = t
         enddo   ! i1
      enddo   ! i2
c
      do imes = 1,4
	 if (lmes(imes).gt.0) then
            f = one/ames(imes)**2
            do n = 1,no
               do k = 1,no
                  gi(n,k)      = -dd(n,k)*f
		  gg(n,k,imes) = zero
               enddo   ! k
               gi(n,n)      = one + gi(n,n)
	           gg(n,n,imes) = fac
	        enddo   ! n
            call lingd(NOX,NOX,NO,NO,gi,gg(1,1,imes),d,ifl)
            if (lpr) then
               nx = 5
	           call aprint(1,1,5,NO,nx,nx,gg(1,1,imes),
     &               ' ',' ','Meson-Propagator')
            endif   ! lpr
	     endif   ! lmes
      enddo   ! imes
c
c
      if (lpr) then
      write(l6,*) ' ****** END GREEMES ********************************'
      endif
c
      return
c-end-GREEMES
      end

c=====================================================================c

      subroutine forces(lpr)

c=====================================================================c
c
c---- options
c     center of mass:
c        icm = 0   hb0 = hbc**2/(2*amu) ecm = 3/4*hom
c              1   hb0 = hb0*(1-1/A)    ecm = 3/4*hom
c              2   hb0 = hb0            ecm = <P**2>/2M
c
c     model-type   DD     density-dependent meson-coupling
c                  PC     point-coupling
c
c---------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
      logical lpr
      character parname*10                     ! common partyp

      common /dforce/ a_s,a_v,a_ts,a_tv,b_s,b_v,b_ts,b_tv,
     &                c_s,c_v,c_ts,c_tv,d_s,d_v,d_ts,d_tv,dsat
      common /masses/ amu,amsig,amome,amdel,amrho
      common /coupld/ ddsig,ddome,dddel,ddrho
      common /couplg/ ggsig,ggome,ggdel,ggrho,lmes(4)
      common /couplm/ gsig,gome,gdel,grho
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /partyp/ parname
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /tmrpar/ gl(2),ga
c
      if (lpr) then
      write(l6,*) ' ****** BEGIN FORCES *******************************'
      endif
c
c=====================================================================c
c---- forces
c=====================================================================c
       if (parname.eq.'DD-ME2') then
         amu    =  939.d0                 ! MeV
         amsig  =  550.12380d0            ! MeV
         amome  =  783.d0                 ! MeV
         amdel  =    1.d-10
	     amrho  =  763.d0                 ! MeV
         gsig   =   10.5396d0
         gome   =   13.0189d0
         gdel   =   zero
	     grho   =    3.6836d0
         b_s    =    1.0943d0
	     c_s    =    1.7057d0
         c_v    =    1.4620d0
         a_tv   =    0.5647d0
         dsat   =    0.152d0
c
         d_s  = one/sqrt(3.d0*c_s)
         a_s  = (one+c_s*(one+d_s)**2)/(one+b_s*(one+d_s)**2)
         d_v  = one/sqrt(3.d0*c_v)
         facs =two*a_s*(b_s-c_s)*(one-3.d0*c_s*(one+d_s)**2)/
     &                          (one+c_s*(1+d_s)**2)**3
         faco = (one-3.d0*c_v*(one+d_v)**2)/(one+c_v*(1+d_v)**2)**3
         x    = facs/(two*faco)
         fac1 = x+c_v*(one+c_v*(one+d_v)**2)
         fac2 = one+c_v*(one+d_v)**2-x*(one+d_v)**2
         b_v = fac1/fac2
         a_v =(one+c_v*(one+d_v)**2)/(one+b_v*(one+d_v)**2)
c
         a_ts   =    zero
         b_ts   =    zero
         c_ts   =    zero
         d_ts   =    zero
c
         ipc    =  0
	     icm    =  2
	     idd    =  2
c===============================================================
       elseif (parname.eq.'DD-PC1') then
c---------------------------------------------------------------
c        G(x) = a + (b + c*x) * exp(-d*x)
c---------------------------------------------------------------

         dsat   =  0.152d0              ! fm^-3
         amu    =  939.d0               ! MeV
c
c        scalar-scalar
         a_s    = -10.0462d0            ! fm^-2
         b_s    =  -9.1504d0            ! fm^-2
         c_s    =  -6.4273d0            ! fm^-2
         d_s    =  +1.3724d0
c
c        scalar-vector
         a_v    =  +5.9195d0            ! fm^-2
         b_v    =  +8.8637d0            ! fm^-2
         c_v    =   0.00000d0           ! fm^-2
         d_v    =  +0.6584d0
c
c        vector-scalar
         a_ts   =    zero
         b_ts   =    zero
         c_ts   =    zero
         d_ts   =    zero
c
c        vector-vector
         a_tv   =   0.0000d0            ! fm^-2
         b_tv   =   1.8360d0            ! fm^-2
         c_tv   =   0.0000d0            ! fm^-2
         d_tv   =   0.6403d0
c
c----- derivative terms                 ! MeV^-4
         ddsig  =  -0.8149d0
         ddome  =   zero
         dddel  =   zero
         ddrho  =   zero
c
c----------------------------------------------------
c
c        gg = 1
         ggsig = one
         ggome = one
         ggdel = one
         ggrho = one
         amsig = zero
         amome = zero
         amdel = zero
         amrho = zero

         icm    = 2
         idd    = 2
         ipc    = 1
      else
c=====================================================================c
          stop 'This type of force is not defined'
      endif
c=====================================================================c
      amu      = amu/hbc
      amsig    = amsig/hbc
      amome    = amome/hbc
      amdel    = amdel/hbc
      amrho    = amrho/hbc

      if (ipc.eq.0) then
         ggsig = -(gsig/(amsig+1.d-10))**2
         ggome = +(gome/(amome+1.d-10))**2
         ggdel = +(gdel/(amdel+1.d-10))**2  !SIGN?!
         ggrho = +(grho/(amrho+1.d-10))**2
         if (abs(gsig).gt.1.d-5) lmes(1) = 1
         if (abs(gome).gt.1.d-5) lmes(2) = 1
         if (abs(gdel).gt.1.d-5) lmes(3) = 1
         if (abs(grho).gt.1.d-5) lmes(4) = 1
      endif

c=====================================================================c
c---- Separable pairing force
c=====================================================================c
c
      gl(1) = -728.d0
      gl(2) = -728.d0
c      ga    = 0.644204936     ! fm
      ga    = 0.415d0
c




      !ab
c=======================================================================c
c=======================================================================c
c=======================================================================c
      intape = 555;
      open( intape , file = 'dirqfam.dat' , status = 'old' );
      read( intape , * );read( intape , * );read( intape , * );
      read( intape , * );read( intape , * );read( intape , * );
      read( intape , * );read( intape , * );read( intape , * );
      read( intape , * );read( intape , * );read( intape , * );
      read( intape , * );read( intape , * );read( intape , * );
      read( intape , * );read( intape , * );read( intape , * );
      read( intape , * );read( intape , * );read( intape , * );
      read( intape , * );
      read( intape , '(18x,i9)' ) i_pairing;

      if( i_pairing .eq. 0 ) then
          gl(1) = 0.D0;
          gl(2) = 0.D0;
      endif

      close(intape);
c=======================================================================c
c=======================================================================c
c=======================================================================c
      !ab





c---- printout of force:
      if (lpr) call pripar
c
      if (lpr) then
      write(l6,*) ' ****** END FORCES *********************************'
      endif
c
      return
c-end-FORCES
      end
c=====================================================================c

      subroutine pripar

c=====================================================================c
c
c     prints parameters of the Lagrangian
c
c---------------------------------------------------------------------c
      use dirqfampar;
c
      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
c
      character parname*10                     ! common partyp
c
      common /partyp/ parname
      common /dforce/ a_m(4),b_m(4),c_m(4),d_m(4),dsat
      common /mathco/ zero,one,two,half,third,pi
      common /masses/ amu,amsig,amome,amdel,amrho
      common /pair  / del(2),spk(2),spk0(2)
      common /coupld/ ddmes(4)
      common /couplg/ ggmes(4),lmes(4)
      common /couplm/ gsig,gome,gdel,grho
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /tmrpar/ gl(2),ga
c
      write(l6,'(/,a,a10)') ' NEDF Parameters: ',parname
c
      if (ipc.eq.0) then
	 write(l6,'(a,f8.3)') ' AMU   = ',amu*hbc
         write(l6,'(a,f8.3,a,f10.4,a,f10.4)') ' msig  = ',amsig*hbc,
     &                 '  gsig = ',gsig,'  Gsig = ',ggmes(1)
         write(l6,'(a,f8.3,a,f10.4,a,f10.4)') ' mome  = ',amome*hbc,
     &                 '  gome = ',gome,'  Gome = ',ggmes(2)
         write(l6,'(a,f8.3,a,f10.4,a,f10.4)') ' mdel  = ',amdel*hbc,
     &                 '  gdel = ',gdel,'  Gdel = ',ggmes(3)
         write(l6,'(a,f8.3,a,f10.4,a,f10.4)') ' mrho  = ',amrho*hbc,
     &                 '  grho = ',grho,'  Grho = ',ggmes(4)
c
      elseif (ipc.eq.1) then
         write(l6,'(11x,a,4x,a,4x,a,4x,a)') 'scasca','scavec',
     &                                         'vecsca','vecvec'
         write(l6,'(a,4f10.6)') ' GG   = ',ggmes
         write(l6,'(a,4f10.6)') ' DD   = ',ddmes
      else
         stop 'in PRIPAR: ipc not properly defined'
      endif   ! ipc=1
      write(l6,*) ' '
      write(l6,'(a)') ' Density dependence parameters:'
      write(l6,'(11x,a,4x,a,4x,a,4x,a)') 'scasca','scavec',
     &                                   'vecsca','vecvec'
      write(l6,'(a,4f10.6)') ' a    = ',(a_m(m),m=1,4)
      write(l6,'(a,4f10.6)') ' b    = ',(b_m(m),m=1,4)
      write(l6,'(a,4f10.6)') ' c    = ',(c_m(m),m=1,4)
      write(l6,'(a,4f10.6)') ' d    = ',(d_m(m),m=1,4)
      write(l6,'(a,f10.6)') ' dsat = ',dsat

      write(l6,*) ' '
      write(l6,'(a)') ' TMR pairing: Tian,Ma,Ring, PRB 676, 44 (2009)'
      write(l6,'(a,f11.5,a)') ' TMR pairing strength    gl =: ',gl(1),
     &              ' [MeV*fm^3]'
      write(l6,'(a,f11.5,a)') ' TMR pairing width        a =: ',
     &                           sqrt(ga),' [fm]'
      return
C-end-PRIPAR
      end
C======================================================================c

      subroutine gamma(it)

c======================================================================c
c
c     calculats the Dirac-Matrix in the Hartee-equation
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use blokap;
      use bloosc;
      use gamgam;
      use potpot;
      use single;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'

!ab   character tb*6                                            ! blokap
c
      common /baspar/ hom,hb0,b0
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
!ab   common /gamgam/ hh(NHHX,NB2X)
      common /masses/ amu,ames(4)
      common /mathco/ zero,one,two,half,third,pi
      common /physco/ hbc,alphi,r0
!ab   common /potpot/ vps(MG,2),vms(MG,2)
!ab   common /single/ sp(NFGX,NBX)
c
      emcc2 = 2*amu*hbc
      f = hbc/b0
c
      do ib = 1,nb
         nf  = id(ib,1)
         ng  = id(ib,2)
         nh  = nf + ng
         i0f = ia(ib,1)
         i0g = ia(ib,2)
         m   = ib + (it-1)*NBX
c
         do n2 = 1,nf
         do n1 = 1,ng
               hh(nf+n1+(n2-1)*nh,m) = f*sp(n1+(n2-1)*ng,ib)
         enddo
         enddo
         call pot(i0f,nh,nf,vps(1,it),hh(1,m))
         call pot(i0g,nh,ng,vms(1,it),hh(nf+1+nf*nh,m))
         do n = nf+1,nh
            hh(n+(n-1)*nh,m) = hh(n+(n-1)*nh,m) - emcc2
         enddo
c
c     symmetrize HH
         do n2 =    1,nh
         do n1 = n2+1,nh
            hh(n2+(n1-1)*nh,m) = hh(n1+(n2-1)*nh,m)
         enddo
         enddo
      enddo  !ib

      return
C-end-GAMMA
      end

c======================================================================c

      subroutine gaupol(lpr)

c======================================================================c
c
c     calculates the wavefunctions for the zylindrical oscillator.
c     the are given as:
c
c     Phi(zz,rr,phi) = 1/((b0**3 bz*bp*bp)^(1/2) * psi_nz(zz) *
c                      psi_(L,nr)(rr) * exp(i*L*phi)/sqrt(2*pi)
c
c     zz is z-coordinate in units of fm
c     rr is perpendicular coordinate in units of fm
c
c     psi_nz(zz)     = N_nz * H_nz(z) * exp(-z*z/2)
c
c     psi_(L,nr)(rr) = N_(nr,L) * sqrt(2) *
c                               eta^(L/2) * L_nr^L(eta) * exp(-eta/2)
c
c     z = zz/(bz*b0),    r = rr/(bp*b0),       eta = r*r
c
c     N_nz     = 1/sqrt(sqrt(pi) * 2^nz * nz!)
c
c     N_(nr,L) = sqrt( nr! / (nr+L)! )
c
c
c     the contribution to the density from the level i is
c
c     rho_k(zz,rr)= 1/(2*pi b0*3 bz*bp*bp) *
c                     ( psi_nz(zz) * psi_(L,nr)(rr) )^2
c
c---- the z-function at meshpoint xh(ih) is stored in QH(nz,ih)
c     such that QH is normalized in such way that the
c     norm integral reads
c
c     \int dz_(-inf)^(inf) (psi_nz)**2 = 1 = \sum_i QH(nz,i)**2
c
c     this means, that QH contains the following factors:
c
c     a)  the z-part of the wavefunction psi_nz(zz)
c     b)  the Gaussian weight sqrt( WH(i) ):
c         \inf_(-inf)^inf f(z) dz = \sum_i f(x_i) * WH(i)
c
c     having QH(nz,i) we get the z-wavefunction:
c
c     psi_nz(zz) =  QH(nz,i) / sqrt( WH(i) )
c
c---- the r-function at meshpoint XL(il) is stored in QL(nr,L,il)
c     such that QL is normalized in such way that the
c     2-dimensional norm integral reads
c
c     \int_0^inf r dr (phi_(L,nr)**2 = 1 = \sum_i QL(nr,L,i)**2
c
c     this means, that QL contains the following factors:
c
c     a)  the part of the wavefunction psi_(nr,L)(rr)
c     b)  a factor sqrt(1/2) from the transformation from r to eta
c     c)  the Gaussian weight sqrt( WL(i) ):
c         \inf_0^inf f(eta) d(eta) = \sum_i f(XL(i)) * WL(i)
c
c     having QL(nr,L,i) we get the r-wavefunction:
c
c     psi_(nr,L)(rr) =  QL(nr,L,i) * sqrt( 2 / WL(i) )
c
c
c---- the density contribution from the level k
c
c     rho_k(zz,rr)= (QH(nz,ih) * QL(nr,L,il))^2 /
c                         (  pi * WH(ih)*WL(il)* b0**3 * bz*bp*pb)
c
c----------------------------------------------------------------------c
c
c     QH1 contains the z-derivatives in the following form:
c
c     d/dz psi_nz(zz) = QH1(nz,i) / sqrt( WH(i) )
c
c     QL1 contains the r-derivatives in the following form:
c
c     d/dr psi_(nr,L)(rr) = QL1(nr,L,i) * sqrt( 2 / WL(i) )
c
c----------------------------------------------------------------------c
c
c     QHB(nz,i) is the z-function for the expansion of the mesonfields
c     QLB(nr,i) is the r-function for the expansion of the mesonfields
c
c     QHB(nz,i) = psi_nz(zz) / (b0*bz)^(1/2)
c     QLB(nr,i) = psi_(nr,L=0)(rr) / ( sqrt(2*pi) * b0*bp)
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use bospol;
      use bosqua;
      use gaussh;
      use gaussl;
      use gfvmod;
      use herpol;
      use lagpol;
      use propag;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
!ab   common /bospol/  qhb(0:N0BX,0:NGH),qlb(0:N0BX,0:NGL)
      common /basnnn/ n0f,n0b
!ab   common /bosqua/ nzb(NOX),nrb(NOX),NO
!ab   common /gaussh/ xh(0:NGH),wh(0:NGH),zb(0:NGH)
!ab   common /gaussl/ xl(0:NGL),wl(0:NGL),sxl(0:NGL),rb(0:NGL)
!ab   common /gfvsq / sq(0:IGFV)
!ab   common /gfvsqi/ sqi(0:IGFV)
!ab   common /gfvwfi/ wfi(0:IGFV)
!ab   common /herpol/ qh(0:NZX,0:NGH),qh1(0:NZX,0:NGH)
!ab   common /lagpol/ ql(0:2*NRX,0:MLX,0:NGL),ql1(0:2*NRX,0:MLX,0:NGL)
!ab   common /propag/ gg(NOX,NOX,4),psi(NOX,MG)
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /zdimos/ nzm,nrm,mlm
c
      if (lpr) then
      write(l6,*) ' ****** BEGIN GAUPOL *******************************'
      endif
c
c
c======================================================================c
c     z - dependence
c======================================================================c
      w4pii = pi**(-0.25d0)
      do ih = 0,NGH
         z   = xh(ih)
         zz  = z*z
         w0  = w4pii*exp(-half*zz)
c
c---------------------------------------------------------
c        the functions qh,qh1 contain already the measure:
         w1  = w0*sqrt(wh(ih))
         wb  = w0
c---------------------------------------------------------
c
c
c------- basis for the fermions
         qh(0,ih)  = w1
         qh(1,ih)  = sq(2)*w1*z
         qh1(0,ih) = -w1*z
         qh1(1,ih) = sq(2)*w1*(one-zz)
         do n = 2,nzm
           qh(n,ih)  = sqi(n)*(sq(2)*z*qh(n-1,ih)-sq(n-1)*qh(n-2,ih))
           qh1(n,ih) = sq(n+n)*qh(n-1,ih) - z*qh(n,ih)
         enddo   ! n
c
c------- basis for bosons
         qhb(0,ih) = wb
         qhb(1,ih) = sq(2)*wb*z
         do n = 2,N0B
           qhb(n,ih) = sqi(n)*(sq(2)*z*qhb(n-1,ih)-sq(n-1)*qhb(n-2,ih))
         enddo   ! n
c
      enddo   ! ih
c
c
c======================================================================c
c     r - dependence
c======================================================================c
      do il = 0,NGL
         x  = xl(il)
	     ri = one/sxl(il)
         w0 = sq(2)*exp(-half*x)
         wb = sq(2)*exp(-half*x) / sqrt(2*pi)
c
c------- basis for the fermions
         do l = 0,mlm
c
c---------------------------------------------------------
c           the functions ql, ql1 contain already the measure
            w1 = w0*sqrt(half*wl(il)*x**l)
c---------------------------------------------------------
c
            ql(0,l,il)  = w1*wfi(l)
            ql(1,l,il)  = w1*(l+1-x)*wfi(l+1)
            ql1(0,l,il) = w1*(l-x)*wfi(l)
            ql1(1,l,il) = w1*(l*l+l-x*(l+l+3)+x*x)*wfi(l+1)
c            do n = 2,nrm
            do n = 2,2*nrm
               dsq = sqi(n)*sqi(n+l)
               d1  = n+n+l-1-x
               d2  = sq(n-1)*sq(n-1+l)
               d3  = n+n+l-x
               d4  = 2*sq(n)*sq(n+l)
               ql(n,l,il)  = (d1*ql(n-1,l,il)-d2*ql(n-2,l,il))*dsq
               ql1(n,l,il) =  d3*ql(n,l,il)  -d4*ql(n-1,l,il)
            enddo
	    do n = 0,nrm
	       ql1(n,l,il) = ql1(n,l,il)*ri
	    enddo
         enddo
c
c------- basis for the bosons
         qlb(0,il) = wb*wfi(0)
         qlb(1,il) = wb*(1-x)*wfi(1)
         do n = 2,N0B
            qlb(n,il) = ((2*n-1-x)*qlb(n-1,il)-(n-1)*qlb(n-2,il))/n
         enddo
      enddo
c
c
c======================================================================c
c     mixed boson wavefunctions
c======================================================================c
c
      do n = 1,NO
         nz = nzb(n)
         nr = nrb(n)
         do il = 0,NGL
         do ih = 0,NGH
	        psi(n,1+ih+il*(NGH+1)) = qhb(nz,ih)*qlb(nr,il)
         enddo   ! NGH
         enddo   ! NGL
      enddo   ! no
c
c
c======================================================================c
c     check for orthogonality
c======================================================================c
      if (lpr) then
c
c        printout
         ix = 3
         do n = 0,nzm
            write(l6,'(a,i2,a)') ' QH(nz=',n,',ih=1...)'
            write(l6,100) ' H     xh',n,(qh(n,ih),ih=1,ix)
            write(l6,100) ' dH/dx xh',n,(qh1(n,ih),ih=1,ix)
c           write(l6,100) ' H  cb*xh',n,(qhb(n,ih)*sqrt(wh(ih)),ih=1,ix)
            write(l6,100) ' H   b   ',n,(qhb(n,ih),ih=1,ix)
  100       format(a,i3,3f15.8)
            write(l6,*) ' '
         enddo
         do l = 0,mlm
            write(l6,*) ' nr ml    QL(nr,l,il=1,...)'
            do n = 0,nrm
               write(l6,'(i4,i3,3f15.8)') n,l,(ql(n,l,il),il=1,ix)
               write(l6,'(i4,i3,3f15.8)') n,l,(ql1(n,l,il),il=1,ix)
               if (l.eq.0)
     &            write(l6,'(i4,i3,3f15.8)') n,l,(qlb(n,il),il=1,ix)
               write(l6,*) ' '
            enddo
         enddo
c
c        test for hermite integration
         do n1 = 0,nzm
         do n2 = 0,n1
            if (mod(n1-n2,2).eq.0) then
               s  = zero
               sb = zero
               s2 = zero
               do ih = 1,NGH
                  s  = s  + qh(n1,ih)*qh(n2,ih)
                  sb = sb + qhb(n1,ih)*qhb(n2,ih)*wh(ih)
               enddo
               write(l6,110) n1,n2,s,sb
  110          format(' G.Hermit: n1 =',i3,'  n2 =',i3,3f15.8)
            endif
         enddo
         enddo
c
c        test for Laguerre integration
         do l = 0,mlm
            do n1 = 0,nrm
            do n2 = 0,n1
               s  = zero
               sb = zero
               s2 = zero
               do il = 1,NGL
                  s  = s + ql(n1,l,il)*ql(n2,l,il)
                  if (l.eq.0)
     &            sb = sb + qlb(n1,il)*qlb(n2,il)*wl(il)
               enddo
            write(l6,101) l,n1,n2,s,sb*pi
  101       format(' G.Lag.: l ='
     &             ,i2,' n1 =',i3,'  n2 =',i3,3f15.8)
            enddo
            enddo
         enddo
      endif
c
      if (lpr) then
      write(l6,*) ' ****** END GAUPOL *********************************'
      endif
c
      return
c-end-GAUPOL
      end

c======================================================================c
c
      subroutine gdd(lpr)
c
c----------------------------------------------------------------------c
c
c     calculates
c        fmes(x,1)       density-dependent coupling constants
c     and
c        fmes(x,2)       their derivatives
c
c======================================================================c
      use dirqfampar;
      use couplf;
      use dens;

      implicit real*8(a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr

!ab   common /couplf/ ff(MG,4,2)
!ab   common /dens  / ro(MG,4),dro(MG,4)
      common /couplg/ ggmes(4),lmes(4)
      common /dforce/ a_m(4),b_m(4),c_m(4),d_m(4),dsat
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      fun1(x,a,b,c,d) = a*(1+b*(x+d)**2)/(1+c*(x+d)**2)   ! Typel-Wolter
      dun1(x,a,b,c,d) = 2*a*(b-c)*(x+d)/(1+c*(x+d)**2)**2 ! Typel-Wolter

      fun2(x,a,b,c,d) = exp(-a*(x-one))                   ! Typel-Wolter
      dun2(x,a,b,c,d) = -a*exp(-a*(x-one))                ! Typel-Wolter

      fun3(x,a,b,c,d) = a+(b+c*x)*exp(-d*x)               ! DD-PC1
      dun3(x,a,b,c,d) = c*exp(-d*x)-d*(b+c*x)*exp(-d*x)   ! DD-PC1

      if (lpr) then
      write(l6,*) '****** BEGIN GDD ***********************************'
      endif

      if (ipc.eq.0) then
         do i = 1,MG
            x=ro(i,2)/dsat
            do m = 1,2
               ff(i,m,1) = fun1(x,a_m(m),b_m(m),c_m(m),d_m(m))
               ff(i,m,2) = dun1(x,a_m(m),b_m(m),c_m(m),d_m(m))/dsat
            enddo ! m
            do m = 3,4
               ff(i,m,1) = fun2(x,a_m(m),b_m(m),c_m(m),d_m(m))
               ff(i,m,2) = dun2(x,a_m(m),b_m(m),c_m(m),d_m(m))/dsat
            enddo ! m
         enddo   ! i
      elseif (ipc.eq.1) then
         do i = 1,MG
            x=ro(i,2)/dsat
            do m = 1,4
               ff(i,m,1) = fun3(x,a_m(m),b_m(m),c_m(m),d_m(m))
               ff(i,m,2) = dun3(x,a_m(m),b_m(m),c_m(m),d_m(m))/dsat
            enddo ! m
         enddo   ! i
      endif   ! ipc
c
      if (lpr) then
         ix = 5
         do i = 1,ix
            write(6,100) i,(ff(i,m,1),m=1,4)
         enddo   ! i
      endif   ! lpr
c
      if (lpr) then
      write(l6,*) '****** END GDD *************************************'
      endif
c
  100 format(i4,4f12.6)
c
      return
c-end-GDD
      end

c======================================================================c

      subroutine inout(is,lpr)

c======================================================================c
c
c     is = 0  writes nothing to tape
c          1  reads  fields from tape
c          2  writes fields to tape
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use potpot;

      implicit real*8 (a-h,o-z)
!ab   include 'dirqfam.par'
c
      logical lpr
      character*2 nucnam
c
!ab   dimension tz1(2) (not even used...)
c
      common /fermi / ala(2),tz(2)
      common /con_b2/ betac,q0c,cquad,c0,alaq,calcq0,icstr
      common /initia/ inin,inink
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,npr(3),nucnam
!ab   common /potpot/ vps(MG,2),vms(MG,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      if (is.eq.1.and.inin.ne.0) return
c
      if (lpr) then
      write(l6,*) ' ****** BEGIN INOUT ********************************'
      endif
c
c
c---- reading of meson fields from tape:
      if (is.eq.1) then


      !ab
c=======================================================================c
         open( lwin , file = './output/GS_output/dirhb.wel' ,
     &         status = 'old' , form = 'unformatted' );
c=======================================================================c
      !ab


         read(lwin) ng0
         if (ng0.ne.MG) stop 'in INOUT: ngauss wrong'
         read(lwin) ala,alaq
         read(lwin) vps,vms
         close(lwin)
c
         write(l6,*) ' potentials read from tape dirhb.wel:'
      endif
c
c---- writes fields to tape
      if (is.ge.2) then


      !ab
c=======================================================================c
         open( lwou , file = './output/GS_output/dirhb.wel',
     &         status = 'unknown' , form = 'unformatted' );
c=======================================================================c
      !ab


         write(lwou) MG
         write(lwou) ala,alaq
         write(lwou) vps,vms
         close(lwou)
c
         if(lpr) write(l6,*) ' potentials written to tape dirhb.wel'
c
      endif
c
      if (lpr) then
      write(l6,*) ' ****** END INOUT **********************************'
      endif
c
c
      return
c-end-INOUT
      end

c=====================================================================c

      subroutine iter(lpr)

c======================================================================c
c
c     main iteration for the spherical Dirac program
c
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)
c
!ab   logical lpr,lprx
      logical lpr;
      character*2 nucnam
      character*14 text3
!ab   character*27 text1,text2
      !ab
      character*30 text1;
      character*28 text2;
      !ab
c
      common /erwar / ea,rms,betg,gamg
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nneu,npro,nmas,nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /pair  / del(2),spk(2),spk0(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
c
      data spk_min/0.1/
c
c
      text1 = ': Iteration interrupted after '
      text2 = ': Iteration converged after '
      text3 = ' steps   si = '
c
c      CALL zeit(1,'ITER_____',6,.false.)
c
      write(l6,*) '****** BEGIN ITER **********************************'

      ii=0
      do it=1,2
         call gamma(it)
      enddo
      call broyden(.false.)
c
      do ite = 1,maxi
         ii = ite
         write(l6,102) ii,'.It. si = ',si,'  E/A = ',ea,
     &                    ' R = ',rms,' b = ',betg,'  mix =',xmix
         write( 6,102) ii,'.It. si = ',si,'  E/A = ',ea,
     &                    ' R = ',rms,' b = ',betg,'  mix =',xmix
c
c
c
c------- loop over neutrons and protons
         do it = 1,2
c
c---------- calculation of the mean field Gamma
c           call gamma(it)
c
c---------- diagonalization of the Dirac-Bogolibov equation
            call dirhb(it,.false.)
c
c---------- calculation of densities in oscillator basis
            call denssh(it,.false.)
c
         enddo   ! it
c
c------- calculation of new densities in r-space
         call densit(.false.)
c
c------- new coupling constants
         call gdd(.false.)
c
c------- calculation of new fields
         call field(.false.)
c
c------- calculation of the Coulomb potential
         call coulom(.false.)
c
c------- calculation of expectation values
         call expect(.false.)
c
c------- potentials in r-space
         call cstrpot(.false.)
         call poten(.false.)
c
c------- pairing field
         do it = 1,2
            call gamma(it)
            call delta(it,.false.)
         enddo

         call broyden(.false.)
c
c------- check for convergence
         if (ii.gt.2) then
            ic = itestc()
            if (ic.eq.1) goto 20
            if (ic.eq.2) goto 30
         endif
c
      enddo   ! ite
   20 write(6,100) nucnam,nmas,text1,ii,text3,si
      if (l6.ne.6) write(l6,100) nucnam,nmas,text1,ii,text3,si
      goto 40
c
   30 write(6,101) nucnam,nmas,text2,ii,text3,si
      if (l6.ne.6) write(l6,100) nucnam,nmas,text2,ii,text3,si
c
   40 write(l6,*) '****** END ITER ************************************'
      write( 6,*) '****** END ITER ************************************'
c     read*
c
c      CALL zeit(2,'ITER_____',6,.false.)
c
  100 format(1x,68(1h*),/,2x,a2,i4,a27,i4,a14,f17.10,/,1x,68(1h*))
  101 format(1x,a2,i4,a27,i4,a14,f17.10)
  102 format(i3,a,f10.6,3(a,f7.3),a,f5.2)
c
      return
c-end-ITER
      end

c=====================================================================c

      subroutine plot(lpr)

c=====================================================================c
C
C     prepares plot of densities in coordinate space
C
c---------------------------------------------------------------------c
      use dirqfampar;
      use plotmod;
      use dens;
      use gaussh;
      use gaussl;

      implicit real*8 (a-h,o-z)
!ab   include 'dirqfam.par'
      logical lpr
c
!ab   dimension pn(NOX1),zp(0:NGL,0:NGH,3)
c
      common /baspar/ hom,hb0,b0
      common /defbas/ beta0,q,bp,bz
!ab   common /dens/ ro(0:ngh,0:ngl,4),dro(0:ngh,0:ngl,4)
!ab   common /gaussh/ xh(0:NGH),wh(0:NGH),zb(0:NGH)
!ab   common /gaussl/ xl(0:NGL),wl(0:NGL),sxl(0:NGL),rb(0:NGL)
      common /masses/ amu,ames(4)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
c
      if (lpr)
     &write(l6,*) ' ****** BEGIN PLOT ********************************'
c
c     number of points for the plot
      mxplz  = 160
      mxplr  = 80
c     plot step in (fm)
      stplz = 0.1
      stplr = 0.1
      am = amu*hbc

c---- plot for densities:


      !ab
c=======================================================================c
      open( lplo , file = './output/GS_output/dirhb.plo',
     &      status = 'unknown' );
c=======================================================================c
      !ab


!ab   call splin2(rb,zb,NGL,NGH,ro(0,0,2),zp)
      call splin2(rb,zb,NGL,NGH,ro(1,2),zp);
      r=0.d0
      do ir=0,mxplr
         z=-8.d0
         do iz=0,mxplz
!ab         call splint2(r,abs(z),NGL,NGH,rb,zb,ro(0,0,2),zp,rv)
            call splint2(r,abs(z),NGL,NGH,rb,zb,ro(1,2),zp,rv);
            write(lplo,100) r,z,rv
            z=z+stplz
         enddo
        r=r+stplr
       enddo
  100       format(2f10.3,f15.6)
103          format(f10.3,2f15.6)
      close(lplo)
c
      if (lpr)
     &write(l6,*) ' ****** END PLOT **********************************'
      return
C-end-PLOT
      end

c======================================================================c

      subroutine poten(lpr)

c======================================================================c
c
c     CALCULATION OF THE POTENTIALS AT GAUSS-MESHPOINTS
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use constr;
      use coulmb;
      use couplf;
      use dens;
      use fields;
      use potpot;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
      dimension glt(4)
c
      common /baspar/ hom,hb0,b0
!ab   common /constr/ vc(MG,2)
      common /con_b2/ betac,q0c,cquad,c0,alaq,calcq0,icstr
!ab   common /coulmb/ cou(MG),drvp(MG)
      common /coupld/ ddmes(4)
!ab   common /couplf/ ff(MG,4,2)
      common /couplg/ ggmes(4),lmes(4)
!ab   common /dens  / ro(MG,4),dro(MG,4)
!ab   common /fields/ phi(MG,4)
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
!ab   common /potpot/ vps(MG,2),vms(MG,2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp

      if (lpr) then
      write(l6,*) ' ****** BEGIN POTEN ********************************'
      endif
c
      do i = 1,MG
c------- meson-fields
         if (ipc.eq.0) then
            do m = 1,4
               glt(m) = ggmes(m)*ff(i,m,1)*phi(i,m)
            enddo   ! m
c
c           rearangement field
	    re = zero
	    do m = 1,4
               re = re + ggmes(m)*ff(i,m,2)*phi(i,m)*ro(i,m)
	    enddo   ! m
            glt(2) = glt(2) + re
c
c------- point-coupling models
	     elseif (ipc.eq.1) then
	        do m = 1,4
               glt(m) = ggmes(m)*ff(i,m,1)*ro(i,m)
	        enddo   ! m
c
c           derivative terms
	        do m = 1,4
               glt(m) = glt(m) + ddmes(m)*dro(i,m)
	        enddo   ! m
c
c           rearangement field
	        re = zero
	        do m = 1,4
               re = re + ggmes(m)*ff(i,m,2)*ro(i,m)**2
	        enddo   ! m
            glt(2) = glt(2) + half*re
c
         else
            stop 'in POTEN: ipc not properly defined'
         endif   ! ipc
c
         s1 = hbc*(glt(1) - glt(3))                 ! neutron scalar
         s2 = hbc*(glt(1) + glt(3))                 ! proton  scalar
         v1 = hbc*(glt(2) - glt(4))                 ! neutron vector
         v2 = hbc*(glt(2) + glt(4) + cou(i))        ! proton  vector
c
c------- constraining potential
         if (icstr.gt.0) then
            v1 = v1 + vc(i,1)
            v2 = v2 + vc(i,2)
         endif   ! icstr
c
         vps(i,1) = v1 + s1
         vps(i,2) = v2 + s2
         vms(i,1) = v1 - s1
         vms(i,2) = v2 - s2

      enddo   ! i

      if (lpr) then
      write(l6,*) ' ****** END POTEN **********************************'
      endif
c
      return
c-end-POTEN
      end
c======================================================================c

      subroutine pot(i0,nh,n,v,aa)

c======================================================================c
c
c     calculates the potentials V+S and V-S in an axial oscillator basis
c
c     i0   initial point for the quantum number array
c     n    dimension
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use quaosc;
      use herpol;
      use lagpol;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
!ab   character tt*8                                            ! quaosc
c
      dimension aa(nh,nh),v(0:NGH,0:NGL)
c
!ab   common /quaosc/ nt,nz(NTX),nr(NTX),ml(NTX),ms(NTX),np(NTX),tt(NTX)
      common /mathco/ zero,one,two,half,third,pi
!ab   common /herpol/ qh(0:nzx,0:NGH),qh1(0:nzx,0:NGH)
!ab   common /lagpol/ ql(0:2*NRX,0:MLX,0:NGL),ql1(0:2*NRX,0:MLX,0:NGL)
c
      do i2 = 1,n
         nz2 = nz(i0+i2)
         nr2 = nr(i0+i2)
         ml2 = ml(i0+i2)
         do i1 = i2,n
            nz1 = nz(i0+i1)
            nr1 = nr(i0+i1)
            ml1 = ml(i0+i1)
c
            if (ml1.eq.ml2) then
               t = zero
               do il = 0,NGL
                  s = zero
                  do ih = 0,NGH
                     s = s + v(ih,il)*qh(nz1,ih)*qh(nz2,ih)
                  enddo
                  t = t + s*ql(nr1,ml1,il)*ql(nr2,ml2,il)
               enddo
               aa(i1,i2) = t
               aa(i2,i1) = t
            else
               aa(i1,i2) = zero
            endif
c
         enddo
      enddo
c
      return
c-end-POT
      end

c======================================================================c

      subroutine prep(lpr)

c======================================================================c
c
c     preparations
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use blokap;
      use bloosc;
      use gaucor;
      use gaussh;
      use gaussl;

      implicit real*8 (a-h,o-z)
!ab   include 'dirqfam.par'
c
      logical lpr
c
      character tp*1,tis*1,tit*8,tl*1                           ! textex
      character nucnam*2                                        ! nucnuc
!ab   character tb*6                                            ! blokap
c
      common /baspar/ hom,hb0,b0
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
      common /defbas/ beta0,q,bp,bz
!ab   common /gaucor/ ww(MG)
!ab   common /gaussh/ xh(0:NGH),wh(0:NGH),zb(0:NGH)
!ab   common /gaussl/ xl(0:NGL),wl(0:NGL),sxl(0:NGL),rb(0:NGL)
      common /iterat/ si,siold,epsi,xmix,xmix0,xmax,maxi,ii,inxt,iaut
      common /masses/ amu,ames(4)
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nneu,npro,nmas,nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /con_b2/ betac,q0c,cquad,c0,alaq,calcq0,icstr
      common /pair  / del(2),spk(2),spk0(2)
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
c
c
      if (lpr) then
      endif
      write(l6,*) '****** BEGIN PREP **********************************'
c
c
c---- basis parameters
      hb0 = hbc/(two*amu)
      hom = 41.0*amas**(-third)
      if (icm.eq.1) hb0 = hb0*(one - one/amas)
      b0 = sqrt(two*hb0/hom)

      bp    = q**(-1.d0/6.d0)
      bz    = q**(+1.d0/3.d0)
      bpp   = (b0*bp)**2
c
      write(l6,*) ' '
      write(l6,100) ' hom   =',hom
      write(l6,100) ' hb0   =',hb0
      write(l6,100) ' b0    =',b0
      write(l6,100) ' beta0 =',beta0
      write(l6,100) ' bz    =',bz
      write(l6,100) ' bp    =',bp
c

      if (icstr.gt.0) then
         r02=(1.2*amas**third)**2
         q0c = 3*amas*r02*betac/4/pi
         f0 = sqrt(16*pi/5)
         q0c = q0c*f0  ! (2z^2-x^2-y^2)
         c0  = cquad*hom/(hbc*b0**2)
      endif
c
      do ih = 0,NGH
c
c        z-coordinate in fm
         zb(ih) = xh(ih)*b0*bz
c
      enddo  ! NGH
      do il = 0,NGL
c
c        r-coordinate in fm
         rb(il)  = sxl(il)*b0*bp
         do ih = 0,NGH
            ww(1+ih+il*(NGH+1)) = b0**3 * bz*bp*bp * pi * wh(ih)*wl(il)
         enddo
      enddo
c
c
c
      if(icou.eq.1) write(l6,*) 'With Coulomb force'
      if(icm.eq.2) write(l6,*) 'With microscopic c.m. correction'

c---- is the configuration space large enough ?
      n = 0
      do ib = 1,nb
	     n = n + 2*id(ib,1)
      enddo
      if (nneu.gt.n) stop 'in PREP: N is larger than conf.space'
      if (npro.gt.n) stop 'in PREP: Z is larger than conf.space'
c
c
  100 format(a,4f10.6)
!ab 101 format(a,2i4)
!ab 103 format(a,4f10.3)
c
      if (lpr) then
      write(l6,*) '****** END PREP ************************************'
      endif
c
      return
c-end PREP
      end

c======================================================================c

      subroutine reader(lpr)

c======================================================================c
      use dirqfampar;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
      character parname*10                                      ! partyp
      character nucnam*2                                        ! nucnuc
      character tq*1

      CHARACTER*8 date
      CHARACTER*10 time
      CHARACTER*5 zone
      INTEGER*4 VALUES(8)
c
      common /defbas/ beta0,q,bp,bz
      common /defini/ betai
      common /basnnn/ n0f,n0b
      common /baspar/ hom,hb0,b0
      common /con_b2/ betac,q0c,cquad,c0,alaq,calcq0,icstr
      common /fermi / ala(2),tz(2)
      common /initia/ inin,inink
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,nneu,npro,nmas,nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /partyp/ parname
      common /pair  / del(2),spk(2),spk0(2)
      common /physco/ hbc,alphi,r0
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      open(lin,file='dirqfam.dat',status='old')


      !ab
c=======================================================================c
      open( l6 , file = './output/GS_output/dirhb.out' ,
     &      status = 'unknown' );
c=======================================================================c
      !ab

c
      call date_and_time( date, time, zone, values )

c---- basis parameters
      read(lin,'(10x,2i5)') n0f,n0b
c
c---- deformation-parameters of the basis
      read(lin,'(a1,9x,2f10.4)') tq,q
      w = 3*sqrt(5/(16*pi))
      beta0 = q
      q     = exp(w*beta0)
c
c---- initial deformation-parameters
      read(lin,'(10x,f10.4) ') betai
c
c---- Initialization of wavefunctions:
      read(lin,'(10x,2i5)') inin,inink
c
c---- Nucleus under consideration
      read(lin,'(a2,i4)') nucnam,nmas
      call nucleus(2,npro,nucnam)
      nneu = nmas - npro
      amas = nmas
      tz(1) = nneu
      tz(2) = npro
c
c------------------------------------------------------------
c---- pairing
c------------------------------------------------------------
c---- initial gap parameters
      read(lin,*)
      read(lin,'(10x,2f10.4)' ) del

c------------------------------------------------------------
c---- Effective interaction
c------------------------------------------------------------
      read(lin,*)
      read(lin,'(12x,a10)') parname
c
c------------------------------------------------------------
c---- constraining fields
c------------------------------------------------------------
      read(lin,*)
      read(lin,'(10x,i5)') icstr
      read(lin,'(10x,f10.4) ') betac
      read(lin,'(10x,f10.4) ') cquad


!ab 104 format(10x,4i3)
!ab 117 format(a,4i3)

      close(lin)

      write(6,'(a)')
      write(6,'(a)') '  ******************************************  '
      write(6,'(a)') '  *           Fortran 77 Code              *  '
      write(6,'(a)') '  *     Axially symmetric H.O. basis       *  '
      write(6,'(a)') '  * Dirac-Hartree-Bogoliubov calculation   *  '
      write(6,'(a)') '  *        with density-dependent          *  '
      write(6,'(a)') '  * meson-exchange or point coupling force *  '
      write(6,'(a)') '  *          and separable pairing         *  '
      write(6,'(a)') '  * -------------------------------------- *  '
      write(6,'(a)') '  *      Niksic, Paar, Vretenar, Ring      *  '
      write(6,'(a,i2,a,i2,a,i4,a,a2,a,a2,a,a2,a)')
     &               '  *          ',values(3),'/',values(2),'/',
     &          values(1),'/',time(1:2),':',time(3:4),':',time(5:6),
     &          '           *'
      write(6,'(a,a2,i4,a,i3,a,i3)')
     &              '  *       ',
     &    nucnam,nmas,'  N = ',nmas-npro,'  Z = ',npro
      write(6,'(a,16x,a10)') '  *',parname
      write(6,'(a)') '  ******************************************  '

      write(l6,'(a)')
      write(l6,'(a)') '  ******************************************  '
      write(l6,'(a)') '  *           Fortran 77 Code              *  '
      write(l6,'(a)') '  *     Axially symmetric H.O. basis       *  '
      write(l6,'(a)') '  * Dirac-Hartree-Bogoliubov calculation   *  '
      write(l6,'(a)') '  *        with density-dependent          *  '
      write(l6,'(a)') '  * meson-exchange or point coupling force *  '
      write(l6,'(a)') '  *          and separable pairing         *  '
      write(l6,'(a)') '  * -------------------------------------- *  '
      write(l6,'(a)') '  *      Niksic, Paar, Vretenar, Ring      *  '
      write(l6,'(a,i2,a,i2,a,i4,a,a2,a,a2,a,a2,a)')
     &               '  *          ',values(3),'/',values(2),'/',
     &          values(1),'/',time(1:2),':',time(3:4),':',time(5:6),
     &          '           *'
      write(l6,'(a,a2,i4,a,i3,a,i3)')
     &              '  *       ',
     &    nucnam,nmas,'  N = ',nmas-npro,'  Z = ',npro
      write(l6,'(a,16x,a10)') '  *',parname
      write(l6,'(a)') '  ******************************************  '

      write(l6,*) '****** BEGIN READER ********************************'
      write(l6,'(a,2i5)') ' Number of oscillator shells : ',n0f,n0b
      write(l6,'(a,f10.4) ') ' Deformation bas: beta0      : ',beta0
      write(l6,'(a,f10.4)' ) ' Deformation ini: betai      : ',betai
      write(l6,'(a,2i5)') ' Initialization   inin,inink : ',inin,inink
      write(l6,115) ' Nucleus',': ',nneu,npro,nucnam,nmas
  115 format(a,21x,a,2i5,3x,a,i4)
      write(l6,'(a,2f10.4) ') ' Initial Gap Parameters      : ',del
      write(l6,'(a,a10)') ' Parameter set               : ',parname
      if (icstr.eq.0) then
          write(l6,'(a)') ' Calculation without quadrupole constraint'
      else
          write(l6,'(a)') ' Calculation with quadrupole constraint'
          write(l6,'(a,f10.4)' ) ' contraint beta       : ',betac
         write(l6,'(a,f10.4)' ) ' spring strength      : ',cquad
      endif
      write(l6,*) '****** END READER **********************************'
      if (lpr) then
      endif

      return
c-end-READER
      end

c======================================================================c

      subroutine resu(lpr)

c======================================================================c
      use dirqfampar;
      use blodir;
      use blokap;
      use bloosc;
      use coulmb;
      use fields;
      use quaosc;
      use waveuv;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
      character tp*1,tis*1,tit*8,tl*1                           ! textex
!ab   character tb*6                                            ! blokap
!ab   character tt*8                                            ! quaosc
      character tph*1
c
!ab   common /blodir/ ka(NBX,4),kd(NBX,4)
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
!ab   common /coulmb/ cou(MG),drvp(MG)
      common /fermi / ala(2),tz(2)
!ab   common /fields/ phi(MG,4)
      common /mathco/ zero,one,two,half,third,pi
      common /optopt/ itx,icm,icou,ipc,inl,idd
!ab   common /quaosc/ nt,nnn(NTX,5),tt(NTX)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /textex/ tp(2),tis(2),tit(2),tl(0:30)
!ab   common /waveuv/ fguv(NHBX,KX,4),equ(KX,4)

      data emaxqp/80.d0/

c
c---- prints expectation values
      call expect(.true.)

      if(.not.lpr) return
      write(l6,*) '****** BEGIN RESU **********************************'
c
c
c---- quasi-particle energies
      do it = 1,2
         write(l6,100) tit(it)
         write(l6,101) 'p/h','K pi','[nz,nr,ml]','smax','Eqp','uu','vv'
      do ib = 1,nb
         nf = id(ib,1)
         ng = id(ib,2)
         nh = nf + ng
         k1 = ka(ib,it) + 1
         k2 = ka(ib,it) + kd(ib,it)
         i0 = ia(ib,1)
c
         do k = k1,k2
c
c           search for main oscillator component
            call maxa(nh,fguv(1,k,it),imax,fx)
            su=zero
            sv=zero
            do n = 1,nh
               su = su + fguv(   n,k,it)**2
               sv = sv + fguv(nh+n,k,it)**2
            enddo
            if (sv.lt.half) tph='p'
            if (sv.ge.half) tph='h'
c
            write(l6,102) k,tph,tb(ib),tt(i0+imax),fx,equ(k,it),su,sv
         enddo  ! k
      enddo   ! ib
      enddo   ! it
c
      write(l6,*) ' ****** END RESU ***********************************'
c
  100 format(//,' quasi-particle properties: ',a,/,1x,66(1h-))
  101 format(5x,a,2x,a,a,2x,a,5x,a,7x,a,7x,a)
  102 format(i4,2x,a1,2x,a6,a8,f7.2,5f10.4)
c
      return
c-end-RESU
      end

c======================================================================c

      subroutine singd(lpr)

c======================================================================c
c
c     calculates matrix elements V_nz and V_nr' for TMR-pairing
c     in the axially deformed oscillator basis
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use singdmod;
      use blokap;
      use bloosc;
      use tmrwnn;
      use vvvikf;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
!ab   character tb*6                                            ! blokap
c
!ab   dimension pnosc(0:N0FX),vnz(1:MVX,0:N0FX),vnr(1:MVX,0:N0FX)
!ab   dimension wn(MVX)
c
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /tmrpar/ gl(2),ga
!ab   common /tmrwnn/ wnn(MVX,NNNX),nnmax
!ab   common /vvvikf/ mv,ipos(NBX),nib(MVX),nni(2,MVX)
c
      data eps/1.d-20/
c
      if (lpr) then
      write(l6,*) '****** BEGIN VPAIR TMR *****************************'
      endif
c
c
c----------------------------------------------------------------------c
c     calculate VNZ and VNR
c----------------------------------------------------------------------c
      call vvnz(nnzm,pnosc,vnz,.false.)
      call vvnp(nnrm,pnosc,vnr,.false.)
c
c----------------------------------------------------------------------c
c    calculate the single particle matrix elements WNN
c----------------------------------------------------------------------c
c     write(6,*) 'nnzm =',nnzm
c     write(6,*) 'nnrm =',nnrm
      nn = 0
      do nnz = 0,nnzm
      do nnr = 0,nnrm
         nn = nn + 1
         if (nn.gt.NNNX) stop 'in VPAIR: NNNX too small'
         smax = zero
         do i = 1,mv
            wn(i) = vnz(i,nnz)*vnr(i,nnr)
            smax = max(smax,abs(wn(i)))
         enddo   ! i
         if (smax.lt.eps) then
            nn = nn - 1
         else
            do i = 1,mv
               wnn(i,nn) = wn(i)
            enddo   ! i
         endif  ! smax < eps
      enddo   ! nnz
      enddo   ! nnr
      nnmax = nn
c     write(6,*)  'VPAIR: nnmax = ',nnmax,(nnzm+1)*(nnrm+1)
c
c---- Printout WNN
      if (lpr) then
         ix = mv
         ix = 10
         na = 1
         nx = 5
         nx = min(nnmax,nx)
         write(6,*) 'WNN',MVX,NNNX
         write(6,*) 'WNN',mv,nnmax
         write(6,*) 'WNN',ix,na,nx
         do i = 1,ix
            write(6,100) i,(wnn(i,n),n=na,nx)
         enddo   ! i
      endif   ! lpr
c
      if (lpr) then
      write(l6,*) '****** END VPAIR ***********************************'
      endif
c
  100 format(i5,10f8.4)
c
      return
c-end-VPAIR
      end
c======================================================================c

      subroutine vvnz(nnzm,pnosc,vnz,lpr)

c======================================================================c
c
c     calculates matrix elements V_nz for TMR-pairing
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use blokap;
      use bloosc;
      use gfvmod;
      use quaosc;
      use vvvikf;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
!ab   character tt*8                                            ! quaosc
!ab   character tb*6                                            ! blokap
c
      dimension pnosc(0:N0FX),vnz(1:MVX,0:N0FX)
c
      common /basnnn/ n0f,n0b
      common /baspar/ hom,hb0,b0
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
      common /defbas/ beta0,q,bp,bz
!ab   common /gfvsq / sq(0:IGFV)
!ab   common /quaosc/ nt,nz(NTX),nr(NTX),ml(NTX),ms(NTX),np(NTX),tt(NTX)
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /tmrpar/ gl(2),ga
!ab   common /vvvikf/ mv,ipos(NBX),nib(MVX),nni(2,MVX)
c
      data eps/1.d-20/
c
      if (lpr) then
      write(l6,*) '****** BEGIN VVNZ TMR ******************************'
      endif
c
      b0z  = b0*bz
c
      call pnoscz(n0f,sqrt(ga),b0z,pnosc)
c
c----------------------------------------------------------------------c
c    calculate VNZ
c----------------------------------------------------------------------c
      call mzero(MVX,mv,n0f+1,vnz)
      il = 0
      do ib = 1,nb
         nf = id(ib,1)
         i0f= ia(ib,1)
      do n2   = 1,nf
         nz2 = nz(i0f+n2)
      do n1 = n2,nf
         nz1 = nz(i0f+n1)
         il  = il + 1
c
         if (mod(nz1+nz2,2).eq.0) then
            do nn = 0,nz1+nz2,2
               n  = nz1 + nz2 - nn
               nh  = n/2
               nnh = nn/2
               if (nh .gt.n0f) stop 'in VPAIR: n too large'
               if (nnh.gt.n0f) stop 'in VPAIR: nn too large'
c   talmos1 is Eq. 17 from the paper
               vnz(il,nnh) = pnosc(nh)*talmos1(nz1,nz2,nn,n)
            enddo   ! nn
         endif
c
      enddo   ! n1
      enddo   ! n2
      enddo   ! ib
      call icheck_n1nen2(il,mv,'in VPAIR: mv wrong')
c
c----------------------------------------------------------------------c
c---- calculate the maximal nnz:  nnzm
c----------------------------------------------------------------------c
      nnzm =0
      do nn = 0,n0f
         nx = 0
         do i = 1,mv
            if (abs(vnz(i,nn)).gt.eps) then
               nx = nn
            endif
         enddo   ! i
         nnzm = max(nnzm,nx)
      enddo   ! nn
c
c----------------------------------------------------------------------c
c---- Printout VNZ
c----------------------------------------------------------------------c
      if (lpr) then
         na = 0
         nx = 6
         nx = min(nx,nnzm)
         ix = 20
         ix = min(ix,mv)
         write(6,*) 'VNZ',ix,mv,na,nx
         write(6,101) (n,n=na,nx,2)
         do i = 1,ix
            write(6,100) i,(vnz(i,n),n=na,nx)
         enddo   ! i
      endif
c
      if (lpr) then
      write(l6,*) '****** END VVNZ ************************************'
      endif
c
  100 format(i4,20f10.6)
  101 format(4x,9i10)
c
      return
c-end-VVNZ
      end
c======================================================================c

      subroutine pnoscz(nm,a,b,pnosc)

c======================================================================c
      use dirqfampar;
      use gfvmod;

      implicit real*8 (a-h,o-z)
!ab   include 'dirqfam.par'
c
      dimension pnosc(0:nm)
c
      common /baspar/ hom,hb0,b0
!ab   common /gfvsq / sq(0:IGFV)
!ab   common /gfvwf / wf(0:IGFV)
!ab   common /gfvfi / fi(0:IGFV)
      common /mathco/ zero,one,two,half,third,pi
      common /tmrpar/ gl(2),ga
c
      alpha = a/b
c      br = b*sq(2)
c
      s0 = (2*pi)**0.25d0
      s0 = s0*sqrt(b/(a*a+b*b))
      s1 = (a*a-b*b)/(a*a+b*b)
      do nh = 0,nm
         n = 2*nh
	 pnosc(nh) = s0 * s1**nh * wf(n) / 2**nh * fi(nh)

      enddo    ! n
c     write(6,100) 'PNZ',nm,(pnosc(n),n=0,5)
c 100 format(a,i3,20f10.6)
c
      return
c-end-PNOSCZ
      end
c======================================================================c

      subroutine vvnp(nnrm,pnocz,vn,lpr)

c======================================================================c
c
c     matrix elements V_np (perpendicular direction) for TMR-pairing
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use vvnpmod;
      use blokap;
      use bloosc;
      use gfvmod;
      use quaosc;
      use vvvikf;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
!ab   character tt*8                                            ! quaosc
!ab   character tb*6                                            ! blokap
c
!ab   dimension pnosc(0:N0FX)
      dimension vn(1:MVX,0:N0FX)
      dimension pnocz(0:N0FX)
c
      common /basnnn/ n0f,n0b
      common /baspar/ hom,hb0,b0
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
      common /defbas/ beta0,q,bp,bz
!ab   common /gfvsq / sq(0:IGFV)
!ab   common /quaosc/ nt,nz(NTX),nr(NTX),ml(NTX),ms(NTX),np(NTX),tt(NTX)
      common /mathco/ zero,one,two,half,third,pi
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
      common /tmrpar/ gl(2),ga
!ab   common /vvvikf/ mv,ipos(NBX),nib(MVX),nni(2,MVX)
c
      data eps/1.d-20/
c
      if (lpr) then
      write(l6,*) '****** BEGIN VVNP **********************************'
      endif
c

      b0p  = b0*bp
c
      call pnoscp(n0f,sqrt(ga),b0p,pnosc)
c
      call mzero(MVX,mv,n0f+1,vn)
      il = 0
      do ib = 1,nb
         nf = id(ib,1)
         i0f= ia(ib,1)
      do n2 = 1,nf
         nr2 = nr(i0f+n2)
         ml2 = ml(i0f+n2)
      do n1 = n2,nf
         nr1 = nr(i0f+n1)
         ml1 = ml(i0f+n1)
         il  = il + 1
c
         if (ml1.eq.ml2) then
            do nn = 0,nr1+nr2+ml1
               n = nr1+nr2+ml1-nn
               if (nn.gt.n0f) stop 'in VPAIR: nn too large'
               if (n .gt.n0f) stop 'in VPAIR: n  too large'
               vn(il,nn) = pnosc(n)*talmos2(ml1,nr1,nr2,0,nn,n)
            enddo   ! nn
         endif
c
      enddo   ! n1
      enddo   ! n2
      enddo   ! ib
      call icheck_n1nen2(il,mv,'in VVNP: mv wrong')
c
c----------------------------------------------------------------------c
c---- calculate the maximal nnr:  nnrm
c----------------------------------------------------------------------c
      nnrm = 0
      do nn = 0,n0f
         nx = 0
         do i = 1,mv
            if (abs(vn(i,nn)).gt.eps) then
               nx = nn
            endif
         enddo   ! i
         nnrm = max(nnrm,nx)
      enddo   ! nn
c
c----------------------------------------------------------------------c
c---- Printout VN
c----------------------------------------------------------------------c
      if (lpr) then
         ix = 20
         ix = min(ix,mv)
         na = 0
         nx = 6
         nx = min(nx,nnr0)
         write(6,*) 'VNR',ix,mv,na,nx
         write(6,101) (n,n=na,nx)
         do i = 1,ix
            write(6,100) i,(vn(i,n),n=na,nx)
         enddo   ! i
      endif
c
      if (lpr) then
      write(l6,*) '****** END VVNP ************************************'
      endif
c
  100 format(i4,9f10.6)
  101 format(4x,9i10)
c
      return
c-end-VVNP
      end
c======================================================================c

      subroutine pnoscp(nm,a,b,pnosc)

c======================================================================c
      use dirqfampar;
      use gfvmod;

      implicit real*8 (a-h,o-z)
!ab   include 'dirqfam.par'
c
      dimension pnosc(0:nm)
c
      common /basnnn/ n0f,n0b
      common /baspar/ hom,hb0,b0
      common /defbas/ beta0,q,bp,bz
!ab   common /gfvsq / sq(0:IGFV)
!ab   common /gfvwfi/ wfi(0:IGFV)
      common /mathco/ zero,one,two,half,third,pi
!ab   common /rmesh / rm(MR),wr(MR)
c
      call icheck_n1len2(nm,n0f,' in PNOSCP: nm too large')
c
c      br    = b*sq(2)
c      alpha = a/br
      do n = 0,nm
c         s = zero
c         do i = 1,NGL
c            x = sxl(i)
c            y = xl(i)/(4*alpha**2)
c            s = s + osc2(n,0,x)*exp(-y)*wl(i)
c         enddo   ! i
c         pnosc(n) = s/(8*pi*br*alpha**2)
          f1=b/(b**2+a**2)
          f2=((b**2-a**2)/(b**2+a**2))**n
          pnosc(n) = f1*f2/2/pi
      enddo    ! n
c
!ab 100 format(a,i3,3f10.6)
c
      return
c-end-PNOSCP
      end
c=======================================================================c
c
      real*8 function talmos1(n1,n2,n3,n4)
c
c=======================================================================c
c
c	1-dimensional Talmi-Moshinsky bracket: <n1 , n2 | n3 , n4 >
c
c	quantum number n_i start from zero: n = 0,1,2,....
c
c-----------------------------------------------------------------------c
      use dirqfampar;
      use gfvmod;

      implicit real *8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
!ab   common /gfviv / iv(-IGFV:IGFV)
!ab   common /gfvwf / wf(0:IGFV)
!ab   common /gfvwfi/ wfi(0:IGFV)
!ab   common /bin0  / bin(0:IGFV,0:IGFV)
      common /mathco/ zero,one,two,half,third,pi
c
      talmos1 = zero
      if (n1+n2.ne.n3+n4) return
      f = wf(n1)*wfi(n3)*wf(n2)*wfi(n4)/sqrt(two**(n3+n4))
      s = zero
      do m4 = 0,n4
         m3 = n1 - n4 + m4
         if (m3.ge.0.and.n3.ge.m3) then
            s = s + iv(m4)*bin(n3,m3)*bin(n4,m4)
         endif
      enddo
      talmos1 = f*s
c
      return
c-ene-TALMOS1
      end
c=======================================================================c
c
      real*8 function talmos2(im1,in1,in2,im3,in3,in4)
c
c=======================================================================c
c
c	2d-moshinsky bracket:
c
c	<n1 m1, n2 -m1 | n3 m3, n4 -m3>
c
c	radial quantum number start from zero: n=0,1,2,....
c	orbital angular momentum m1=-m2>0,m3=-m4>0
c
c-----------------------------------------------------------------------c
C
c     iv(n)  =  (-1)**n
c     fak(n) =  n!
c     fi(n)  =  1/n!
c     wf(n)  =  sqrt(n!)
c     wfi(n) =  1/sqrt(n!)
C
C-----------------------------------------------------------------------c
      use dirqfampar;
      use gfvmod;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
!ab   common /gfviv / iv(-IGFV:IGFV)
!ab   common /gfvfak/ fak(0:IGFV)
!ab   common /gfvfi / fi(0:IGFV)
!ab   common /gfvwf / wf(0:IGFV)
!ab   common /gfvwfi/ wfi(0:IGFV)
      common /mathco/ zero,one,two,half,third,pi
c
      im2 = -im1
      im4 = -im3
c
      talmos2 = zero
c
      nn12 = 2*in1+iabs(im1)+2*in2+iabs(im2)
      nn34 = 2*in3+iabs(im3)+2*in4+iabs(im4)
      if (im1.lt.0)     stop 'in TALMOS2: m1 < 0'
      if (im3.lt.0)     stop 'in TALMOS2: m3 < 0'
      if (nn12.ne.nn34) return
c
      s34 = zero
c
      if (im1.lt.im2) then
         n1 = in2
         n2 = in1
         m1 = im2
         m2 = im1
      else
         n1 = in1
         n2 = in2
	 m1 = im1
         m2 = im2
      endif
      n3 = in3
      n4 = in4
      m3 = im3
      m4 = im4
      nn1 = 2*n1 + abs(m1)
      nn2 = 2*n2 + abs(m2)
      nn3 = 2*n3 + abs(m3)
      nn4 = 2*n4 + abs(m4)
      if (m3.gt.m4) then
         if (n1+n2.ne.n3+n4+m2-m4) return
      else
         if (n1+n2.ne.n3+n4+m2-m3) return
      endif
      prout = iv(n3+n4-n1-n2)/sqrt(two**(nn3+nn4))*
     &	wf(n1)*wf(n1+abs(m1))*
     &	wf(n2)*wf(n2+abs(m2))*
     &	wfi(n3)*wfi(n3+abs(m3))*
     &	wfi(n4)*wfi(n4+abs(m4))
c
      sn3 = zero
      sn4 = zero
      sm3 = zero
      sm4 = zero
c
      do i3 = 0,n3
      do j3 = 0,n3
      do k3 = 0,n3
c
         l3 = n3 - i3 - j3 - k3
c
      do i4 = 0,n4
      do j4 = 0,n4
      do 20 k4 = 0,n4
c
         l4 = n4 - i4 - j4 - k4
c
         if (l3.lt.0.or.l4.lt.0) goto 20
c
         do it3 = 0,abs(m3)
         do 10 it4 = 0,abs(m4)
c
         if (m3.gt.m4) then
            if (i3+i4+j3+j4+it3.ne.n2) goto 10
            if (j3+j4.ne.m2+k3+k4-it3+it4) goto 10
            if (l3+l4.ne.n3+n4-n2-k3-k4+it3) goto 10
            sn3 = fak(n3)*fi(l3)*fi(i3)*fi(j3)*fi(k3)
            sn4 = iv(j4+k4)*fak(n4)*fi(l4)*fi(i4)*fi(j4)*fi(k4)
            sm3 = fak(abs(m3))*fi(it3)*fi(abs(m3)-it3)
            sm4 = iv(it4)*fak(abs(m4))*fi(it4)*fi(abs(m4)-it4)
            s34 = s34 + sn3*sn4*sm3*sm4
c
         else
            if (i3+i4+j3+j4+it4.ne.n2) goto 10
            if (j3+j4.ne.m2+k3+k4+it3-it4) goto 10
            if (l3+l4.ne.n3+n4-n2-k3-k4+it4) goto 10
            sn3 = fak(n3)*fi(l3)*fi(i3)*fi(j3)*fi(k3)
            sn4 = iv(j4+k4)*fak(n4)*fi(l4)*fi(i4)*fi(j4)*fi(k4)
            sm3 = fak(abs(m3))*fi(it3)*fi(abs(m3)-it3)
            sm4 = iv(it4)*fak(abs(m4))*fi(it4)*fi(abs(m4)-it4)
            s34 = s34 + sn3*sn4*sm3*sm4
c
         endif ! m3,m4

   10 continue   ! it4 modified by tamara
      enddo   ! it3
c
   20 continue   ! k2
      enddo   ! j2
      enddo   ! i2
c
      enddo   ! k1
      enddo   ! j1
      enddo   ! i1
c
      talmos2 = s34*prout
c
      return
c-end-TALMOS2
      end

c======================================================================c

      subroutine singf(lpr)

c======================================================================c
c
c     calculates single particle matrix elements for Fermions
c     in the zylindrical oscillator basis
c
c----------------------------------------------------------------------c
      use dirqfampar;
      use bloosc;
      use blokap;
      use single;
      use single2;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
!ab   character tb*6                                            ! blokap
c
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
!ab   common /single/ sp(nfgx,NBX)
!ab   common /single2/ sprr(NDDX,2*NBX)
c
      if (lpr) then
      write(l6,*) ' ****** BEGIN SINGF ********************************'
      endif
c
      do ib = 1,nb
         nf = id(ib,1)
         ng = id(ib,2)
c
c        SIGMA*P
c----------------
         call sigp(nf,ng,ib,sp(1,ib),lpr)
         call sigrr(nf,1,ib,sprr(1,ib),lpr)
         call sigrr(ng,2,ib,sprr(1,nbx+ib),lpr)
      enddo
C
      if (lpr) then
      write(l6,*) ' ****** END SINGF **********************************'
      endif
c
      return
c-end-SINGF
      end
c=====================================================================c

      subroutine sigp(nf,ng,ib,aa,lpr)

c=====================================================================c
      use dirqfampar;
      use blokap;
      use bloosc;
      use quaosc;
      use gaussl;
      use gfvmod;
      use herpol;
      use lagpol;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
!ab   character tt*8                                            ! quaosc
!ab   character tb*6                                            ! blokap
c
      dimension aa(ng,nf)
c
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
      common /defbas/ beta0,q,bp,bz
!ab   common /quaosc/ nt,nz(NTX),nr(NTX),ml(NTX),ms(NTX),np(NTX),tt(NTX)
!ab   common /gaussl/ xl(0:NGL),wl(0:NGL),sxl(0:NGL),rb(0:NGL)
!ab   common /gfviv / iv(-IGFV:IGFV)
!ab   common /gfvsq / sq(0:IGFV)
!ab   common /gfvsqi/ sqi(0:IGFV)
      common /mathco/ zero,one,two,half,third,pi
!ab   common /herpol/ qh(0:NZX,0:NGH),qh1(0:NZX,0:NGH)
!ab   common /lagpol/ ql(0:2*NRX,0:mlx,0:NGL),ql1(0:2*NRX,0:mlx,0:NGL)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      fz  = one/bz
      fp  = one/bp
      k   = iabs(kb(ib))
      i0f = ia(ib,1)
      i0g = ia(ib,2)
      do i2 = 1,nf
         nz2 = nz(i0f+i2)
         nr2 = nr(i0f+i2)
         ml2 = ml(i0f+i2)
         ms2 = 2*(k - ml2) - 1
         do i1 = 1,ng
            nz1 = nz(i0g+i1)
            nr1 = nr(i0g+i1)
            ml1 = ml(i0g+i1)
            ms1 = 2*(k - ml1) - 1
c
            s = zero
c
c           s1 = s2
            if (ms1.eq.ms2) then
               if (nr1.eq.nr2) then
                  if (nz1.eq.nz2+1) s = -ms1*fz*sq(nz1)*sqi(2)
                  if (nz1.eq.nz2-1) s = +ms1*fz*sq(nz2)*sqi(2)
               endif
c
c           s1 ne s2
            else
               if (nz1.eq.nz2) then
                  if (ml1.eq.ml2+1) then
                     m = -ml2
                  else
                     m = +ml2
                  endif
                  do il = 1,NGL
                     s = s + ql(nr1,ml1,il)*
     &                (ql1(nr2,ml2,il) + m*ql(nr2,ml2,il)/sxl(il))
                  enddo
                  s = s*fp
               endif
            endif
            aa(i1,i2) = -s
c
            nn1 = nz1 + 2*nr1 + iabs(ml1)
            nn2 = nz2 + 2*nr2 + iabs(ml2)
         enddo   ! i1
c        check for matrix elements connecting N->N+1
         i = 0
         do i1 = 1,ng
            if (abs(aa(i1,i2)).gt.1.d-5) i = 1
         enddo   ! i1
         if (i.eq.0) then
            do i1 = 1,ng
               write(6,*) i1,i2,aa(i1,i2)
            enddo   ! i1
            stop 'in SIGP: N->N+1 connection wrong'
         endif
      enddo   ! i2
c
      if (lpr) then
         write(l6,'(//,a)') tb(ib)
         call aprint(1,3,6,ng,ng,nf,aa,tt(i0g+1),tt(i0f+1),'Sigma * P')
         read*
      endif
c
      return
c-end-SIGP
      end
c=====================================================================c

      subroutine sigrr(nd,idd,ib,aa,lpr)

c=====================================================================c
      use dirqfampar;
      use blokap;
      use bloosc;
      use quaosc;
      use gaussl;
      use gfvmod;
      use herpol;
      use lagpol;

      implicit real*8 (a-h,o-z)
c
!ab   include 'dirqfam.par'
c
      logical lpr
c
!ab   character tt*8                                            ! quaosc
!ab   character tb*6                                            ! blokap
c
      dimension aa(nd,nd)
c
!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
      common /defbas/ beta0,q,bp,bz
      common /baspar/ hom,hb0,b0
!ab   common /quaosc/ nt,nz(NTX),nr(NTX),ml(NTX),ms(NTX),np(NTX),tt(NTX)
!ab   common /gaussl/ xl(0:NGL),wl(0:NGL),sxl(0:NGL),rb(0:NGL)
!ab   common /gfviv / iv(-IGFV:IGFV)
!ab   common /gfvsq / sq(0:IGFV)
!ab   common /gfvsqi/ sqi(0:IGFV)
!ab   common /gfvwf / wf(0:IGFV)
!ab   common /gfvwfi/ wfi(0:IGFV)
      common /mathco/ zero,one,two,half,third,pi
!ab   common /herpol/ qh(0:NZX,0:NGH),qh1(0:NZX,0:NGH)
!ab   common /lagpol/ ql(0:2*NRX,0:mlx,0:NGL),ql1(0:2*NRX,0:mlx,0:NGL)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
      k   = iabs(kb(ib))
      i0d = ia(ib,idd)
      do i2 = 1,nd
         nz2 = nz(i0d+i2)
         nr2 = nr(i0d+i2)
         ml2 = ml(i0d+i2)
         ms2 = 2*(k - ml2) - 1
         fz2 = pi**(-0.25d0)*2.d0**(-nz2)*wfi(nz2)
         fr2 = wf(nr2)*wfi(nr2+ml2)
         do i1 = 1,nd
            nz1 = nz(i0d+i1)
            nr1 = nr(i0d+i1)
            ml1 = ml(i0d+i1)
            ms1 = 2*(k - ml1) - 1
            fz1 = pi**(-0.25d0)*2.d0**(-nz1)*wfi(nz1)
            fr1 = wfi(nr1)*wf(nr1+ml1)

            fac = fz2*fr2*fz1*fr1
c
            szz = zero
            srp = zero
c
c           s1 = s2
            if (ms1.eq.ms2) then
               if (nr1.eq.nr2) then
                  if (nz1.eq.nz2+2) szz = sq(nz2+1)*sq(nz2+2)/2
                  if (nz1.eq.nz2-2) szz = sq(nz2)*sq(nz2-1)/2
                  if (nz1.eq.nz2)   szz = (nz1+0.5d0)
               endif
               if (nz1.eq.nz2) then

               !ab
c=======================================================================c
                  if (nr1.eq.nr2+1) then
                      srp = -DSQRT(DBLE(nr1*(nr1+ml1)));
                  endif
                  if (nr1.eq.nr2-1) then
                      srp = -DSQRT(DBLE((nr1+1)*(nr1+ml1+1)));
                  endif
c=======================================================================c
               !ab

                  if (nr1.eq.nr2)   srp = (nr1+1)
               endif
            endif
            aa(i1,i2) = b0**2*(szz+srp)
c
            nn1 = nz1 + 2*nr1 + iabs(ml1)
            nn2 = nz2 + 2*nr2 + iabs(ml2)
         enddo   ! i1

      enddo   ! i2
c

      if (lpr) then
!ab      write(l6,'(//,a)') tb(ib)
!ab      call aprint(1,3,6,ng,ng,nf,aa,tt(i0g+1),tt(i0f+1),'Sigma * P')
!ab      read*
      endif

c
      return
c-end-SIGR
      end

c=================== BEGIN SPLINE LIB =========================================
c
c****************************************************************
c
c    To calculate the interpolation value at (xx,yy)
c    before this one, you should run splin2 at first
c
c     (xx,yy) :  the gaussian mesh points
c     m,n     :  dimension of r_perp and z
c     x,y,z,zp:  input function, zp from subroutine splin2
c     ef      :  the output function
c
c****************************************************************
c
      subroutine splint2(xx,yy,m,n,x,y,z,zp,ef)
c
      use dirqfampar;

      implicit real*8(a-h,o-z)
!ab   include 'dirqfam.par'
c
      dimension x(0:m),y(0:n),z(0:m,0:n),zp(0:m,0:n,3)
c
c---- find the location of (xx,yy) ------------------------
c
      ilo = 0
      ihi = m
1     if ((ihi-ilo).gt.1) then
         k = (ihi+ilo)/2
         if (x(k).gt.xx) then
            ihi = k
         else
            ilo = k
         endif
         goto 1
      endif
c
      jlo = 0
      jhi = n
2     if ((jhi-jlo).gt.1) then
         k = (jhi+jlo)/2
         if (y(k).gt.yy) then
            jhi = k
         else
            jlo = k
         endif
         goto 2
      endif
c
c---- cal. the interpolation value at (xx,yy) -------------------
c
      vmlo   = fu(yy,y(jlo),y(jhi),z(ilo,jlo),z(ilo,jhi),
     *            zp(ilo,jlo,2),zp(ilo,jhi,2))
      vmhi   = fu(yy,y(jlo),y(jhi),z(ihi,jlo),z(ihi,jhi),
     *            zp(ihi,jlo,2),zp(ihi,jhi,2))
      vm2xlo = fu(yy,y(jlo),y(jhi),zp(ilo,jlo,1),zp(ilo,jhi,1),
     *            zp(ilo,jlo,3),zp(ilo,jhi,3))
      vm2xhi = fu(yy,y(jlo),y(jhi),zp(ihi,jlo,1),zp(ihi,jhi,1),
     *            zp(ihi,jlo,3),zp(ihi,jhi,3))
c
      ef = fu(xx,x(ilo),x(ihi),vmlo,vmhi,vm2xlo,vm2xhi)
c
      return
      end
c****************************************************************
c
c    To do the 2D spline interpolation
c
c     ng :  the index for different parameters
c     x,y:  data points
c     m,n:  dimension of \beta and \gamma
c     z  :  the 2D function (the parameters, like V_{coll} ...
c     zp :  The derivative of z
c
c****************************************************************
c
      subroutine splin2(x,y,m,n,z,zp)
c
      use dirqfampar;

      implicit real*8(a-h,o-z)
!ab   include 'dirqfam.par'
c
      dimension x(0:m),y(0:n),z(0:m,0:n),zp(0:m,0:n,3)
      dimension temp1(0:m),temp2(0:n)
c
c------- do spline for r_perp direction -------------------
c
         do 30 j=0, n
c
            do 10 i=0, m
               temp1(i)=z(i,j)
   10       continue
c
            call derlim2(x,temp1,m,derivp,derivk)
            call spline(x,temp1,temp2,m,derivp,derivk)
c
            do 20 i=0, m
               zp(i,j,1)=temp2(i)
   20       continue
c
   30    continue
c-------------------------------------
c
c------- do spline for z direction ------------------
c
         do 60 i=0, m
c
            do 40 j=0, n
               temp1(j)=z(i,j)
   40       continue
c
            call derlim2(y,temp1,n,derivp,derivk)
            call spline(y,temp1,temp2,n,derivp,derivk)
c
            do 50 j=0, n
               zp(i,j,2)=temp2(j)
   50       continue
c
   60    continue
c-------------------------------------
c
c------- do spline of d^2f/dg^2 for \gamma direction ----
c
         do 90 i=1, m
c
            do 70 j=1, n
               temp1(j)=zp(i,j,1)
   70       continue
c
            call derlim2(y,temp1,n,derivp,derivk)
            call spline(y,temp1,temp2,n,derivp,derivk)
c
            do 80 j=1, n
               zp(i,j,3)=temp2(j)
   80       continue
c
   90    continue
      return
      end
c
c
c***********************************************************
c
c    To cal. the first-order derivative
c
c    derivp: df/dx at x(1)
c    derivk: df/dx at x(N)
c
c***********************************************************
c
      subroutine derlim2(x,tem,m,derivp,derivk)
c
      implicit real*8(a-h,o-z)
c
      dimension  x(0:m),tem(0:m)
c
      t0= (tem(1)-tem(0))/(x(1)-x(0))
      t1=(tem(2)-tem(1))/(x(2)-x(1))
      t2=(t1-t0)/(x(2)-x(0))
      derivp=2*t2*x(0) + (t0-t2*(x(0)+x(1)))
c
      t1= (tem(m-2)-tem(m-1))/(x(m-2)-x(m-1))
      t2=(tem(m-1)-tem(m))/(x(m-1)-x(m))
      t3=(t2-t1)/(x(m)-x(m-2))
      derivk=2*t3*x(m) + (t1-t3*(x(m-1)+x(m-2)))
c
      return
      end
c
c
c======================================================================c

      subroutine spline(x,y,y2,n,yp0,ypn)

c======================================================================c
c
c     SPLINE-Routine of "Numerical Recipies" p.88
c
c input:
c     X,Y       tabulated function (0..N)
c     YP0,YPN   first derivatives at point 0 and N
c               of the interpolating spline-function
c               (if larger then 1.e30, natural spline: y''= 0)
c output:
c     Y2        second derivatives of the interpolating function
c               is used as input for function SPLINT
c
c----------------------------------------------------------------------c
      implicit real*8(a-h,o-z)
      parameter (nmax=500)
      dimension x(0:n),y(0:n),y2(0:n),u(0:nmax)
c
      if (nmax.lt.n) stop ' in SPLINE: nmax too small'
      if (yp0.gt.0.999d30) then
         y2(0) = 0.0
         u(0)  = 0.0
      else
         y2(0) = -0.5d0
         u(0)  = (3.d0/(x(1)-x(0))) * ((y(1)-y(0))/(x(1)-x(0))-yp0)
      endif
      do 11 i = 1,n-1
         sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
         p     = sig*y2(i-1) + 2.d0
         y2(i) = (sig - 1.d0)/p
         u(i)  = (6.d0*( (y(i+1)-y(i))/(x(i+1)-x(i)) -
     &                   (y(i)-y(i-1))/(x(i)-x(i-1)) )/
     &                   (x(i+1)-x(i-1)) - sig*u(i-1))/p
   11    continue
      if (ypn.gt..999d30) then
         qn = 0.0
         un = 0.0
      else
         qn = 0.5d0
         un = (3.d0/(x(n)-x(n-1))) * (ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do 12 k = n-1,0,-1
         y2(k) = y2(k)*y2(k+1)+u(k)
   12 continue
c
      return
c-end-SPLINE
      end
c======================================================================c

      subroutine splint(is,xa,ya,y2a,n,x,y,y1,y2)

c======================================================================c
c
c     SPLINT-Routine of "Numerical Recipies" p.89
c
c input:
c     XA,YA     tabulated function (0:N)
c     Y2A       first derivatives (output von SPLINE)
c     X         given value on the abscissa
c
c output:
c   is = 0:  Y  interpolated value  Y(x)
c   is = 1;  y1 in addition interpolated value of the derivativ dY/dx
c   is = 2;  y2 in addition interpolated value of 2nd derivativ d2Y/dx2
c
c----------------------------------------------------------------------c
      implicit real*8(a-h,o-z)
      dimension xa(0:n),ya(0:n),y2a(0:n)
      data sixth/0.1666666666666666667d0/
c
      klo = 0
      khi = n
    1 if (khi-klo.gt.1) then
         k = (khi+klo)/2
         if (xa(k).gt.x) then
            khi = k
         else
            klo = k
         endif
         goto 1
      endif
      h = xa(khi)-xa(klo)
      if (h.eq.0.0) stop ' in SPINT: bad xa input '
      hi = 1.d0/h
      a  = (xa(khi)-x)*hi
      b  = (x-xa(klo))*hi
c
c     value of the function
      y = a*ya(klo)+b*ya(khi)+
     &    ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)*sixth
c
c     first derivative
      if (is.lt.1) return
      y1 = hi*(-ya(klo)+ya(khi)) +
     &    (-(3*a**2-1)*y2a(klo)+(3*b**2-1)*y2a(khi))*h*sixth
c
c     second derivative
      if (is.lt.2) return
      y2 = a*y2a(klo) + b*y2a(khi)
c
      return
c-end-SPLINT
      end
c*************************************************************
c
c    To calculate the value and first derivative of
c    of interpolation at v
c
c*************************************************************
c
      double precision function fu(v,alo,ahi,blo,bhi,slo,shi)
      double precision v,alo,ahi,blo,bhi,slo,shi
      double precision h,a,b,sk1,sk2
      h = ahi-alo
      if (h.eq.0.) then
       write(*,*) 'pause in routine fu'
       write(*,*) ' ... bad XA input'
      endif
      a = (ahi-v)/h
      b = (v-alo)/h
      sk1=a*blo+b*bhi
      sk2=((a*a*a-a)*slo+(b*b*b-b)*shi)*(h*h)/6.0
      fu=sk1+sk2
      return
      end
c
c-----------------------------------------------------------
c
      double precision function difu(v,alo,ahi,blo,bhi,slo,shi)
      double precision v,alo,ahi,blo,bhi,slo,shi
      double precision h,a,b,sk1,sk2
      h = ahi-alo
      if (h.eq.0.) then
       write(*,*) 'pause in routine difu'
       write(*,*) ' ... bad XA input'
      endif
      a = (ahi-v)/h
      b = (v-alo)/h
      sk1=(bhi-blo)/h
      sk2=(-(3*a*a-1)*slo+(3*b*b-1)*shi)*h/6.0
      difu=sk1+sk2
      return
      end
c

c======================================================================c

      subroutine start(lpr)
c======================================================================c
c
c     initializes potentials (inin=1) and pairing tensor (inink=1)
c----------------------------------------------------------------------c
      use dirqfampar;

      implicit real*8 (a-h,o-z)
!ab   include 'dirqfam.par'

      logical lpr
      common /initia/ inin,inink

      if(inin.eq.1) call startpot(lpr)
      if(inink.eq.1) call startdel(lpr)

      return
      end

c======================================================================c
      subroutine startpot(lpr)
c======================================================================c
      use dirqfampar;
      use coulmb;
      use gaussh;
      use gaussl;
      use gfvmod;
      use potpot;
      use deldel;

      implicit real*8 (a-h,o-z)
!ab   include 'dirqfam.par'
c
      logical lpr
c
      character nucnam*2                                 ! common nucnuc
c
      dimension vso(2),r0v(2),av(2),rso(2),aso(2)
      dimension rrv(2),rls(2),vp(2),vls(2)
c
      common /baspar/ hom,hb0,b0
!ab   common /coulmb/ cou(MG),drvp(MG)
      common /defbas/ beta0,q,bp,bz
      common /defini/ betai
!ab   common /gaussh/ xh(0:NGH),wh(0:NGH),zb(0:NGH)
!ab   common /gaussl/ xl(0:NGL),wl(0:NGL),sxl(0:NGL),rb(0:NGL)
!ab   common /gfvsq / sq(0:IGFV)
      common /initia/ inin,inink
      common /mathco/ zero,one,two,half,third,pi
      common /nucnuc/ amas,npr(3),nucnam
      common /optopt/ itx,icm,icou,ipc,inl,idd
      common /physco/ hbc,alphi,r0
!ab   common /potpot/ vps(MG,2),vms(MG,2)
!ab   common /deldel/ de(NHHX,NB2X)
      common /pair  / del(2),spk(2),spk0(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp
c
c=======================================================================
c     Saxon-Woods parameter von Koepf und Ring, Z.Phys. (1991)
      data v0/-71.28/,akv/0.4616/
      data r0v/1.2334,1.2496/,av/0.615,0.6124/
      data vso/11.1175,8.9698/
      data rso/1.1443,1.1401/,aso/0.6476,0.6469/

c
      if (lpr) then
      write(l6,*) ' ***** BEGIN START ******************************'
      endif

c---- saxon-woods potential
      if (lpr) then
         write(l6,100) ' v0     = ',v0
         write(l6,100) ' kappa  = ',akv
         write(l6,100) ' lambda = ',vso
         write(l6,100) ' r0     = ',r0v
         write(l6,100) ' a      = ',av
         write(l6,100) ' r0-so  = ',rso
         write(l6,100) ' a-so   = ',aso
         write(l6,100) ' betai  = ',betai
      endif
c
      betas = betai * half* dsqrt(5.d0/(4.d0*pi))
      fac =  one + betas
      fac = (one + betas*cos(120.d0*pi/180.d0))*fac
      fac = (one + betas*cos(-120.d0*pi/180.d0))*fac
      fac = fac**(-third)
      do ih = 0,NGH
         z  = zb(ih)
         zz = z**2
      do il = 0,NGL
         ihl = 1+ih + il*(NGH+1)
         rr = (rb(il)**2 + zz)
         r  = sqrt(rr)
c
C------- Woods Saxon
         ctet = zz/rr
         p20  = 3*ctet - one
         facb = fac*(one + betas*p20)
         do it = 1,2
            ita = 3-it
            rrv(it) = r0v(it)*amas**third
            rls(it) = rso(it)*amas**third
            vp(it)  = v0*(one - akv*(npr(it)-npr(ita))/amas)
            vls(it) = vp(it) * vso(it)
c
            argv=(r - rrv(it)*facb) / av(it)
            if (argv.lt.65.d0) then
                u = vp(it) /(one + exp(argv))
            else
                u = 0.d0
            endif
            argo=(r - rls(it)*facb) / aso(it)
            if (argo.lt.65.d0) then
                w  = -vls(it)/(one + exp(argo ))
            else
                w = 0.d0
            endif
c
            vps(ihl,it) = u
            vms(ihl,it) = w
c
         enddo   ! it
c
c------- Coulomb
         cou(ihl) = zero
         if (icou.ne.0) then
            rc = rrv(2)
            if (r.lt.rc) then
               c = half*(3/rc-r*r/(rc**3))
            else
               c = one/r
            endif
            cou(ihl)   = c*npr(2)/alphi
            vps(ihl,2) = vps(ihl,2) + cou(ihl)*hbc
            vms(ihl,2) = vms(ihl,2) + cou(ihl)*hbc
         endif
      enddo   ! il
      enddo   ! ih



  100 format(a,2f10.4)
c
      if (lpr) then
      write(l6,*) ' ***** END START ********************************'
      endif
c
      return
c-end STARTPOT
      end
c======================================================================c
      subroutine startdel(lpr)
c======================================================================c
      use dirqfampar;
      use blokap;
      use bloosc;
      use deldel;

      implicit real*8 (a-h,o-z)
!ab   include 'dirqfam.par'
c
      logical lpr
!ab   character tb*6

!ab   common /blokap/ nb,kb(NBX),mb(NBX),tb(NBX)
!ab   common /bloosc/ ia(NBX,2),id(NBX,2)
      common /mathco/ zero,one,two,half,third,pi
!ab   common /deldel/ de(NHHX,NB2X)
      common /pair  / del(2),spk(2),spk0(2)
      common /tapes / l6,lin,lou,lwin,lwou,lplo,laka,lvpp

      call mzero(NHHX,NHHX,NB2X,de)

      do ib = 1,nb
         nf  = id(ib,1)
         ng  = id(ib,2)
         nh = nf + ng
         m1 = ib
         m2 = ib + NBX
c
         do it = 1,2
            if (abs(del(it)).lt.1.d-5) then
               del(it)=1.d-5
               spk(it) = zero
            else
               spk(it) = one
            endif
         enddo   ! it

         ne = nf
         delt = del(1)
         do n = 1,ne
            de(n+(n-1)*nh,m1) = delt
         enddo   ! n

         ne = nf
         delt = del(2)
         do n = 1,ne
            de(n+(n-1)*nh,m2) = delt
         enddo   ! n
      enddo    ! ib

      return
c-end STARTDEL
      end

