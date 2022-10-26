      module dirhbpar
      implicit none;

!-----maximal number for GFV
      integer :: IGFV;

!-----number of r-meshpoints (4n+1 points)
      double precision :: RMAX;
      integer :: MR;

!-----number of q-meshpoints (4n+1 points)
      double precision :: QMAX;
      integer :: MQ;

!-----number of gauss-meshpoints
      integer :: NGH;
      integer :: NGL;
      integer :: NGLEG;

!---- maximal oscillator quantum number for fermions
      integer :: N0FX;
      integer :: nxx;

!-----maximal number of (k,parity)-blocks
      integer :: NBX;

!-----max. number of all levels for protons or neutrons
      integer :: NTX;

!-----max. number of eigenstates for protons or neutrons
      integer :: KX;

!-----max. nz-quantum number of fermions
      integer :: NZX

!-----max. nr-quantum number of fermions
      integer :: NRX;

!-----max. ml-quantum number of fermions
      integer :: MLX;

!-----maximal dimension F of one k-block
      integer :: NFX;

!-----maximal dimension G of one k-block
      integer :: NGX;

!-----maximum of nf and ng in all blocks
      integer :: NDX;

!-----oscillator quantum number for bosons
      integer :: N0BX;

!-----number of bosons
      integer :: nbxx;
      integer :: NOX;

!-----for the plot
      integer :: NOX1;

!-----auxiliary
      integer :: NGH2;
      integer :: NB2X;
      integer :: NHX;
      integer :: NDDX;
      integer :: NHBX;
      integer :: NHBQX;
      integer :: NFFX;
      integer :: NFGX;
      integer :: NHHX;
      integer :: MG;
      integer :: N02;
      integer :: NNNX;

!-----working space
      integer :: MVX;
      integer :: MVTX;

      end module dirhbpar
