      module dirqfampar
      implicit none;

c-----maximal number for GFV
      integer IGFV;

c-----number of r-meshpoints (4n+1 points)
      double precision RMAX;
      integer MR;

c-----number of q-meshpoints (4n+1 points)
      double precision QMAX;
      integer MQ;

c-----number of gauss-meshpoints
      integer NGH;
      integer NGL;
      integer NGLEG;

c---- maximal oscillator quantum number for fermions
      integer N0FX;
      integer nxx;

c-----maximal number of (k,parity)-blocks
      integer NBX;

c-----max. number of all levels for protons or neutrons
      integer NTX;

c-----max. number of eigenstates for protons or neutrons
      integer KX;

c-----max. nz-quantum number of fermions
      integer NZX

c-----max. nr-quantum number of fermions
      integer NRX;

c-----max. ml-quantum number of fermions
      integer MLX;

c-----maximal dimension F of one k-block
      integer NFX;

c-----maximal dimension G of one k-block
      integer NGX;

c-----maximum of nf and ng in all blocks
      integer NDX;

c-----oscillator quantum number for bosons
      integer N0BX;

c-----number of bosons
      integer nbxx;
      integer NOX;

c-----for the plot
      integer NOX1;

c-----auxiliary
      integer NGH2;
      integer NB2X;
      integer NHX;
      integer NDDX;
      integer NHBX;
      integer NHBQX;
      integer NFFX;
      integer NFGX;
      integer NHHX;
      integer MG;
      integer N02;
      integer NNNX;

c-----working space
      integer MVX;
      integer MVTX;

c-----QFAM parameters
      integer NBSX;
      integer J_MAX;
      integer NFAM_BROYD;
      integer MFAM_BROYD;
      integer NWMAX;
      integer KTRUNC;
      integer NCOORD;
      integer NMESMAX;
      integer NHMAX;

      end module dirqfampar
