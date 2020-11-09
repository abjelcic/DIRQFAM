      module simplex
      implicit none;

      integer N_total;
      integer N_blocks;
      integer,   dimension(:),   allocatable :: ia_spx;
      integer,   dimension(:),   allocatable :: id_spx;
      integer,   dimension(:),   allocatable :: nf_size;
      integer,   dimension(:),   allocatable :: ng_size;
      integer,   dimension(:,:), allocatable :: nz_spx;
      integer,   dimension(:,:), allocatable :: nr_spx;
      integer,   dimension(:,:), allocatable :: ml_spx;
      character, dimension(:,:), allocatable :: fg_spx;

      contains
          subroutine alloc_simplex()
          use dirqfampar
          implicit none;

          allocate(      ia_spx(NBX) );
          allocate(      id_spx(NBX) );
          allocate(     nf_size(NBX) );
          allocate(     ng_size(NBX) );
          allocate( nz_spx(NBSX,NBX) );
          allocate( nr_spx(NBSX,NBX) );
          allocate( ml_spx(NBSX,NBX) );
          allocate( fg_spx(NBSX,NBX) );

          end

          subroutine dealloc_simplex()
          implicit none;

          deallocate( ia_spx  );
          deallocate( id_spx  );
          deallocate( nf_size );
          deallocate( ng_size );
          deallocate( nz_spx  );
          deallocate( nr_spx  );
          deallocate( ml_spx  );
          deallocate( fg_spx  );

          end

      end module simplex

      module v_matrix
      implicit none;

      double complex, dimension(:,:,:,:), allocatable :: v;

      contains
          subroutine alloc_v_matrix()
          use dirqfampar
          implicit none;

          allocate( v(NBSX,NBSX,NBX,2) );

          end

          subroutine dealloc_v_matrix()
          implicit none;

          deallocate( v );

          end

      end module v_matrix

      module v_energy
      implicit none;

      double precision, dimension(:,:,:), allocatable :: E_fam_v;

      contains
          subroutine alloc_v_energy()
          use dirqfampar
          implicit none;

          allocate( E_fam_v(NBSX,NBX,2) );

          end

          subroutine dealloc_v_energy()
          implicit none;

          deallocate( E_fam_v );

          end

      end module v_energy

      module u_matrix
      implicit none;

      double complex, dimension(:,:,:,:), allocatable :: u;

      contains
          subroutine alloc_u_matrix()
          use dirqfampar
          implicit none;

          allocate( u(NBSX,NBSX,NBX,2) );

          end

          subroutine dealloc_u_matrix()
          implicit none;

          deallocate( u );

          end

      end module u_matrix

      module u_energy
      implicit none;

      double precision, dimension(:,:,:), allocatable :: E_fam_u;

      contains
          subroutine alloc_u_energy()
          use dirqfampar
          implicit none;

          allocate( E_fam_u(NBSX,NBX,2) );

          end

          subroutine dealloc_u_energy()
          implicit none;

          deallocate( E_fam_u );

          end

      end module u_energy

      module fam
      implicit none;

      double precision omega;
      double precision gamma_smear;

      double precision omega_start;
      double precision omega_end;
      double precision delta_omega;
      double precision omega_print;

      integer          i_calculation_type;
      integer          i_coulomb;
      integer          i_pairing;

      integer          J_multipole;
      integer          K_multipole;
      integer          ISO;

      integer          tape_strength
      integer          tape_rhov;
      integer          tape_nuclocfunc;

      contains
          subroutine alloc_fam()
          use dirqfampar
          implicit none;

          end

          subroutine dealloc_fam()
          implicit none;

          end

      end module fam

      module ddpc1ddme2
      implicit none;

      double precision a_s , b_s , c_s , d_s;
      double precision a_v , b_v , c_v , d_v;
      double precision a_tv, b_tv, c_tv, d_tv;
      double precision del_s;

      double precision a_sig, b_sig, c_sig, d_sig, g0_sig, m_sig;
      double precision a_ome, b_ome, c_ome, d_ome, g0_ome, m_ome;
      double precision a_rho,                      g0_rho, m_rho;
      double precision rho_sat;

      contains
          subroutine alloc_ddpc1ddme2()
          use dirqfampar
          implicit none;

          end

          subroutine dealloc_ddpc1ddme2()
          implicit none;

          end

      end module ddpc1ddme2

      module pairparams
      implicit none;

      double precision G_pairing, a_pairing;

      contains
          subroutine alloc_pairparams()
          use dirqfampar
          implicit none;

          end

          subroutine dealloc_pairparams()
          implicit none;

          end

      end module pairparams

      module fam_energies
      implicit none;

      double precision, dimension(:,:), allocatable :: E_fam;

      contains
          subroutine alloc_fam_energies()
          use dirqfampar
          implicit none;

          allocate( E_fam(NTX,2) );

          end

          subroutine dealloc_fam_energies()
          implicit none;

          deallocate( E_fam );

          end

      end module fam_energies

      module nnz_blocks
      implicit none;

      logical, dimension(:,:), allocatable :: dh_nnz;
      logical, dimension(:,:), allocatable :: dDelta_nnz;
      logical, dimension(:,:), allocatable :: dkappa_nnz;
      logical, dimension(:,:), allocatable :: f_nnz;

      contains
          subroutine alloc_nnz_blocks()
          use dirqfampar
          implicit none;

          allocate(     dh_nnz(NBX,NBX) );
          allocate( dDelta_nnz(NBX,NBX) );
          allocate( dkappa_nnz(NBX,NBX) );
          allocate(      f_nnz(NBX,NBX) );

          end

          subroutine dealloc_nnz_blocks()
          implicit none;

          deallocate( dh_nnz     );
          deallocate( dDelta_nnz );
          deallocate( dkappa_nnz );
          deallocate( f_nnz      );

          end

      end module nnz_blocks

      module quadrature
      implicit none;

      double precision, dimension(:),   allocatable :: zb_fam;
      double precision, dimension(:),   allocatable :: wz;
      double precision, dimension(:),   allocatable :: rb_fam;
      double precision, dimension(:),   allocatable :: wr;
      double precision, dimension(:,:), allocatable :: wzwr;

      contains
          subroutine alloc_quadrature()
          use dirqfampar
          implicit none;

          allocate( zb_fam(NGH)     );
          allocate(     wz(NGH)     );
          allocate( rb_fam(NGL)     );
          allocate(     wr(NGL)     );
          allocate(   wzwr(NGH,NGL) );

          end

          subroutine dealloc_quadrature()
          implicit none;

          deallocate( zb_fam );
          deallocate( wz     );
          deallocate( rb_fam );
          deallocate( wr     );
          deallocate( wzwr   );

          end

      end module quadrature

      module basis
      implicit none;

      double precision, dimension(:,:), allocatable ::  phi_z;
      double precision, dimension(:,:), allocatable :: dphi_z;
      double precision, dimension(:,:), allocatable ::  phi_r;
      double precision, dimension(:,:), allocatable :: dphi_r;

      contains
          subroutine alloc_basis()
          use dirqfampar
          implicit none;

          allocate(  phi_z(-NGH:NGH,NTX) );
          allocate( dphi_z(-NGH:NGH,NTX) );
          allocate(  phi_r(   1:NGL,NTX) );
          allocate( dphi_r(   1:NGL,NTX) );

          end

          subroutine dealloc_basis()
          implicit none;

          deallocate(  phi_z );
          deallocate( dphi_z );
          deallocate(  phi_r );
          deallocate( dphi_r );

          end

      end module basis

      module wbasis
      implicit none;

      double precision, dimension(:,:,:), allocatable :: wPhi;

      contains
          subroutine alloc_wbasis()
          use dirqfampar
          implicit none;

          allocate( wPhi(NGH,NGL,NTX) );

          end

          subroutine dealloc_wbasis()
          implicit none;

          deallocate( wPhi );

          end

      end module wbasis

      module PHI
      implicit none;

      integer k_PHI;
      double precision, dimension(:,:), allocatable :: PHI_U;
      double precision, dimension(:,:), allocatable :: PHI_SVt;

      contains
          subroutine alloc_PHI()
          use dirqfampar
          implicit none;

          allocate( PHI_U(NTX,KTRUNC)      );
          allocate( PHI_SVt(KTRUNC,NCOORD) );

          end

          subroutine dealloc_PHI()
          implicit none;

          deallocate( PHI_U   );
          deallocate( PHI_SVt );

          end

      end module PHI

      module gs_dens
      implicit none;

      double precision, dimension(:,:,:), allocatable :: rhov_GS;
      double precision, dimension(:,:),   allocatable :: rhos_GS;

      contains
          subroutine alloc_gs_dens()
          use dirqfampar
          implicit none;

          allocate( rhov_GS(-NGH:NGH,1:NGL,2) );
          allocate( rhos_GS(-NGH:NGH,1:NGL  ) );

          end

          subroutine dealloc_gs_dens()
          implicit none;

          deallocate( rhov_GS );
          deallocate( rhos_GS );

          end

      end module gs_dens

      module gs_mesons
      implicit none;

      double precision, dimension(:,:), allocatable :: sig_GS;
      double precision, dimension(:,:), allocatable :: ome_GS;
      double precision, dimension(:,:), allocatable :: rho_GS;

      contains
          subroutine alloc_gs_mesons()
          use dirqfampar
          implicit none;

          allocate( sig_GS(-NGH:NGH,1:NGL) );
          allocate( ome_GS(-NGH:NGH,1:NGL) );
          allocate( rho_GS(-NGH:NGH,1:NGL) );

          end

          subroutine dealloc_gs_mesons()
          implicit none;

          deallocate( sig_GS );
          deallocate( ome_GS );
          deallocate( rho_GS );

          end

      end module gs_mesons

      module fmatrix
      implicit none;

      double complex, dimension(:,:,:), allocatable :: f1_JK;
      double complex, dimension(:,:,:), allocatable :: f2_JK;

      contains
          subroutine alloc_fmatrix()
          use dirqfampar
          implicit none;

          allocate( f1_JK(NTX,NTX,2) );
          allocate( f2_JK(NTX,NTX,2) );

          end

          subroutine dealloc_fmatrix()
          implicit none;

          deallocate( f1_JK );
          deallocate( f2_JK );

          end

      end module fmatrix

      module f02f20matrix
      implicit none;

      double complex, dimension(:,:,:), allocatable :: f20;
      double complex, dimension(:,:,:), allocatable :: f02;

      contains
          subroutine alloc_f02f20matrix()
          use dirqfampar
          implicit none;

          allocate( f20(NTX,NTX,2) );
          allocate( f02(NTX,NTX,2) );

          end

          subroutine dealloc_f02f20matrix()
          implicit none;

          deallocate( f20 );
          deallocate( f02 );

          end

      end module f02f20matrix

      module spurious
      implicit none;

      double complex lamR;
      double complex lamP;
      double complex, dimension(:),     allocatable:: RcmPcm_commutator;
      double complex, dimension(:,:,:), allocatable:: r20;
      double complex, dimension(:,:,:), allocatable:: p20;

      contains
          subroutine alloc_spurious()
          use dirqfampar
          implicit none;

          allocate( RcmPcm_commutator(2) );
          allocate(       r20(NTX,NTX,2) );
          allocate(       p20(NTX,NTX,2) );

          end

          subroutine dealloc_spurious()
          implicit none;

          deallocate( RcmPcm_commutator );
          deallocate( r20               );
          deallocate( p20               );

          end

      end module spurious

      module mesmat
      implicit none;

      integer,          dimension(:),     allocatable :: nP;
      double precision, dimension(:,:,:), allocatable :: Psig;
      double precision, dimension(:,:,:), allocatable :: Pome;
      double precision, dimension(:,:,:), allocatable :: Prho;

      contains
          subroutine alloc_mesmat()
          use dirqfampar
          implicit none;

          allocate(                nP(0:J_MAX+1) );
          allocate( Psig(NHMAX,NCOORD,0:J_MAX+1) );
          allocate( Pome(NHMAX,NCOORD,0:J_MAX+1) );
          allocate( Prho(NHMAX,NCOORD,0:J_MAX+1) );

          end

          subroutine dealloc_mesmat()
          implicit none;

          deallocate( nP   );
          deallocate( Psig );
          deallocate( Pome );
          deallocate( Prho );

          end

      end module mesmat

      module fam_green
      implicit none;

      double precision, dimension(:,:,:,:), allocatable :: G1;
      double precision, dimension(:,:,:,:), allocatable :: G2;
      double precision, dimension(:,:,:,:), allocatable :: G3;

      contains
          subroutine alloc_fam_green()
          use dirqfampar
          implicit none;

          allocate( G1(-NGH:NGH,1:NGL,-NGH:NGH,1:NGL) );
          allocate( G2(-NGH:NGH,1:NGL,-NGH:NGH,1:NGL) );
          allocate( G3(-NGH:NGH,1:NGL,-NGH:NGH,1:NGL) );

          end

          subroutine dealloc_fam_green()
          implicit none;

          deallocate( G1 );
          deallocate( G2 );
          deallocate( G3 );

          end

      end module fam_green

      module Wpairing
      implicit none;

      double precision, dimension(:), allocatable :: W;

      contains
          subroutine alloc_Wpairing()
          use dirqfampar
          implicit none;

          allocate( W( NWMAX ) );

          end

          subroutine dealloc_Wpairing()
          implicit none;

          deallocate( W );

          end

      end module Wpairing

      module dh
      implicit none;

      double complex, dimension(:,:,:), allocatable :: dh_1;
      double complex, dimension(:,:,:), allocatable :: dh_2;

      contains
          subroutine alloc_dh()
          use dirqfampar
          implicit none;

          allocate( dh_1(NTX,NTX,2) );
          allocate( dh_2(NTX,NTX,2) );

          end

          subroutine dealloc_dh()
          implicit none;

          deallocate( dh_1 );
          deallocate( dh_2 );

          end

      end module dh

      module dDelta
      implicit none;

      double complex, dimension(:,:,:), allocatable :: dDelta_pl;
      double complex, dimension(:,:,:), allocatable :: dDelta_mi;

      contains
          subroutine alloc_dDelta()
          use dirqfampar
          implicit none;

          allocate( dDelta_pl(NTX,NTX,2) );
          allocate( dDelta_mi(NTX,NTX,2) );

          end

          subroutine dealloc_dDelta()
          implicit none;

          deallocate( dDelta_pl );
          deallocate( dDelta_mi );

          end

      end module dDelta

      module ddens
      implicit none;

      double complex, dimension(:,:,:), allocatable :: drho_v;
      double complex, dimension(:,:),   allocatable :: drho_s;

      contains
          subroutine alloc_ddens()
          use dirqfampar
          implicit none;

          allocate( drho_v(-NGH:NGH,1:NGL,2) );
          allocate( drho_s(-NGH:NGH,1:NGL  ) );

          end

          subroutine dealloc_ddens()
          implicit none;

          deallocate( drho_v );
          deallocate( drho_s );

          end

      end module ddens

      module fam_iter
      implicit none;

      integer          iter;
      integer          iter_max;
      double precision error;
      double precision tol;

      contains
          subroutine alloc_fam_iter()
          use dirqfampar
          implicit none;

          end

          subroutine dealloc_fam_iter()
          implicit none;

          end

      end module fam_iter

      module h02h20matrix
      implicit none;

      double complex, dimension(:,:,:), allocatable :: h20;
      double complex, dimension(:,:,:), allocatable :: h02;

      contains
          subroutine alloc_h02h20matrix()
          use dirqfampar
          implicit none;

          allocate( h20(NTX,NTX,2) );
          allocate( h02(NTX,NTX,2) );

          end

          subroutine dealloc_h02h20matrix()
          implicit none;

          deallocate( h20 );
          deallocate( h02 );

          end

      end module h02h20matrix

      module xyfam
      implicit none;

      double complex, dimension(:,:,:), allocatable :: x_fam;
      double complex, dimension(:,:,:), allocatable :: y_fam;

      contains
          subroutine alloc_xyfam()
          use dirqfampar
          implicit none;

          allocate( x_fam(NTX,NTX,2) );
          allocate( y_fam(NTX,NTX,2) );

          end

          subroutine dealloc_xyfam()
          implicit none;

          deallocate( x_fam );
          deallocate( y_fam );

          end

      end module xyfam

      module drho
      implicit none;

      double complex, dimension(:,:,:), allocatable :: drho_1;
      double complex, dimension(:,:,:), allocatable :: drho_2;

      contains
          subroutine alloc_drho()
          use dirqfampar
          implicit none;

          allocate( drho_1(NTX,NTX,2) );
          allocate( drho_2(NTX,NTX,2) );

          end

          subroutine dealloc_drho()
          implicit none;

          deallocate( drho_1 );
          deallocate( drho_2 );

          end

      end module drho

      module dkappa
      implicit none;

      double complex, dimension(:,:,:), allocatable :: dkappa_pl;
      double complex, dimension(:,:,:), allocatable :: dkappa_mi;

      contains
          subroutine alloc_dkappa()
          use dirqfampar
          implicit none;

          allocate( dkappa_pl(NTX,NTX,2) );
          allocate( dkappa_mi(NTX,NTX,2) );

          end

          subroutine dealloc_dkappa()
          implicit none;

          deallocate( dkappa_pl );
          deallocate( dkappa_mi );

          end

      end module dkappa

      module dcurr
      implicit none;

      double complex, dimension(:,:,:), allocatable :: dj_r;
      double complex, dimension(:,:,:), allocatable :: dj_p;
      double complex, dimension(:,:,:), allocatable :: dj_z;
      double complex, dimension(:,:,:), allocatable :: dj_1;
      double complex, dimension(:,:,:), allocatable :: dj_2;
      double complex, dimension(:,:,:), allocatable :: dj_3;

      contains
          subroutine alloc_dcurr()
          use dirqfampar
          implicit none;

          allocate( dj_r(-NGH:NGH,1:NGL,2) );
          allocate( dj_p(-NGH:NGH,1:NGL,2) );
          allocate( dj_z(-NGH:NGH,1:NGL,2) );
          allocate( dj_1(-NGH:NGH,1:NGL,2) );
          allocate( dj_2(-NGH:NGH,1:NGL,2) );
          allocate( dj_3(-NGH:NGH,1:NGL,2) );

          end

          subroutine dealloc_dcurr()
          implicit none;

          deallocate( dj_r );
          deallocate( dj_p );
          deallocate( dj_z );
          deallocate( dj_1 );
          deallocate( dj_2 );
          deallocate( dj_3 );

          end

      end module dcurr

      module dlaplace
      implicit none;

      double complex, dimension(:,:), allocatable :: ldrho_vp;
      double complex, dimension(:,:), allocatable :: ldrho_s;
      double complex, dimension(:,:), allocatable :: ldj_1p;
      double complex, dimension(:,:), allocatable :: ldj_2p;
      double complex, dimension(:,:), allocatable :: ldj_3p;

      contains
          subroutine alloc_dlaplace()
          use dirqfampar
          implicit none;

          allocate( ldrho_vp(-NGH:NGH,1:NGL) );
          allocate(  ldrho_s(-NGH:NGH,1:NGL) );
          allocate(   ldj_1p(-NGH:NGH,1:NGL) );
          allocate(   ldj_2p(-NGH:NGH,1:NGL) );
          allocate(   ldj_3p(-NGH:NGH,1:NGL) );

          end

          subroutine dealloc_dlaplace()
          implicit none;

          deallocate( ldrho_vp );
          deallocate( ldrho_s  );
          deallocate( ldj_1p   );
          deallocate( ldj_2p   );
          deallocate( ldj_3p   );

          end

      end module dlaplace

      module famlaplacianmod
      implicit none;

      integer NSH;
      integer NSIZE;
      integer,          dimension(:),     allocatable :: NZZ;
      integer,          dimension(:),     allocatable :: NRR;
      double complex,   dimension(:),     allocatable :: c;
      double precision, dimension(:,:),   allocatable :: phiz;
      double precision, dimension(:,:,:), allocatable :: phirK;

      contains
          subroutine alloc_famlaplacianmod()
          use dirqfampar
          implicit none;

          NSH   = 2*(N0FX+1);
          NSIZE = ( (NSH+1)*(NSH+3) + 1 - MOD(NSH,2) )/4;
          allocate( NZZ(NSIZE) );
          allocate( NRR(NSIZE) );
          allocate(   c(NSIZE) );
          allocate( phiz (-NGH:NGH,0:(NSH)             ) );
          allocate( phirK(   1:NGL,0:(NSH/2),0:J_MAX+1 ) );

          end

          subroutine dealloc_famlaplacianmod()
          implicit none;

          deallocate( NZZ   );
          deallocate( NRR   );
          deallocate( c     );
          deallocate( phiz  );
          deallocate( phirK );

          end

      end module famlaplacianmod

      module dmesons
      implicit none;

      double complex, dimension(:,:), allocatable :: dsig;
      double complex, dimension(:,:), allocatable :: dome_0;
      double complex, dimension(:,:), allocatable :: dome_r;
      double complex, dimension(:,:), allocatable :: dome_p;
      double complex, dimension(:,:), allocatable :: dome_z;
      double complex, dimension(:,:), allocatable :: drho_0;
      double complex, dimension(:,:), allocatable :: drho_r;
      double complex, dimension(:,:), allocatable :: drho_p;
      double complex, dimension(:,:), allocatable :: drho_z;

      contains
          subroutine alloc_dmesons()
          use dirqfampar
          implicit none;

          allocate(   dsig(-NGH:NGH,1:NGL) );
          allocate( dome_0(-NGH:NGH,1:NGL) );
          allocate( dome_r(-NGH:NGH,1:NGL) );
          allocate( dome_p(-NGH:NGH,1:NGL) );
          allocate( dome_z(-NGH:NGH,1:NGL) );
          allocate( drho_0(-NGH:NGH,1:NGL) );
          allocate( drho_r(-NGH:NGH,1:NGL) );
          allocate( drho_p(-NGH:NGH,1:NGL) );
          allocate( drho_z(-NGH:NGH,1:NGL) );

          end

          subroutine dealloc_dmesons()
          implicit none;

          deallocate( dsig   );
          deallocate( dome_0 );
          deallocate( dome_r );
          deallocate( dome_p );
          deallocate( dome_z );
          deallocate( drho_0 );
          deallocate( drho_r );
          deallocate( drho_p );
          deallocate( drho_z );

          end

      end module dmesons

      module dcoul
      implicit none;

      double complex, dimension(:,:), allocatable :: dVCou_0;
      double complex, dimension(:,:), allocatable :: dVCou_r;
      double complex, dimension(:,:), allocatable :: dVCou_p;
      double complex, dimension(:,:), allocatable :: dVCou_z;

      contains
          subroutine alloc_dcoul()
          use dirqfampar
          implicit none;

          allocate( dVCou_0(-NGH:NGH,1:NGL) );
          allocate( dVCou_r(-NGH:NGH,1:NGL) );
          allocate( dVCou_p(-NGH:NGH,1:NGL) );
          allocate( dVCou_z(-NGH:NGH,1:NGL) );

          end

          subroutine dealloc_dcoul()
          implicit none;

          deallocate( dVCou_0 );
          deallocate( dVCou_r );
          deallocate( dVCou_p );
          deallocate( dVCou_z );

          end

      end module dcoul

      module dpot
      implicit none;

      double complex, dimension(:,:,:), allocatable :: dVpS;
      double complex, dimension(:,:,:), allocatable :: dVmS;
      double complex, dimension(:,:,:), allocatable :: dSig_z;
      double complex, dimension(:,:,:), allocatable :: dSig_r;
      double complex, dimension(:,:,:), allocatable :: dSig_p;

      contains
          subroutine alloc_dpot()
          use dirqfampar
          implicit none;

          allocate(   dVpS(-NGH:NGH,1:NGL,2) );
          allocate(   dVmS(-NGH:NGH,1:NGL,2) );
          allocate( dSig_z(-NGH:NGH,1:NGL,2) );
          allocate( dSig_r(-NGH:NGH,1:NGL,2) );
          allocate( dSig_p(-NGH:NGH,1:NGL,2) );

          end

          subroutine dealloc_dpot()
          implicit none;

          deallocate( dVpS   );
          deallocate( dVmS   );
          deallocate( dSig_z );
          deallocate( dSig_r );
          deallocate( dSig_p );

          end

      end module dpot

      module fambroydenmod
      implicit none;

      integer nn;
      integer mm;
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
          subroutine alloc_fambroydenmod()
          use dirqfampar
          implicit none;

          nn = NFAM_BROYD;
          mm = MFAM_BROYD;

          allocate( bbeta(mm,mm) );
          allocate(    df(nn,mm) );
          allocate(    dv(nn,mm) );
          allocate(    bwork(mm) );
          allocate(     curv(nn) );
          allocate(      vin(nn) );
          allocate(      vou(nn) );
          allocate(   ibwork(mm) );

          end

          subroutine dealloc_fambroydenmod()
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

      end module fambroydenmod

      module nuclearlocfunc
      implicit none;

      double precision, dimension(:,:,:), allocatable :: C0;
      double complex,   dimension(:,:,:), allocatable :: dC;

      contains
          subroutine alloc_nuclearlocfunc()
          use dirqfampar
          implicit none;

          allocate( C0(-NGH:NGH,1:NGL,2) );
          allocate( dC(-NGH:NGH,1:NGL,2) );

          end

          subroutine dealloc_nuclearlocfunc()
          implicit none;

          deallocate( C0 );
          deallocate( dC );

          end

      end module nuclearlocfunc
