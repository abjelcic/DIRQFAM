      module dataTypes
          implicit none;

	      type , public :: integer1Darray
	          integer , dimension(:) , allocatable :: index;
	      end type integer1Darray

	      type , public :: character1Darray
	          character , dimension(:) , allocatable :: index;
	      end type character1Darray

          type , public :: real1Darray
              double precision , dimension(:) , allocatable :: vec;
          end type real1Darray

          type , public :: real2Darray
              double precision , dimension(:,:) , allocatable :: mat;
          end type real2Darray

          type , public :: real3Darray
              double precision , dimension(:,:,:) , allocatable :: arr;
          end type real3Darray

          type , public :: complex2Darray
              double complex , dimension(:,:) , allocatable :: mat;
          end type complex2Darray


          type , public :: realBlockVector
              type(real1Darray) , dimension(:) , allocatable :: blocks;
          end type realBlockVector

          type , public :: realBlockMatrix
              type(real2Darray) , dimension(:,:) , allocatable :: blocks;
              logical           , dimension(:,:) , allocatable :: nnzblocks;
          end type realBlockMatrix

          type , public :: complexBlockMatrix
              type(complex2Darray) , dimension(:,:) , allocatable :: blocks;
              logical              , dimension(:,:) , allocatable :: nnzblocks;
          end type complexBlockMatrix

      end module dataTypes

      module fam_input
          implicit none;


          integer                        :: n0f;
          integer                        :: n0b;

          double precision               :: beta0;

          character(:)     , allocatable :: nucleusName;
          integer                        :: nucleusZ;
          integer                        :: nucleusN;

          character(:)     , allocatable :: LagrangianModelName;

		  integer                        :: calculation_type;

		  integer                        :: include_coulomb;
		  integer                        :: include_pairing;
          integer                        :: NGH;
          integer                        :: NGL;
		  double precision               :: gamma_smear;
          double precision               :: selfConsistencyTolerance;
          integer                        :: NoArnoldiVectors;

          integer          , parameter   :: J_MAX = 5;
		  integer                        :: J_multipole;
		  integer                        :: K_multipole;
		  integer                        :: Isospin;

		  double precision               :: omega_start;
		  double precision               :: omega_end;
		  double precision               :: delta_omega;

		  double precision               :: omega_print;

		  double precision               :: omega_center;
		  double precision               :: omega_radius;
		  integer                        :: NoContourPoints;

      end module fam_input

      module simplex
          use dataTypes;
          implicit none;

		  integer 							     			  :: N_total;
		  integer 								 			  :: N_blocks;
		  integer 			     , dimension(:) , allocatable :: ia_spx;
		  integer 			     , dimension(:) , allocatable :: id_spx;
		  integer 			     , dimension(:) , allocatable :: nf_size;
		  integer 			     , dimension(:) , allocatable :: ng_size;
		  type(integer1Darray)   , dimension(:) , allocatable :: nz_spx;
		  type(integer1Darray)   , dimension(:) , allocatable :: nr_spx;
		  type(integer1Darray)   , dimension(:) , allocatable :: ml_spx;
		  type(character1Darray) , dimension(:) , allocatable :: fg_spx;
		  integer                , dimension(:) , allocatable :: getBlock;
		  integer                , dimension(:) , allocatable :: getIndexWithinBlock;

      contains

		  ! Arrays of simplex module are allocated in base_simplex subroutine.

          subroutine dealloc_simplex()
              implicit none;
		      integer :: ib;

			  deallocate( ia_spx  );
			  deallocate( id_spx  );
			  deallocate( nf_size );
			  deallocate( ng_size );

			  do ib = 1 , size(nz_spx)
				  deallocate( nz_spx(ib)%index );
			  end do
			  deallocate( nz_spx );

			  do ib = 1 , size(nr_spx)
				  deallocate( nr_spx(ib)%index );
			  end do
			  deallocate( nr_spx );

			  do ib = 1 , size(ml_spx)
				  deallocate( ml_spx(ib)%index );
			  end do
			  deallocate( ml_spx );

			  do ib = 1 , size(fg_spx)
				  deallocate( fg_spx(ib)%index );
			  end do
			  deallocate( fg_spx );

			  deallocate( getBlock            );
			  deallocate( getIndexWithinBlock );

	          return;
          end subroutine dealloc_simplex

      end module simplex

      module u_matrix
          use dataTypes;
          implicit none;

          type(complexBlockMatrix) , dimension(2) :: u;

      contains

          subroutine alloc_u_matrix()
              use simplex;
              implicit none;
              integer :: it, ib1, ib2;

              do it = 1 , 2

                  allocate( u(it)%nnzblocks(1:N_blocks,1:N_blocks) );

                  u(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  do ib1 = 1 , N_blocks
                      do ib2 = 1 , N_blocks
                          if( ib1 == ib2 ) &
                              u(it)%nnzblocks(ib1,ib2) = .true.;
                      end do
                  end do


                  allocate( u(it)%blocks(1:N_blocks,1:N_blocks) );
                  do ib1 = 1 , N_blocks
                      do ib2 = 1 , N_blocks
                          if( u(it)%nnzblocks(ib1,ib2) ) &
                              allocate( u(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );
                      end do
                  end do

              end do

		      return;
          end subroutine alloc_u_matrix

          subroutine dealloc_u_matrix()
              implicit none;
              integer :: it, ib1, ib2;

              do it = 1 , 2

                  do ib1 = 1 , size(u(it)%nnzblocks,1)
                      do ib2 = 1 , size(u(it)%nnzblocks,2)
                          if( u(it)%nnzblocks(ib1,ib2) ) &
                              deallocate( u(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( u(it)%blocks );

                  deallocate( u(it)%nnzblocks );

              end do

		      return;
          end subroutine dealloc_u_matrix

      end module u_matrix

      module u_energy
          use dataTypes;
          implicit none;

          type(realBlockVector) , dimension(2) :: E_fam_u;

      contains

          subroutine alloc_u_energy()
              use simplex;
              implicit none;
              integer :: it, ib;

              do it = 1 , 2
                  allocate( E_fam_u(it)%blocks(1:N_blocks) );
                  do ib = 1 , N_blocks
                      allocate( E_fam_u(it)%blocks(ib)%vec(1:id_spx(ib)) );
                  end do
              end do

              return;
          end subroutine alloc_u_energy

          subroutine dealloc_u_energy()
              implicit none;
              integer :: it, ib;

              do it = 1 , 2
                  do ib = 1 , size(E_fam_u(it)%blocks)
                      deallocate( E_fam_u(it)%blocks(ib)%vec );
                  end do
                  deallocate( E_fam_u(it)%blocks );
              end do

              return;
          end subroutine dealloc_u_energy

      end module u_energy

      module v_matrix
          use dataTypes;
          implicit none;

          type(complexBlockMatrix) , dimension(2) :: v;

      contains

          subroutine alloc_v_matrix()
              use simplex;
              implicit none;
              integer :: it, ib1, ib2;

              do it = 1 , 2

                  allocate( v(it)%nnzblocks(1:N_blocks,1:N_blocks) );

                  v(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  do ib1 = 1 , N_blocks
                      do ib2 = 1 , N_blocks
                          if( ib1 == ib2 ) &
                              v(it)%nnzblocks(ib1,ib2) = .true.;
                      end do
                  end do


                  allocate( v(it)%blocks(1:N_blocks,1:N_blocks) );
                  do ib1 = 1 , N_blocks
                      do ib2 = 1 , N_blocks
                          if( v(it)%nnzblocks(ib1,ib2) ) &
                              allocate( v(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );
                      end do
                  end do

              end do

		      return;
          end subroutine alloc_v_matrix

          subroutine dealloc_v_matrix()
              implicit none;
              integer :: it, ib1, ib2;

              do it = 1 , 2

                  do ib1 = 1 , size(v(it)%nnzblocks,1)
                      do ib2 = 1 , size(v(it)%nnzblocks,2)
                          if( v(it)%nnzblocks(ib1,ib2) ) &
                              deallocate( v(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( v(it)%blocks );

                  deallocate( v(it)%nnzblocks );

              end do

		      return;
          end subroutine dealloc_v_matrix

      end module v_matrix

      module v_energy
          use dataTypes;
          implicit none;

          type(realBlockVector) , dimension(2) :: E_fam_v;

      contains

          subroutine alloc_v_energy()
              use simplex;
              implicit none;
              integer :: it, ib;

              do it = 1 , 2
                  allocate( E_fam_v(it)%blocks(1:N_blocks) );
                  do ib = 1 , N_blocks
                      allocate( E_fam_v(it)%blocks(ib)%vec(1:id_spx(ib)) );
                  end do
              end do

              return;
          end subroutine alloc_v_energy

          subroutine dealloc_v_energy()
              implicit none;
              integer :: it, ib;

              do it = 1 , 2
                  do ib = 1 , size(E_fam_v(it)%blocks)
                      deallocate( E_fam_v(it)%blocks(ib)%vec );
                  end do
                  deallocate( E_fam_v(it)%blocks );
              end do

              return;
          end subroutine dealloc_v_energy

      end module v_energy

      module fam_energies
          use dataTypes;
          implicit none;

          type(realBlockVector) , dimension(2) :: E_fam;

      contains

          subroutine alloc_fam_energies()
              use simplex;
              implicit none;
              integer :: it, ib;

              do it = 1 , 2
                  allocate( E_fam(it)%blocks(1:N_blocks) );
                  do ib = 1 , N_blocks
                      allocate( E_fam(it)%blocks(ib)%vec(1:id_spx(ib)) );
                  end do
              end do

              return;
          end subroutine alloc_fam_energies

          subroutine dealloc_fam_energies()
              implicit none;
              integer :: it, ib;

              do it = 1 , 2
                  do ib = 1 , size(E_fam(it)%blocks)
                      deallocate( E_fam(it)%blocks(ib)%vec );
                  end do
                  deallocate( E_fam(it)%blocks );
              end do

              return;
          end subroutine dealloc_fam_energies

      end module fam_energies

      module ddpc1ddme2
          implicit none;

          ! DD-PC1 parameters.
          double precision :: a_s      = - 10.0462d0; ![fm^2]
          double precision :: b_s      = -  9.1504d0; ![fm^2]
          double precision :: c_s      = -  6.4273d0; ![fm^2]
          double precision :: d_s      = +  1.3724d0; ![    ]
          double precision :: a_v      = +  5.9195d0; ![fm^2]
          double precision :: b_v      = +  8.8637d0; ![fm^2]
          double precision :: c_v      = +  0.0000d0; ![fm^2]
          double precision :: d_v      = +  0.6584d0; ![    ]
          double precision :: a_tv     = +  0.0000d0; ![fm^2]
          double precision :: b_tv     = +  1.8360d0; ![fm^2]
          double precision :: c_tv     = +  0.0000d0; ![fm^2]
          double precision :: d_tv     = +  0.6403d0; ![    ]
          double precision :: delta_s  = -  0.8149d0; ![fm^4]

          ! DD-ME2 parameters.
          double precision :: a_sigma  = +   1.388148180042941d0; ![       ]
          double precision :: b_sigma  = +   1.094300000000000d0; ![       ]
          double precision :: c_sigma  = +   1.705700000000000d0; ![       ]
          double precision :: d_sigma  = +   0.442066950716287d0; ![       ]
          double precision :: a_omega  = +   1.389291065981963d0; ![       ]
          double precision :: b_omega  = +   0.923974841420762d0; ![       ]
          double precision :: c_omega  = +   1.462000000000000d0; ![       ]
          double precision :: d_omega  = +   0.477491545490170d0; ![       ]
          double precision :: a_rho    = +   0.564700000000000d0; ![       ]
          double precision :: g0_sigma = +  10.539600000000000d0; ![       ]
          double precision :: g0_omega = +  13.018900000000000d0; ![       ]
          double precision :: g0_rho   = +   3.683600000000000d0; ![       ]
          double precision :: m_sigma  = + 550.123800000000000d0; ![MeV/c^2]
          double precision :: m_omega  = + 783.000000000000000d0; ![MeV/c^2]
          double precision :: m_rho    = + 763.000000000000000d0; ![MeV/c^2]

          double precision :: rho_sat = +  0.1520d0; ![fm^-3]

      contains

          pure double precision function alpha_s( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  alpha_s = a_s + ( b_s + c_s*x )*exp(-d_s*x);

              return;
          end function alpha_s

          pure double precision function dalpha_s( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  dalpha_s = ( c_s - d_s*( b_s + c_s*x ) )*exp(-d_s*x)*(1/rho_sat);

              return;
          end function dalpha_s

          pure double precision function ddalpha_s( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  ddalpha_s = ( -2*c_s*d_s + d_s**2*( b_s + c_s*x ) )*exp(-d_s*x)*(1/rho_sat**2);

              return;
          end function ddalpha_s

          pure double precision function alpha_v( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  alpha_v = a_v + ( b_v + c_v*x )*exp(-d_v*x);

              return;
          end function alpha_v

          pure double precision function dalpha_v( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  dalpha_v = ( c_v - d_v*( b_v + c_v*x ) )*exp(-d_v*x)*(1/rho_sat);

              return;
          end function dalpha_v

          pure double precision function ddalpha_v( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  ddalpha_v = ( -2*c_v*d_v + d_v**2*( b_v + c_v*x ) )*exp(-d_v*x)*(1/rho_sat**2);

              return;
          end function ddalpha_v

          pure double precision function alpha_tv( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  alpha_tv = a_tv + ( b_tv + c_tv*x )*exp(-d_tv*x);

              return;
          end function alpha_tv

          pure double precision function dalpha_tv( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  dalpha_tv = ( c_tv - d_tv*( b_tv + c_tv*x ) )*exp(-d_tv*x)*(1/rho_sat);

              return;
          end function dalpha_tv

          pure double precision function ddalpha_tv( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  ddalpha_tv = ( -2*c_tv*d_tv + d_tv**2*( b_tv + c_tv*x ) )*exp(-d_tv*x)*(1/rho_sat**2);

              return;
          end function ddalpha_tv

          pure double precision function g_sigma( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  g_sigma = g0_sigma * a_sigma * ( 1 + b_sigma*(x+d_sigma)**2 ) / ( 1 + c_sigma*(x+d_sigma)**2 );

              return;
          end function g_sigma

          pure double precision function dg_sigma( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  dg_sigma = g0_sigma * 2*a_sigma*(b_sigma-c_sigma) * (x+d_sigma) / ( 1 + c_sigma*(x+d_sigma)**2 )**2 * (1/rho_sat);

              return;
          end function dg_sigma

          pure double precision function ddg_sigma( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  ddg_sigma = g0_sigma * 2*a_sigma*(b_sigma-c_sigma) *  ( 1 - 3*c_sigma*(x+d_sigma)**2 ) / ( 1 + c_sigma*(x+d_sigma)**2 )**3 * (1/rho_sat**2);

              return;
          end function ddg_sigma

          pure double precision function g_omega( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  g_omega = g0_omega * a_omega * ( 1 + b_omega*(x+d_omega)**2 ) / ( 1 + c_omega*(x+d_omega)**2 );

              return;
          end function g_omega

          pure double precision function dg_omega( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  dg_omega = g0_omega * 2*a_omega*(b_omega-c_omega) * (x+d_omega) / ( 1 + c_omega*(x+d_omega)**2 )**2 * (1/rho_sat);

              return;
          end function dg_omega

          pure double precision function ddg_omega( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  ddg_omega = g0_omega * 2*a_omega*(b_omega-c_omega) *  ( 1 - 3*c_omega*(x+d_omega)**2 ) / ( 1 + c_omega*(x+d_omega)**2 )**3 * (1/rho_sat**2);

              return;
          end function ddg_omega

          pure double precision function g_rho( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  g_rho = g0_rho * exp(-a_rho*(x-1));

              return;
          end function g_rho

          pure double precision function dg_rho( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  dg_rho = g0_rho * exp(-a_rho*(x-1)) * (-a_rho/rho_sat);

              return;
          end function dg_rho

          pure double precision function ddg_rho( rho_v )
              implicit none;
              double precision , intent(in) :: rho_v;
              double precision              :: x;

                  x = rho_v / rho_sat;

                  ddg_rho = g0_rho * exp(-a_rho*(x-1)) * (-a_rho/rho_sat)**2;

              return;
          end function ddg_rho

      end module ddpc1ddme2

      module pairparams
         implicit none;

          double precision :: G_pairing = + 728.000000000d0; ![MeV/fm^3]
          double precision :: a_pairing = +   0.644204936d0; ![   fm   ]

      end module pairparams

      module quadrature
          implicit none;

          double precision , dimension(:)   , allocatable :: zb_fam;
          double precision , dimension(:)   , allocatable :: wz;
          double precision , dimension(:)   , allocatable :: rb_fam;
          double precision , dimension(:)   , allocatable :: wr;
          double precision , dimension(:,:) , allocatable :: wzwr;

      contains

          subroutine alloc_quadrature()
              use fam_input;
              implicit none;

              allocate( zb_fam(1:NGH)       );
              allocate(     wz(1:NGH)       );
              allocate( rb_fam(1:NGL)       );
              allocate(     wr(1:NGL)       );
              allocate(   wzwr(1:NGH,1:NGL) );

              return;
          end subroutine alloc_quadrature

          subroutine dealloc_quadrature()
              implicit none;

              deallocate( zb_fam );
              deallocate( wz     );
              deallocate( rb_fam );
              deallocate( wr     );
              deallocate( wzwr   );

              return;
          end subroutine dealloc_quadrature

      end module quadrature

      module basis
          use dataTypes;
          implicit none;

          double precision                               ::     bz;
          double precision                               ::     bp;
          type(real2Darray) , dimension(:) , allocatable ::  phi_z;
          type(real2Darray) , dimension(:) , allocatable :: dphi_z;
          type(real2Darray) , dimension(:) , allocatable ::  phi_r;
          type(real2Darray) , dimension(:) , allocatable :: dphi_r;

      contains

          subroutine alloc_basis()
              use fam_input;
              use simplex;
              implicit none;
              integer :: ib;

              allocate( phi_z(1:N_blocks) );
              do ib = 1 , N_blocks
                  allocate( phi_z(ib)%mat( -NGH:+NGH , 1:id_spx(ib) ) );
              end do

              allocate( dphi_z(1:N_blocks) );
              do ib = 1 , N_blocks
                  allocate( dphi_z(ib)%mat( -NGH:+NGH , 1:id_spx(ib) ) );
              end do

              allocate( phi_r(1:N_blocks) );
              do ib = 1 , N_blocks
                  allocate( phi_r(ib)%mat( 1:NGL , 1:id_spx(ib) ) );
              end do

              allocate( dphi_r(1:N_blocks) );
              do ib = 1 , N_blocks
                  allocate( dphi_r(ib)%mat( 1:NGL , 1:id_spx(ib) ) );
              end do

              return;
          end subroutine alloc_basis

          subroutine dealloc_basis()
              implicit none;
              integer :: ib;

              do ib = 1 , size(phi_z)
                  deallocate( phi_z(ib)%mat );
              end do
              deallocate( phi_z );

              do ib = 1 , size(dphi_z)
                  deallocate( dphi_z(ib)%mat );
              end do
              deallocate( dphi_z );

              do ib = 1 , size(phi_r)
                  deallocate( phi_r(ib)%mat );
              end do
              deallocate( phi_r );

              do ib = 1 , size(dphi_r)
                  deallocate( dphi_r(ib)%mat );
              end do
              deallocate( dphi_r );

              return;
          end subroutine dealloc_basis

      end module basis

      module wbasis
          use dataTypes;
          implicit none;

          type(real3Darray) , dimension(:) , allocatable :: wPhi;

      contains

          subroutine alloc_wbasis()
              use fam_input;
              use simplex;
              implicit none;
              integer :: ib;

              allocate( wPhi(1:N_blocks) );
              do ib = 1 , N_blocks
                  allocate( wPhi(ib)%arr( 1:NGH , 1:NGL , 1:id_spx(ib) ) );
              end do

              return;
          end subroutine alloc_wbasis

          subroutine dealloc_wbasis()
              implicit none;
              integer :: ib;

              do ib = 1 , size(wPhi)
                  deallocate( wPhi(ib)%arr );
              end do
              deallocate( wPhi );

              return;
          end subroutine dealloc_wbasis

      end module wbasis

      module PHI
          use dataTypes;
          implicit none;

          double precision      , dimension(:,:) , allocatable :: U_PHI;
          double precision      , dimension(:,:) , allocatable :: V_PHI;

          integer               ,                  parameter   :: nbsize = 64;
          double precision      , dimension(:,:) , allocatable :: AU;
          double precision      , dimension(:,:) , allocatable :: UtAU;
          double precision      , dimension(:,:) , allocatable :: UtAUV;
          double precision      , dimension(:)   , allocatable :: VtUtAUV;

          type(realBlockMatrix)                                :: Ar_tmp;
          type(realBlockMatrix)                                :: Ai_tmp;


      contains

          subroutine alloc_PHI()
              use fam_input;
              use simplex;
              implicit none;
              integer   :: ktrunc, ncoord;
              integer   :: ib1, ib2, i1, i2;
              character :: fg1, fg2;
              integer   :: ml1, ml2;

              ! ktrunc is actually maximal block dimension of block matrices U,V.
              ! In fact, this is the dimension of the omega^pi = 1/2^+ block.
              ktrunc = ((n0f+2)*(n0f+3))/2;
              ncoord = (2*NGH)*NGL;
              allocate(   U_PHI( 1:N_total , 1:ktrunc ) );
              allocate(   V_PHI( 1:ktrunc  , 1:ncoord ) );
              allocate(      AU( 1:N_total , 1:ktrunc ) );
              allocate(    UtAU( 1:ktrunc  , 1:ktrunc ) );
              allocate(   UtAUV( 1:ktrunc  , 1:nbsize ) );
              allocate( VtUtAUV( 1:ncoord             ) );


              allocate( Ar_tmp%nnzblocks(1:N_blocks,1:N_blocks) );
              allocate( Ai_tmp%nnzblocks(1:N_blocks,1:N_blocks) );

              Ar_tmp%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
              Ai_tmp%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
              do ib2 = 1 , N_blocks
                  do ib1 = 1 , ib2
                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , merge( i2 , id_spx(ib1) , ib1==ib2 )

                              fg1 = fg_spx(ib1)%index(i1);
                              ml1 = ml_spx(ib1)%index(i1);

                              fg2 = fg_spx(ib2)%index(i2);
                              ml2 = ml_spx(ib2)%index(i2);

                              if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then
                                  Ar_tmp%nnzblocks(ib1,ib2) = .true.;
                                  Ai_tmp%nnzblocks(ib1,ib2) = .true.;
                              end if
                              if( fg1/=fg2 .and. ( abs(ml1-ml2)==K_multipole .or. abs(ml1+ml2+1)==K_multipole ) ) then
                                  Ar_tmp%nnzblocks(ib1,ib2) = .true.;
                                  Ai_tmp%nnzblocks(ib1,ib2) = .true.;
                              end if

                          end do
                      end do
                  end do
              end do


              allocate( Ar_tmp%blocks(1:N_blocks,1:N_blocks) );
              allocate( Ai_tmp%blocks(1:N_blocks,1:N_blocks) );
              do ib2 = 1 , N_blocks
                  do ib1 = 1 , ib2

                      if( Ar_tmp%nnzblocks(ib1,ib2) ) &
                          allocate( Ar_tmp%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                      if( Ai_tmp%nnzblocks(ib1,ib2) ) &
                          allocate( Ai_tmp%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                  end do
              end do


              return;
          end subroutine alloc_PHI

          subroutine dealloc_PHI()
              implicit none;
              integer :: ib1, ib2;

              deallocate( U_PHI );
              deallocate( V_PHI );
              deallocate( AU    );
              deallocate( UtAU  );
              deallocate( UtAUV );


              do ib2 = 1 , size(Ar_tmp%nnzblocks,2)
                  do ib1 = 1 , size(Ar_tmp%nnzblocks,1)
                      if( Ar_tmp%nnzblocks(ib1,ib2) ) deallocate( Ar_tmp%blocks(ib1,ib2)%mat );
                  end do
              end do
              deallocate( Ar_tmp%blocks    );
              deallocate( Ar_tmp%nnzblocks );

              do ib2 = 1 , size(Ai_tmp%nnzblocks,2)
                  do ib1 = 1 , size(Ai_tmp%nnzblocks,1)
                      if( Ai_tmp%nnzblocks(ib1,ib2) ) deallocate( Ai_tmp%blocks(ib1,ib2)%mat );
                  end do
              end do
              deallocate( Ai_tmp%blocks    );
              deallocate( Ai_tmp%nnzblocks );


              return;
          end subroutine dealloc_PHI

      end module PHI

      module gs_dens
          implicit none;

          double precision , dimension(:,:,:) , allocatable :: rhov_GS;
          double precision , dimension(:,:)   , allocatable :: rhos_GS;

      contains

          subroutine alloc_gs_dens()
              use fam_input;
              implicit none;

              allocate( rhov_GS( -NGH:+NGH , 1:NGL , 1:2 ) );
              allocate( rhos_GS( -NGH:+NGH , 1:NGL       ) );

              return;
          end subroutine alloc_gs_dens

          subroutine dealloc_gs_dens()
              implicit none;

              deallocate( rhov_GS );
              deallocate( rhos_GS );

              return;
          end subroutine dealloc_gs_dens

      end module gs_dens

      module gs_mesons
          implicit none;

          double precision , dimension(:,:) , allocatable :: sigma_GS;
          double precision , dimension(:,:) , allocatable :: omega_GS;
          double precision , dimension(:,:) , allocatable ::   rho_GS;

      contains

          subroutine alloc_gs_mesons()
              use fam_input;
              implicit none;

              allocate( sigma_GS( -NGH:+NGH , 1:NGL ) );
              allocate( omega_GS( -NGH:+NGH , 1:NGL ) );
              allocate(   rho_GS( -NGH:+NGH , 1:NGL ) );

              return;
          end subroutine alloc_gs_mesons

          subroutine dealloc_gs_mesons()
              implicit none;

              deallocate( sigma_GS );
              deallocate( omega_GS );
              deallocate(   rho_GS );

              return;
          end subroutine dealloc_gs_mesons

      end module gs_mesons

      module tempBlockMatrix
          use dataTypes;
          implicit none;

          ! This module is used in all multiplications with u,v matrices
          ! as an intermediate step. For example, to calculate u*x*v', we first
          ! have to calculate and store tmpMat=(u*x) and then calculate (u*x)*v'.
          type(complexBlockMatrix) :: tmpMat;

      contains

          subroutine alloc_tempBlockMatrix()
              use fam_input;
              use simplex;
              implicit none;
              integer   :: ib1, ib2, i1, i2;
              character :: fg1, fg2;
              integer   :: ml1, ml2;


                  allocate( tmpMat%nnzblocks(1:N_blocks,1:N_blocks) );

                  tmpMat%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)

                                  fg1 = fg_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then
                                      tmpMat%nnzblocks(ib1,ib2) = .true.;
                                  end if
                                  if( fg1/=fg2 .and. ( abs(ml1-ml2)==K_multipole .or. abs(ml1+ml2+1)==K_multipole ) ) then
                                      tmpMat%nnzblocks(ib1,ib2) = .true.;
                                  end if

                              end do
                          end do
                      end do
                  end do


                  allocate( tmpMat%blocks(1:N_blocks,1:N_blocks) );
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks
                          if( tmpMat%nnzblocks(ib1,ib2) ) &
                              allocate( tmpMat%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );
                      end do
                  end do


              return;
          end subroutine alloc_tempBlockMatrix

          subroutine dealloc_tempBlockMatrix()
              implicit none;
              integer :: ib1, ib2;


                  do ib2 = 1 , size(tmpMat%nnzblocks,2)
                      do ib1 = 1 , size(tmpMat%nnzblocks,1)
                          if( tmpMat%nnzblocks(ib1,ib2) ) deallocate( tmpMat%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( tmpMat%blocks    );
                  deallocate( tmpMat%nnzblocks );


              return;
          end subroutine dealloc_tempBlockMatrix

      end module tempBlockMatrix

      module fmatrix
          use dataTypes;
          implicit none;

          type(complexBlockMatrix) , dimension(2) :: f1_JK;
          type(complexBlockMatrix) , dimension(2) :: f2_JK;

      contains

          subroutine alloc_fmatrix()
              use fam_input;
              use simplex;
              implicit none;
              integer   :: it, ib1, ib2, i1, i2;
              character :: fg1, fg2;
              integer   :: ml1, ml2;


              do it = 1 , 2

                  allocate( f1_JK(it)%nnzblocks(1:N_blocks,1:N_blocks) );
                  allocate( f2_JK(it)%nnzblocks(1:N_blocks,1:N_blocks) );

                  f1_JK(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  f2_JK(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)

                                  fg1 = fg_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then
                                      f1_JK(it)%nnzblocks(ib1,ib2) = .true.;
                                      f2_JK(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if

                              end do
                          end do
                      end do
                  end do


                  allocate( f1_JK(it)%blocks(1:N_blocks,1:N_blocks) );
                  allocate( f2_JK(it)%blocks(1:N_blocks,1:N_blocks) );
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks

                          if( f1_JK(it)%nnzblocks(ib1,ib2) ) &
                              allocate( f1_JK(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                          if( f2_JK(it)%nnzblocks(ib1,ib2) ) &
                              allocate( f2_JK(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                      end do
                  end do

              end do

              return;
          end subroutine alloc_fmatrix

          subroutine dealloc_fmatrix()
              implicit none;
              integer :: it, ib1, ib2;

              do it = 1 , 2

                  do ib2 = 1 , size(f1_JK(it)%nnzblocks,2)
                      do ib1 = 1 , size(f1_JK(it)%nnzblocks,1)
                          if( f1_JK(it)%nnzblocks(ib1,ib2) ) deallocate( f1_JK(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( f1_JK(it)%blocks    );
                  deallocate( f1_JK(it)%nnzblocks );

                  do ib2 = 1 , size(f2_JK(it)%nnzblocks,2)
                      do ib1 = 1 , size(f2_JK(it)%nnzblocks,1)
                          if( f2_JK(it)%nnzblocks(ib1,ib2) ) deallocate( f2_JK(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( f2_JK(it)%blocks    );
                  deallocate( f2_JK(it)%nnzblocks );

              end do

              return;
          end subroutine dealloc_fmatrix

      end module fmatrix

      module f02f20matrix
          use dataTypes;
          implicit none;

          type(complexBlockMatrix) , dimension(2) :: f20;
          type(complexBlockMatrix) , dimension(2) :: f02;

      contains

          subroutine alloc_f02f20matrix()
              use fam_input;
              use simplex;
              implicit none;
              integer   :: it, ib1, ib2, i1, i2;
              character :: fg1, fg2;
              integer   :: ml1, ml2;


              do it = 1 , 2

                  allocate( f20(it)%nnzblocks(1:N_blocks,1:N_blocks) );
                  allocate( f02(it)%nnzblocks(1:N_blocks,1:N_blocks) );

                  f20(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  f02(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)

                                  fg1 = fg_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then
                                      f20(it)%nnzblocks(ib1,ib2) = .true.;
                                      f02(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if

                              end do
                          end do
                      end do
                  end do


                  allocate( f20(it)%blocks(1:N_blocks,1:N_blocks) );
                  allocate( f02(it)%blocks(1:N_blocks,1:N_blocks) );
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks

                          if( f20(it)%nnzblocks(ib1,ib2) ) &
                              allocate( f20(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                          if( f02(it)%nnzblocks(ib1,ib2) ) &
                              allocate( f02(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                      end do
                  end do

              end do

              return;
          end subroutine alloc_f02f20matrix

          subroutine dealloc_f02f20matrix()
              implicit none;
              integer :: it, ib1, ib2;

              do it = 1 , 2

                  do ib2 = 1 , size(f20(it)%nnzblocks,2)
                      do ib1 = 1 , size(f20(it)%nnzblocks,1)
                          if( f20(it)%nnzblocks(ib1,ib2) ) deallocate( f20(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( f20(it)%blocks    );
                  deallocate( f20(it)%nnzblocks );

                  do ib2 = 1 , size(f02(it)%nnzblocks,2)
                      do ib1 = 1 , size(f02(it)%nnzblocks,1)
                          if( f02(it)%nnzblocks(ib1,ib2) ) deallocate( f02(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( f02(it)%blocks    );
                  deallocate( f02(it)%nnzblocks );

              end do

              return;
          end subroutine dealloc_f02f20matrix

      end module f02f20matrix

      module spurious
          use dataTypes;
          implicit none;

          double complex           , dimension(2) :: RcmPcm_commutator;
          type(complexBlockMatrix) , dimension(2) :: r20;
          type(complexBlockMatrix) , dimension(2) :: p20;

      contains

          subroutine alloc_spurious()
              use fam_input;
              use simplex;
              implicit none;
              integer   :: it, ib1, ib2, i1, i2;
              character :: fg1, fg2;
              integer   :: nz1, nz2;
              integer   :: nr1, nr2;
              integer   :: ml1, ml2;


              do it = 1 , 2

                  allocate( r20(it)%nnzblocks(1:N_blocks,1:N_blocks) );
                  allocate( p20(it)%nnzblocks(1:N_blocks,1:N_blocks) );

                  r20(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  p20(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)

                                  fg1 = fg_spx(ib1)%index(i1);
                                  nz1 = nz_spx(ib1)%index(i1);
                                  nr1 = nr_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  nz2 = nz_spx(ib2)%index(i2);
                                  nr2 = nr_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( K_multipole == 0 ) then
                                      if( fg1==fg2 .and. nr1==nr2 .and. abs(ml1-ml2)==K_multipole ) then
                                          r20(it)%nnzblocks(ib1,ib2) = .true.;
                                          p20(it)%nnzblocks(ib1,ib2) = .true.;
                                      end if
                                  end if
                                  if( K_multipole == 1 ) then
                                      if( fg1==fg2 .and. nz1==nz2 .and. abs(ml1-ml2)==K_multipole ) then
                                          r20(it)%nnzblocks(ib1,ib2) = .true.;
                                          p20(it)%nnzblocks(ib1,ib2) = .true.;
                                      end if
                                  end if

                              end do
                          end do
                      end do
                  end do


                  allocate( r20(it)%blocks(1:N_blocks,1:N_blocks) );
                  allocate( p20(it)%blocks(1:N_blocks,1:N_blocks) );
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks

                          if( r20(it)%nnzblocks(ib1,ib2) ) &
                              allocate( r20(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                          if( p20(it)%nnzblocks(ib1,ib2) ) &
                              allocate( p20(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                      end do
                  end do

              end do

              return;
          end subroutine alloc_spurious

          subroutine dealloc_spurious()
              implicit none;
              integer :: it, ib1, ib2;

              do it = 1 , 2

                  do ib2 = 1 , size(r20(it)%nnzblocks,2)
                      do ib1 = 1 , size(r20(it)%nnzblocks,1)
                          if( r20(it)%nnzblocks(ib1,ib2) ) deallocate( r20(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( r20(it)%blocks    );
                  deallocate( r20(it)%nnzblocks );

                  do ib2 = 1 , size(p20(it)%nnzblocks,2)
                      do ib1 = 1 , size(p20(it)%nnzblocks,1)
                          if( p20(it)%nnzblocks(ib1,ib2) ) deallocate( p20(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( p20(it)%blocks    );
                  deallocate( p20(it)%nnzblocks );

              end do

              return;
          end subroutine dealloc_spurious

      end module spurious

      module mesmat
          use dataTypes;
          implicit none;

          integer           , dimension(:) , allocatable :: nP;
          type(real2Darray) , dimension(:) , allocatable :: Psig;
          type(real2Darray) , dimension(:) , allocatable :: Pome;
          type(real2Darray) , dimension(:) , allocatable :: Prho;

      contains

          subroutine alloc_mesmat()
              use fam_input;
              implicit none;
              integer :: K, N;
              integer :: ncoord;

                  ! Number of (nz>=0,nr>=0) pairs such that nz + 2*nr + K <= n0b.
                  allocate( nP(0:J_MAX+1) );
                  do K = 0 , J_MAX+1
                      call assert( n0b-K>0 , 'n0b is too small!' );

                      N     = n0b-K;
                      nP(K) = ( (N+1)*(N+3) + 1 - mod(N,2) )/4;
                  end do

                  allocate( Psig(0:J_MAX+1) );
                  allocate( Pome(0:J_MAX+1) );
                  allocate( Prho(0:J_MAX+1) );
                  do K = 0 , J_MAX+1
                      ncoord = (2*NGH)*NGL;

                      allocate( Psig(K)%mat( 1:nP(K) , 1:ncoord ) );
                      allocate( Pome(K)%mat( 1:nP(K) , 1:ncoord ) );
                      allocate( Prho(K)%mat( 1:nP(K) , 1:ncoord ) );
                  end do

              return;
          end subroutine alloc_mesmat

          subroutine dealloc_mesmat()
              use fam_input;
              implicit none;
              integer :: K;

                  do K = 0 , J_MAX+1
                      deallocate( Psig(K)%mat );
                      deallocate( Pome(K)%mat );
                      deallocate( Prho(K)%mat );
                  end do
                  deallocate( nP   );
                  deallocate( Psig );
                  deallocate( Pome );
                  deallocate( Prho );

              return;
          end subroutine dealloc_mesmat

      end module mesmat

      module fam_green
          implicit none;

          double precision , dimension(:,:,:,:) , allocatable :: G1;
          double precision , dimension(:,:,:,:) , allocatable :: G2;
          double precision , dimension(:,:,:,:) , allocatable :: G3;

      contains

          subroutine alloc_fam_green()
              use fam_input;
              implicit none;

                  allocate( G1( -NGH:+NGH , 1:NGL , -NGH:+NGH , 1:NGL ) );
                  allocate( G2( -NGH:+NGH , 1:NGL , -NGH:+NGH , 1:NGL ) );
                  allocate( G3( -NGH:+NGH , 1:NGL , -NGH:+NGH , 1:NGL ) );

              return;
          end subroutine alloc_fam_green

          subroutine dealloc_fam_green()
              implicit none;

                  deallocate( G1 );
                  deallocate( G2 );
                  deallocate( G3 );

              return;
          end subroutine dealloc_fam_green

      end module fam_green

      module Wpairing
          implicit none;

          double precision , dimension(:) , allocatable :: W;

      contains

          subroutine alloc_Wpairing()
              use fam_input;
              use simplex;
              implicit none;
              integer   :: ib1, ib2, i1, i2;
              character :: fg1, fg2;
              integer   :: nz1, nz2;
              integer   :: nr1, nr2;
              integer   :: ml1, ml2;
              integer   :: N_r, N_z;
              integer   :: Wsize;

              Wsize = 0;
              do ib2 = 1 , N_blocks
                  do ib1 = 1 , N_blocks
                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , id_spx(ib1)
                              if( (ib1<ib2) .or. ( (ib1==ib2) .and. (i1<=i2) ) ) then ! Loop over the upper triangle.

                                  fg1 = fg_spx(ib1)%index(i1);
                                  nz1 = nz_spx(ib1)%index(i1);
                                  nr1 = nr_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  nz2 = nz_spx(ib2)%index(i2);
                                  nr2 = nr_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( all([fg1,fg2]==['f','f']) .and. abs(ml1-ml2)==K_multipole ) then
                                      do N_r = 0 , ( 2*nr1+abs(ml1) + 2*nr2+abs(ml2) - K_multipole )/2
                                          do N_z = 0 , nz1+nz2
                                              if( mod(N_z,2) == mod(nz1+nz2,2) ) then
                                                  Wsize = Wsize + 1;
                                              end if
                                          end do
                                      end do
                                  end if

                              end if
                          end do
                      end do
                  end do
              end do

              allocate( W(1:Wsize) );

              return;
          end subroutine alloc_Wpairing

          subroutine dealloc_Wpairing()
              implicit none;

                  deallocate( W );

              return;
          end subroutine dealloc_Wpairing

      end module Wpairing

      module dh
          use dataTypes;
          implicit none;

          type(complexBlockMatrix) , dimension(2) :: dh1;
          type(complexBlockMatrix) , dimension(2) :: dh2;

      contains

          subroutine alloc_dh()
              use fam_input;
              use simplex;
              implicit none;
              integer   :: it, ib1, ib2, i1, i2;
              character :: fg1, fg2;
              integer   :: ml1, ml2;


              do it = 1 , 2

                  allocate( dh1(it)%nnzblocks(1:N_blocks,1:N_blocks) );
                  allocate( dh2(it)%nnzblocks(1:N_blocks,1:N_blocks) );

                  dh1(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  dh2(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)

                                  fg1 = fg_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then
                                      dh1(it)%nnzblocks(ib1,ib2) = .true.;
                                      dh2(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if
                                  if( fg1/=fg2 .and. ( abs(ml1-ml2)==K_multipole .or. abs(ml1+ml2+1)==K_multipole ) ) then
                                      dh1(it)%nnzblocks(ib1,ib2) = .true.;
                                      dh2(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if

                              end do
                          end do
                      end do
                  end do


                  allocate( dh1(it)%blocks(1:N_blocks,1:N_blocks) );
                  allocate( dh2(it)%blocks(1:N_blocks,1:N_blocks) );
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks

                          if( dh1(it)%nnzblocks(ib1,ib2) ) &
                              allocate( dh1(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                          if( dh2(it)%nnzblocks(ib1,ib2) ) &
                              allocate( dh2(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                      end do
                  end do

              end do

              return;
          end subroutine alloc_dh

          subroutine dealloc_dh()
              implicit none;
              integer :: it, ib1, ib2;

              do it = 1 , 2

                  do ib2 = 1 , size(dh1(it)%nnzblocks,2)
                      do ib1 = 1 , size(dh1(it)%nnzblocks,1)
                          if( dh1(it)%nnzblocks(ib1,ib2) ) deallocate( dh1(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( dh1(it)%blocks    );
                  deallocate( dh1(it)%nnzblocks );

                  do ib2 = 1 , size(dh2(it)%nnzblocks,2)
                      do ib1 = 1 , size(dh2(it)%nnzblocks,1)
                          if( dh2(it)%nnzblocks(ib1,ib2) ) deallocate( dh2(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( dh2(it)%blocks    );
                  deallocate( dh2(it)%nnzblocks );

              end do

              return;
          end subroutine dealloc_dh

      end module dh

      module dDelta
          use dataTypes;
          implicit none;

          type(complexBlockMatrix) , dimension(2) :: dDelta1_pl;
          type(complexBlockMatrix) , dimension(2) :: dDelta1_mi;

      contains

          subroutine alloc_dDelta()
              use fam_input;
              use simplex;
              implicit none;
              integer   :: it, ib1, ib2, i1, i2;
              character :: fg1, fg2;
              integer   :: ml1, ml2;


              do it = 1 , 2

                  allocate( dDelta1_pl(it)%nnzblocks(1:N_blocks,1:N_blocks) );
                  allocate( dDelta1_mi(it)%nnzblocks(1:N_blocks,1:N_blocks) );

                  dDelta1_pl(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  dDelta1_mi(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)

                                  fg1 = fg_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( all([fg1,fg2]==['f','f']) .and. abs(ml1-ml2)==K_multipole ) then
                                      dDelta1_pl(it)%nnzblocks(ib1,ib2) = .true.;
                                      dDelta1_mi(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if

                              end do
                          end do
                      end do
                  end do


                  allocate( dDelta1_pl(it)%blocks(1:N_blocks,1:N_blocks) );
                  allocate( dDelta1_mi(it)%blocks(1:N_blocks,1:N_blocks) );
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks

                          if( dDelta1_pl(it)%nnzblocks(ib1,ib2) ) &
                              allocate( dDelta1_pl(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                          if( dDelta1_mi(it)%nnzblocks(ib1,ib2) ) &
                              allocate( dDelta1_mi(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                      end do
                  end do

              end do

              return;
          end subroutine alloc_dDelta

          subroutine dealloc_dDelta()
              implicit none;
              integer :: it, ib1, ib2;

              do it = 1 , 2

                  do ib2 = 1 , size(dDelta1_pl(it)%nnzblocks,2)
                      do ib1 = 1 , size(dDelta1_pl(it)%nnzblocks,1)
                          if( dDelta1_pl(it)%nnzblocks(ib1,ib2) ) deallocate( dDelta1_pl(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( dDelta1_pl(it)%blocks    );
                  deallocate( dDelta1_pl(it)%nnzblocks );

                  do ib2 = 1 , size(dDelta1_mi(it)%nnzblocks,2)
                      do ib1 = 1 , size(dDelta1_mi(it)%nnzblocks,1)
                          if( dDelta1_mi(it)%nnzblocks(ib1,ib2) ) deallocate( dDelta1_mi(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( dDelta1_mi(it)%blocks    );
                  deallocate( dDelta1_mi(it)%nnzblocks );

              end do

              return;
          end subroutine dealloc_dDelta

      end module dDelta

      module ddens
          implicit none;

          double complex , dimension(:,:,:) , allocatable :: drho_v;
          double complex , dimension(:,:)   , allocatable :: drho_s;

      contains

          subroutine alloc_ddens()
              use fam_input;
              implicit none;

              allocate( drho_v( -NGH:+NGH , 1:NGL , 1:2 ) );
              allocate( drho_s( -NGH:+NGH , 1:NGL       ) );

              return;
          end subroutine alloc_ddens

          subroutine dealloc_ddens()
              implicit none;

              deallocate( drho_v );
              deallocate( drho_s );

              return;
          end subroutine dealloc_ddens

      end module ddens

      module dh02dh20matrix
          use dataTypes;
          implicit none;

          type(complexBlockMatrix) , dimension(2) :: dh20;
          type(complexBlockMatrix) , dimension(2) :: dh02;

      contains

          subroutine alloc_dh02dh20matrix()
              use fam_input;
              use simplex;
              implicit none;
              integer   :: it, ib1, ib2, i1, i2;
              character :: fg1, fg2;
              integer   :: ml1, ml2;


              do it = 1 , 2

                  allocate( dh20(it)%nnzblocks(1:N_blocks,1:N_blocks) );
                  allocate( dh02(it)%nnzblocks(1:N_blocks,1:N_blocks) );

                  dh20(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  dh02(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)

                                  fg1 = fg_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then
                                      dh20(it)%nnzblocks(ib1,ib2) = .true.;
                                      dh02(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if
                                  if( fg1/=fg2 .and. ( abs(ml1-ml2)==K_multipole .or. abs(ml1+ml2+1)==K_multipole ) ) then
                                      dh20(it)%nnzblocks(ib1,ib2) = .true.;
                                      dh02(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if

                              end do
                          end do
                      end do
                  end do


                  allocate( dh20(it)%blocks(1:N_blocks,1:N_blocks) );
                  allocate( dh02(it)%blocks(1:N_blocks,1:N_blocks) );
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks

                          if( dh20(it)%nnzblocks(ib1,ib2) ) &
                              allocate( dh20(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                          if( dh02(it)%nnzblocks(ib1,ib2) ) &
                              allocate( dh02(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                      end do
                  end do

              end do

              return;
          end subroutine alloc_dh02dh20matrix

          subroutine dealloc_dh02dh20matrix()
              implicit none;
              integer :: it, ib1, ib2;

              do it = 1 , 2

                  do ib2 = 1 , size(dh20(it)%nnzblocks,2)
                      do ib1 = 1 , size(dh20(it)%nnzblocks,1)
                          if( dh20(it)%nnzblocks(ib1,ib2) ) deallocate( dh20(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( dh20(it)%blocks    );
                  deallocate( dh20(it)%nnzblocks );

                  do ib2 = 1 , size(dh02(it)%nnzblocks,2)
                      do ib1 = 1 , size(dh02(it)%nnzblocks,1)
                          if( dh02(it)%nnzblocks(ib1,ib2) ) deallocate( dh02(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( dh02(it)%blocks    );
                  deallocate( dh02(it)%nnzblocks );

              end do

              return;
          end subroutine dealloc_dh02dh20matrix

      end module dh02dh20matrix

      module dh11matrix
          use dataTypes;
          implicit none;

          type(complexBlockMatrix) , dimension(2) :: dh11_1;
          type(complexBlockMatrix) , dimension(2) :: dh11_2;

      contains

          subroutine alloc_dh11matrix()
              use fam_input;
              use simplex;
              implicit none;
              integer   :: it, ib1, ib2, i1, i2;
              character :: fg1, fg2;
              integer   :: ml1, ml2;


              do it = 1 , 2

                  allocate( dh11_1(it)%nnzblocks(1:N_blocks,1:N_blocks) );
                  allocate( dh11_2(it)%nnzblocks(1:N_blocks,1:N_blocks) );

                  dh11_1(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  dh11_2(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)

                                  fg1 = fg_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then
                                      dh11_1(it)%nnzblocks(ib1,ib2) = .true.;
                                      dh11_2(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if
                                  if( fg1/=fg2 .and. ( abs(ml1-ml2)==K_multipole .or. abs(ml1+ml2+1)==K_multipole ) ) then
                                      dh11_1(it)%nnzblocks(ib1,ib2) = .true.;
                                      dh11_2(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if

                              end do
                          end do
                      end do
                  end do


                  allocate( dh11_1(it)%blocks(1:N_blocks,1:N_blocks) );
                  allocate( dh11_2(it)%blocks(1:N_blocks,1:N_blocks) );
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks

                          if( dh11_1(it)%nnzblocks(ib1,ib2) ) &
                              allocate( dh11_1(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                          if( dh11_2(it)%nnzblocks(ib1,ib2) ) &
                              allocate( dh11_2(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                      end do
                  end do

              end do

              return;
          end subroutine alloc_dh11matrix

          subroutine dealloc_dh11matrix()
              implicit none;
              integer :: it, ib1, ib2;

              do it = 1 , 2

                  do ib2 = 1 , size(dh11_1(it)%nnzblocks,2)
                      do ib1 = 1 , size(dh11_1(it)%nnzblocks,1)
                          if( dh11_1(it)%nnzblocks(ib1,ib2) ) deallocate( dh11_1(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( dh11_1(it)%blocks    );
                  deallocate( dh11_1(it)%nnzblocks );

                  do ib2 = 1 , size(dh11_2(it)%nnzblocks,2)
                      do ib1 = 1 , size(dh11_2(it)%nnzblocks,1)
                          if( dh11_2(it)%nnzblocks(ib1,ib2) ) deallocate( dh11_2(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( dh11_2(it)%blocks    );
                  deallocate( dh11_2(it)%nnzblocks );

              end do

              return;
          end subroutine dealloc_dh11matrix

      end module dh11matrix

      module xyfam
          use dataTypes;
          implicit none;

          type(complexBlockMatrix) , dimension(2) :: x_fam;
          type(complexBlockMatrix) , dimension(2) :: y_fam;

      contains

          subroutine alloc_xyfam()
              use fam_input;
              use simplex;
              implicit none;
              integer   :: it, ib1, ib2, i1, i2;
              character :: fg1, fg2;
              integer   :: ml1, ml2;

              do it = 1 , 2

                  allocate( x_fam(it)%nnzblocks(1:N_blocks,1:N_blocks) );
                  allocate( y_fam(it)%nnzblocks(1:N_blocks,1:N_blocks) );

                  x_fam(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  y_fam(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)

                                  fg1 = fg_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then
                                      x_fam(it)%nnzblocks(ib1,ib2) = .true.;
                                      y_fam(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if
                                  if( fg1/=fg2 .and. ( abs(ml1-ml2)==K_multipole .or. abs(ml1+ml2+1)==K_multipole ) ) then
                                      x_fam(it)%nnzblocks(ib1,ib2) = .true.;
                                      y_fam(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if

                              end do
                          end do
                      end do
                  end do


                  allocate( x_fam(it)%blocks(1:N_blocks,1:N_blocks) );
                  allocate( y_fam(it)%blocks(1:N_blocks,1:N_blocks) );
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks

                          if( x_fam(it)%nnzblocks(ib1,ib2) ) &
                              allocate( x_fam(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                          if( y_fam(it)%nnzblocks(ib1,ib2) ) &
                              allocate( y_fam(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                      end do
                  end do

              end do

              return;
          end subroutine alloc_xyfam

          subroutine dealloc_xyfam()
              implicit none;
              integer :: it, ib1, ib2;

              do it = 1 , 2

                  do ib2 = 1 , size(x_fam(it)%nnzblocks,2)
                      do ib1 = 1 , size(x_fam(it)%nnzblocks,1)
                          if( x_fam(it)%nnzblocks(ib1,ib2) ) deallocate( x_fam(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( x_fam(it)%blocks    );
                  deallocate( x_fam(it)%nnzblocks );

                  do ib2 = 1 , size(y_fam(it)%nnzblocks,2)
                      do ib1 = 1 , size(y_fam(it)%nnzblocks,1)
                          if( y_fam(it)%nnzblocks(ib1,ib2) ) deallocate( y_fam(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( y_fam(it)%blocks    );
                  deallocate( y_fam(it)%nnzblocks );

              end do

              return;
          end subroutine dealloc_xyfam

      end module xyfam

      module drho
          use dataTypes;
          implicit none;

          type(complexBlockMatrix) , dimension(2) :: drho1;
          type(complexBlockMatrix) , dimension(2) :: drho2;

      contains

          subroutine alloc_drho()
              use fam_input;
              use simplex;
              implicit none;
              integer   :: it, ib1, ib2, i1, i2;
              character :: fg1, fg2;
              integer   :: ml1, ml2;

              do it = 1 , 2

                  allocate( drho1(it)%nnzblocks(1:N_blocks,1:N_blocks) );
                  allocate( drho2(it)%nnzblocks(1:N_blocks,1:N_blocks) );

                  drho1(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  drho2(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)

                                  fg1 = fg_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then
                                      drho1(it)%nnzblocks(ib1,ib2) = .true.;
                                      drho2(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if
                                  if( fg1/=fg2 .and. ( abs(ml1-ml2)==K_multipole .or. abs(ml1+ml2+1)==K_multipole ) ) then
                                      drho1(it)%nnzblocks(ib1,ib2) = .true.;
                                      drho2(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if

                              end do
                          end do
                      end do
                  end do


                  allocate( drho1(it)%blocks(1:N_blocks,1:N_blocks) );
                  allocate( drho2(it)%blocks(1:N_blocks,1:N_blocks) );
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks

                          if( drho1(it)%nnzblocks(ib1,ib2) ) &
                              allocate( drho1(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                          if( drho2(it)%nnzblocks(ib1,ib2) ) &
                              allocate( drho2(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                      end do
                  end do

              end do

              return;
          end subroutine alloc_drho

          subroutine dealloc_drho()
              implicit none;
              integer :: it, ib1, ib2;

              do it = 1 , 2

                  do ib2 = 1 , size(drho1(it)%nnzblocks,2)
                      do ib1 = 1 , size(drho1(it)%nnzblocks,1)
                          if( drho1(it)%nnzblocks(ib1,ib2) ) deallocate( drho1(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( drho1(it)%blocks    );
                  deallocate( drho1(it)%nnzblocks );

                  do ib2 = 1 , size(drho2(it)%nnzblocks,2)
                      do ib1 = 1 , size(drho2(it)%nnzblocks,1)
                          if( drho2(it)%nnzblocks(ib1,ib2) ) deallocate( drho2(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( drho2(it)%blocks    );
                  deallocate( drho2(it)%nnzblocks );

              end do

              return;
          end subroutine dealloc_drho

      end module drho

      module dkappa
          use dataTypes;
          implicit none;

          type(complexBlockMatrix) , dimension(2) :: dkappa1_pl;
          type(complexBlockMatrix) , dimension(2) :: dkappa1_mi;

      contains

          subroutine alloc_dkappa()
              use fam_input;
              use simplex;
              implicit none;
              integer   :: it, ib1, ib2, i1, i2;
              character :: fg1, fg2;
              integer   :: ml1, ml2;


              do it = 1 , 2

                  allocate( dkappa1_pl(it)%nnzblocks(1:N_blocks,1:N_blocks) );
                  allocate( dkappa1_mi(it)%nnzblocks(1:N_blocks,1:N_blocks) );

                  dkappa1_pl(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  dkappa1_mi(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)

                                  fg1 = fg_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then
                                      dkappa1_pl(it)%nnzblocks(ib1,ib2) = .true.;
                                      dkappa1_mi(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if
                                  if( fg1/=fg2 .and. ( abs(ml1-ml2)==K_multipole .or. abs(ml1+ml2+1)==K_multipole ) ) then
                                      dkappa1_pl(it)%nnzblocks(ib1,ib2) = .true.;
                                      dkappa1_mi(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if

                              end do
                          end do
                      end do
                  end do


                  allocate( dkappa1_pl(it)%blocks(1:N_blocks,1:N_blocks) );
                  allocate( dkappa1_mi(it)%blocks(1:N_blocks,1:N_blocks) );
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks

                          if( dkappa1_pl(it)%nnzblocks(ib1,ib2) ) &
                              allocate( dkappa1_pl(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                          if( dkappa1_mi(it)%nnzblocks(ib1,ib2) ) &
                              allocate( dkappa1_mi(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                      end do
                  end do

              end do

              return;
          end subroutine alloc_dkappa

          subroutine dealloc_dkappa()
              implicit none;
              integer :: it, ib1, ib2;

              do it = 1 , 2

                  do ib2 = 1 , size(dkappa1_pl(it)%nnzblocks,2)
                      do ib1 = 1 , size(dkappa1_pl(it)%nnzblocks,1)
                          if( dkappa1_pl(it)%nnzblocks(ib1,ib2) ) deallocate( dkappa1_pl(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( dkappa1_pl(it)%blocks    );
                  deallocate( dkappa1_pl(it)%nnzblocks );

                  do ib2 = 1 , size(dkappa1_mi(it)%nnzblocks,2)
                      do ib1 = 1 , size(dkappa1_mi(it)%nnzblocks,1)
                          if( dkappa1_mi(it)%nnzblocks(ib1,ib2) ) deallocate( dkappa1_mi(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( dkappa1_mi(it)%blocks    );
                  deallocate( dkappa1_mi(it)%nnzblocks );

              end do

              return;
          end subroutine dealloc_dkappa

      end module dkappa

      module dcurr
          implicit none;

          double complex , dimension(:,:,:) , allocatable :: dj_r;
          double complex , dimension(:,:,:) , allocatable :: dj_p;
          double complex , dimension(:,:,:) , allocatable :: dj_z;
          double complex , dimension(:,:,:) , allocatable :: dj_1;
          double complex , dimension(:,:,:) , allocatable :: dj_2;
          double complex , dimension(:,:,:) , allocatable :: dj_3;

      contains

          subroutine alloc_dcurr()
              use fam_input;
              implicit none;

              allocate( dj_r( -NGH:+NGH , 1:NGL , 1:2 ) );
              allocate( dj_p( -NGH:+NGH , 1:NGL , 1:2 ) );
              allocate( dj_z( -NGH:+NGH , 1:NGL , 1:2 ) );
              allocate( dj_1( -NGH:+NGH , 1:NGL , 1:2 ) );
              allocate( dj_2( -NGH:+NGH , 1:NGL , 1:2 ) );
              allocate( dj_3( -NGH:+NGH , 1:NGL , 1:2 ) );

              return;
          end subroutine alloc_dcurr

          subroutine dealloc_dcurr()
              implicit none;

              deallocate( dj_r );
              deallocate( dj_p );
              deallocate( dj_z );
              deallocate( dj_1 );
              deallocate( dj_2 );
              deallocate( dj_3 );

              return;
          end subroutine dealloc_dcurr

      end module dcurr

      module dlaplace
          implicit none;

          double complex , dimension(:,:) , allocatable :: ldrho_vprot;
          double complex , dimension(:,:) , allocatable :: ldrho_s;
          double complex , dimension(:,:) , allocatable :: ldj_1prot;
          double complex , dimension(:,:) , allocatable :: ldj_2prot;
          double complex , dimension(:,:) , allocatable :: ldj_3prot;

      contains

          subroutine alloc_dlaplace()
              use fam_input;
              implicit none;

              allocate( ldrho_vprot( -NGH:+NGH , 1:NGL ) );
              allocate(     ldrho_s( -NGH:+NGH , 1:NGL ) );
              allocate(   ldj_1prot( -NGH:+NGH , 1:NGL ) );
              allocate(   ldj_2prot( -NGH:+NGH , 1:NGL ) );
              allocate(   ldj_3prot( -NGH:+NGH , 1:NGL ) );

              return;
          end subroutine alloc_dlaplace

          subroutine dealloc_dlaplace()
              implicit none;

              deallocate( ldrho_vprot );
              deallocate( ldrho_s     );
              deallocate( ldj_1prot   );
              deallocate( ldj_2prot   );
              deallocate( ldj_3prot   );

              return;
          end subroutine dealloc_dlaplace

      end module dlaplace

      module famlaplacianmod
          implicit none;

          integer                                           :: Nshells;
          integer          , dimension(:)     , allocatable :: NZZ;
          integer          , dimension(:)     , allocatable :: NRR;
          double complex   , dimension(:)     , allocatable :: c;
          double precision , dimension(:,:)   , allocatable :: phiz;
          double precision , dimension(:,:,:) , allocatable :: phirK;

      contains

          subroutine alloc_famlaplacianmod()
              use fam_input;
              implicit none;
              integer :: nsize;
              integer :: NzMax;
              integer :: NrMax;
              integer :: K_Max;

              ! Linear combinations of: phiz(nz1,z,bz) * phiz(nz2,z,bz) * phir(nr1,ml1,r,bp) * phir(nr2,ml2,r,bp),
              ! for nz1+2*nr1+|ml1|<=n0f+1 and nz2+2*nr2+|ml2|<=n0f+1, and |ml1-ml2|=K, can be written as
              ! linear combinations of phiz(Nz,z,bz/sqrt(2)) * phir(Nr,|K|,r,bp/sqrt(2)),
              ! for Nz+2*Nr+|K|<=2*(n0f+1), or with potentially adding zeros for Nz+2*Nr<=2*(n0f+1).
              !
              ! Linear combinations of: phiz(nz1,z,bz) * phiz(nz2,z,bz) * phir(nr1,ml1,r,bp) * phir(nr2,ml2,r,bp),
              ! for nz1+2*nr1+|ml1|<=n0f+1 and nz2+2*nr2+|ml2|<=n0f+1, and ml1+ml2+1=+K, can be written as
              ! linear combinations of phiz(Nz,z,bz/sqrt(2)) * phir(Nr,|K-1|,r,bp/sqrt(2)),
              ! for Nz+2*Nr+|K-1|<=2*(n0f+1), or with potentially adding zeros for Nz+2*Nr<=2*(n0f+1).
              !
              ! Linear combinations of: phiz(nz1,z,bz) * phiz(nz2,z,bz) * phir(nr1,ml1,r,bp) * phir(nr2,ml2,r,bp),
              ! for nz1+2*nr1+|ml1|<=n0f+1 and nz2+2*nr2+|ml2|<=n0f+1, and ml1+ml2+1=-K, can be written as
              ! linear combinations of phiz(Nz,z,bz/sqrt(2)) * phir(Nr,|K+1|,r,bp/sqrt(2)),
              ! for Nz+2*Nr+|K+1|<=2*(n0f+1), or with potentially adding zeros for Nz+2*Nr<=2*(n0f+1).

              Nshells = 2*(n0f+1);
              ! Number of (Nz>=0,Nr>=0) pairs such that Nz + 2*Nr <= Nshells.
              nsize = ( (Nshells+1)*(Nshells+3) + 1 - mod(Nshells,2) )/4;

              allocate( NZZ(1:nsize) );
              allocate( NRR(1:nsize) );
              allocate(   c(1:nsize) );

              NzMax = Nshells;
              NrMax = Nshells/2;
              K_Max = abs(J_MAX+1);
              allocate( phiz ( -NGH:+NGH  , 0:NzMax           ) );
              allocate( phirK(     1:NGL  , 0:NrMax , 0:K_Max ) );

              return;
          end subroutine alloc_famlaplacianmod

          subroutine dealloc_famlaplacianmod()
              implicit none;

              deallocate( NZZ   );
              deallocate( NRR   );
              deallocate( c     );
              deallocate( phiz  );
              deallocate( phirK );

              return;
          end subroutine dealloc_famlaplacianmod

      end module famlaplacianmod

      module dmesons
          implicit none;

          double complex , dimension(:,:) , allocatable :: dsigma;
          double complex , dimension(:,:) , allocatable :: domega_0;
          double complex , dimension(:,:) , allocatable :: domega_r;
          double complex , dimension(:,:) , allocatable :: domega_p;
          double complex , dimension(:,:) , allocatable :: domega_z;
          double complex , dimension(:,:) , allocatable :: drho_0;
          double complex , dimension(:,:) , allocatable :: drho_r;
          double complex , dimension(:,:) , allocatable :: drho_p;
          double complex , dimension(:,:) , allocatable :: drho_z;

      contains

          subroutine alloc_dmesons()
              use fam_input;
              implicit none;

              allocate( dsigma  ( -NGH:+NGH , 1:NGL ) );
              allocate( domega_0( -NGH:+NGH , 1:NGL ) );
              allocate( domega_r( -NGH:+NGH , 1:NGL ) );
              allocate( domega_p( -NGH:+NGH , 1:NGL ) );
              allocate( domega_z( -NGH:+NGH , 1:NGL ) );
              allocate( drho_0  ( -NGH:+NGH , 1:NGL ) );
              allocate( drho_r  ( -NGH:+NGH , 1:NGL ) );
              allocate( drho_p  ( -NGH:+NGH , 1:NGL ) );
              allocate( drho_z  ( -NGH:+NGH , 1:NGL ) );

              return;
          end subroutine alloc_dmesons

          subroutine dealloc_dmesons()
              implicit none;

              deallocate( dsigma   );
              deallocate( domega_0 );
              deallocate( domega_r );
              deallocate( domega_p );
              deallocate( domega_z );
              deallocate( drho_0   );
              deallocate( drho_r   );
              deallocate( drho_p   );
              deallocate( drho_z   );

              return;
          end subroutine dealloc_dmesons

      end module dmesons

      module dcoulomb
          implicit none;

          double complex , dimension(:,:) , allocatable :: dVCoulomb_0;
          double complex , dimension(:,:) , allocatable :: dVCoulomb_r;
          double complex , dimension(:,:) , allocatable :: dVCoulomb_p;
          double complex , dimension(:,:) , allocatable :: dVCoulomb_z;

      contains

          subroutine alloc_dcoulomb()
              use fam_input;
              implicit none;

              allocate( dVCoulomb_0( -NGH:+NGH , 1:NGL ) );
              allocate( dVCoulomb_r( -NGH:+NGH , 1:NGL ) );
              allocate( dVCoulomb_p( -NGH:+NGH , 1:NGL ) );
              allocate( dVCoulomb_z( -NGH:+NGH , 1:NGL ) );

              return;
          end subroutine alloc_dcoulomb

          subroutine dealloc_dcoulomb()
              implicit none;

              deallocate( dVCoulomb_0 );
              deallocate( dVCoulomb_r );
              deallocate( dVCoulomb_p );
              deallocate( dVCoulomb_z );

              return;
          end subroutine dealloc_dcoulomb

      end module dcoulomb

      module dpotentials
          implicit none;

          double complex , dimension(:,:,:) , allocatable :: dVplusS;
          double complex , dimension(:,:,:) , allocatable :: dVminusS;
          double complex , dimension(:,:,:) , allocatable :: dSigma_z;
          double complex , dimension(:,:,:) , allocatable :: dSigma_r;
          double complex , dimension(:,:,:) , allocatable :: dSigma_p;

      contains

          subroutine alloc_dpotentials()
              use fam_input;
              implicit none;

              allocate( dVplusS ( -NGH:+NGH , 1:NGL , 1:2 ) );
              allocate( dVminusS( -NGH:+NGH , 1:NGL , 1:2 ) );
              allocate( dSigma_z( -NGH:+NGH , 1:NGL , 1:2 ) );
              allocate( dSigma_r( -NGH:+NGH , 1:NGL , 1:2 ) );
              allocate( dSigma_p( -NGH:+NGH , 1:NGL , 1:2 ) );

              return;
          end subroutine alloc_dpotentials

          subroutine dealloc_dpotentials()
              implicit none;

              deallocate( dVplusS  );
              deallocate( dVminusS );
              deallocate( dSigma_z );
              deallocate( dSigma_r );
              deallocate( dSigma_p );

              return;
          end subroutine dealloc_dpotentials

      end module dpotentials

      module nuclearlocfunc
          implicit none;

          double precision , dimension(:,:,:) , allocatable :: C0;
          double complex   , dimension(:,:,:) , allocatable :: dC;

      contains

          subroutine alloc_nuclearlocfunc()
              use fam_input;
              implicit none;

              allocate( C0( -NGH:+NGH , 1:NGL , 1:2 ) );
              allocate( dC( -NGH:+NGH , 1:NGL , 1:2 ) );

              return;
          end subroutine alloc_nuclearlocfunc

          subroutine dealloc_nuclearlocfunc()
              implicit none;

              deallocate( C0 );
              deallocate( dC );

              return;
          end subroutine dealloc_nuclearlocfunc

      end module nuclearlocfunc

      module gmres
          implicit none;

          integer                                       :: xsize;
          integer                                       :: maxIter;
          double complex , dimension(:,:) , allocatable :: Q;
          double complex , dimension(:)   , allocatable :: b, r0, v, x_solution;
          double complex , dimension(:,:) , allocatable :: H, H_tmp;
          double complex , dimension(:)   , allocatable :: beta, beta_tmp, y;

      contains

          subroutine alloc_gmres()
              use fam_input;
              use simplex;
              implicit none;
              integer   :: it, ib1, ib2, i1, i2;
              character :: fg1, fg2;
              integer   :: ml1, ml2;

              ! Determining the dimension of Arnoldi vectors.
              xsize = 0;
              ! Nonzero elements of the upper triangle of dh1 matrix.
              do it = 1 , 2
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)
                                  if( (ib1<ib2) .or. ( (ib1==ib2) .and. (i1<=i2) ) ) then ! Loop over the upper triangle.

                                      fg1 = fg_spx(ib1)%index(i1);
                                      ml1 = ml_spx(ib1)%index(i1);

                                      fg2 = fg_spx(ib2)%index(i2);
                                      ml2 = ml_spx(ib2)%index(i2);

                                      if( fg1==fg2 .and. abs(ml1-ml2)   == K_multipole ) xsize = xsize + 1;
                                      if( fg1/=fg2 .and. abs(ml1-ml2)   == K_multipole ) xsize = xsize + 1;
                                      if( fg1/=fg2 .and. abs(ml1+ml2+1) == K_multipole ) xsize = xsize + 1;

                                  end if
                              end do
                          end do
                      end do
                  end do
              end do
              ! Nonzero elements of the upper triangles of dDelta1_pl and dDelta1_mi matrices.
              do it = 1 , 2
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)
                                  if( (ib1<ib2) .or. ( (ib1==ib2) .and. (i1<=i2) ) ) then ! Loop over the upper triangle.

                                      fg1 = fg_spx(ib1)%index(i1);
                                      ml1 = ml_spx(ib1)%index(i1);

                                      fg2 = fg_spx(ib2)%index(i2);
                                      ml2 = ml_spx(ib2)%index(i2);

                                      if( all([fg1,fg2]==['f','f']) .and. abs(ml1-ml2)==K_multipole ) xsize = xsize + 1; ! dDelta1_pl.
                                      if( all([fg1,fg2]==['f','f']) .and. abs(ml1-ml2)==K_multipole ) xsize = xsize + 1; ! dDelta1_mi.

                                  end if
                              end do
                          end do
                      end do
                  end do
              end do

              maxIter = NoArnoldiVectors;


              allocate(          Q( 1:xsize     , 1:maxIter+1 ) );
              allocate(          b( 1:xsize                   ) );
              allocate(         r0( 1:xsize                   ) );
              allocate(          v( 1:xsize                   ) );
              allocate( x_solution( 1:xsize                   ) );
              allocate(          H( 1:maxIter+1 , 1:maxIter   ) );
              allocate(      H_tmp( 1:maxIter+1 , 1:maxIter   ) );
              allocate(       beta( 1:maxIter+1               ) );
              allocate(   beta_tmp( 1:maxIter+1               ) );
              allocate(          y( 1:maxIter+1               ) );

              return;
          end subroutine alloc_gmres

          subroutine dealloc_gmres()
              implicit none;

              deallocate( Q          );
              deallocate( b          );
              deallocate( r0         );
              deallocate( v          );
              deallocate( x_solution );
              deallocate( H          );
              deallocate( H_tmp      );
              deallocate( beta       );
              deallocate( beta_tmp   );
              deallocate( y          );

              return;
          end subroutine dealloc_gmres

          subroutine extractFromGMRESvector( xin , dh1 , dh2 , dDelta1_pl , dDelta1_mi )
          use dataTypes;
          use fam_input;
          use simplex;
          implicit none;
          double complex           , dimension(:) , intent(in)    :: xin;
          type(complexBlockMatrix) , dimension(2) , intent(inout) :: dh1;
          type(complexBlockMatrix) , dimension(2) , intent(inout) :: dh2;
          type(complexBlockMatrix) , dimension(2) , intent(inout) :: dDelta1_pl;
          type(complexBlockMatrix) , dimension(2) , intent(inout) :: dDelta1_mi;

          integer   :: ix;
          integer   :: it, ib1, ib2, i1, i2;
          character :: fg1, fg2;
          integer   :: ml1, ml2;

          ix = 0;

          do it = 1 , 2

              do ib2 = 1 , N_blocks
                  do ib1 = 1 , ib2
                      if( dh1(it)%nnzblocks(ib1,ib2) ) then

                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , merge( i2 , id_spx(ib1) , ib1==ib2 )

                                  fg1 = fg_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then
                                      ix = ix + 1;
                                      dh1(it)%blocks(ib1,ib2)%mat(i1,i2) = xin(ix);
                                  end if

                                  if( fg1/=fg2 .and. abs(ml1-ml2)==K_multipole ) then
                                      ix = ix + 1;
                                      dh1(it)%blocks(ib1,ib2)%mat(i1,i2) = xin(ix);
                                  end if

                                  if( fg1/=fg2 .and. abs(ml1+ml2+1)==K_multipole ) then
                                      ix = ix + 1;
                                      dh1(it)%blocks(ib1,ib2)%mat(i1,i2) = xin(ix);
                                  end if

                              end do
                          end do

                      end if
                  end do
              end do

              call constructFulldh1FromItsUpperTriangle( dh1(it) );

              call constructFulldh2FromFulldh1( dh1(it) , dh2(it) );

          end do

          do it = 1 , 2

              do ib2 = 1 , N_blocks
                  do ib1 = 1 , ib2
                      if( dDelta1_pl(it)%nnzblocks(ib1,ib2) .and. dDelta1_mi(it)%nnzblocks(ib1,ib2) ) then

                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , merge( i2 , id_spx(ib1) , ib1==ib2 )

                                  fg1 = fg_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( all([fg1,fg2]==['f','f']) .and. abs(ml1-ml2)==K_multipole ) then
                                      ix = ix + 1;
                                      dDelta1_pl(it)%blocks(ib1,ib2)%mat(i1,i2) = xin(ix);

                                      ix = ix + 1;
                                      dDelta1_mi(it)%blocks(ib1,ib2)%mat(i1,i2) = conjg(xin(ix));
                                  end if

                              end do
                          end do

                      end if
                  end do
              end do

              call constructFulldDelta1FromItsUpperTriangle( dDelta1_pl(it) );

              call constructFulldDelta1FromItsUpperTriangle( dDelta1_mi(it) );

          end do

          call assert( ix == size(xin) , 'x size error.' );

          return;
          end subroutine extractFromGMRESvector

          subroutine insertIntoGMRESvector( dh1 , dh2 , dDelta1_pl , dDelta1_mi , xout )
          use dataTypes;
          use fam_input;
          use simplex;
          implicit none;
          type(complexBlockMatrix) , dimension(2) , intent(in)    :: dh1;
          type(complexBlockMatrix) , dimension(2) , intent(in)    :: dh2;
          type(complexBlockMatrix) , dimension(2) , intent(in)    :: dDelta1_pl;
          type(complexBlockMatrix) , dimension(2) , intent(in)    :: dDelta1_mi;
          double complex           , dimension(:) , intent(inout) :: xout;

          integer   :: ix;
          integer   :: it, ib1, ib2, i1, i2;
          character :: fg1, fg2;
          integer   :: ml1, ml2;

          ix = 0;

          do it = 1 , 2
              do ib2 = 1 , N_blocks
                  do ib1 = 1 , ib2
                      if( dh1(it)%nnzblocks(ib1,ib2) ) then

                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , merge( i2 , id_spx(ib1) , ib1==ib2 )

                                  fg1 = fg_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then
                                      ix = ix + 1;
                                      xout(ix) = dh1(it)%blocks(ib1,ib2)%mat(i1,i2);
                                  end if

                                  if( fg1/=fg2 .and. abs(ml1-ml2)==K_multipole ) then
                                      ix = ix + 1;
                                      xout(ix) = dh1(it)%blocks(ib1,ib2)%mat(i1,i2);
                                  end if

                                  if( fg1/=fg2 .and. abs(ml1+ml2+1)==K_multipole ) then
                                      ix = ix + 1;
                                      xout(ix) = dh1(it)%blocks(ib1,ib2)%mat(i1,i2);
                                  end if

                              end do
                          end do

                      end if
                  end do
              end do
          end do

          do it = 1 , 2
              do ib2 = 1 , N_blocks
                  do ib1 = 1 , ib2
                      if( dDelta1_pl(it)%nnzblocks(ib1,ib2) .and. dDelta1_mi(it)%nnzblocks(ib1,ib2) ) then

                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , merge( i2 , id_spx(ib1) , ib1==ib2 )

                                  fg1 = fg_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( all([fg1,fg2]==['f','f']) .and. abs(ml1-ml2)==K_multipole ) then
                                      ix = ix + 1;
                                      xout(ix) = dDelta1_pl(it)%blocks(ib1,ib2)%mat(i1,i2);

                                      ix = ix + 1;
                                      xout(ix) = conjg( dDelta1_mi(it)%blocks(ib1,ib2)%mat(i1,i2) );
                                  end if

                              end do
                          end do

                      end if
                  end do
              end do
          end do

          call assert( ix == size(xout) , 'x size error.' );

          return;
          end subroutine insertIntoGMRESvector

      end module gmres

      module KPMdata
          use dataTypes;
          implicit none;

          type(complexBlockMatrix) , dimension(2) :: alpha_old_x;
          type(complexBlockMatrix) , dimension(2) :: alpha_old_y;
          type(complexBlockMatrix) , dimension(2) :: alpha_new_x;
          type(complexBlockMatrix) , dimension(2) :: alpha_new_y;
          type(complexBlockMatrix) , dimension(2) :: alpha_tmp_x;
          type(complexBlockMatrix) , dimension(2) :: alpha_tmp_y;

      contains

          subroutine alloc_KPMdata()
              use fam_input;
              use simplex;
              implicit none;
              integer   :: it, ib1, ib2, i1, i2;
              character :: fg1, fg2;
              integer   :: ml1, ml2;

              do it = 1 , 2

                  allocate( alpha_old_x(it)%nnzblocks(1:N_blocks,1:N_blocks) );
                  allocate( alpha_old_y(it)%nnzblocks(1:N_blocks,1:N_blocks) );
                  allocate( alpha_new_x(it)%nnzblocks(1:N_blocks,1:N_blocks) );
                  allocate( alpha_new_y(it)%nnzblocks(1:N_blocks,1:N_blocks) );
                  allocate( alpha_tmp_x(it)%nnzblocks(1:N_blocks,1:N_blocks) );
                  allocate( alpha_tmp_y(it)%nnzblocks(1:N_blocks,1:N_blocks) );

                  alpha_old_x(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  alpha_old_y(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  alpha_new_x(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  alpha_new_y(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  alpha_tmp_x(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  alpha_tmp_y(it)%nnzblocks(1:N_blocks,1:N_blocks) = .false.;
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)

                                  fg1 = fg_spx(ib1)%index(i1);
                                  ml1 = ml_spx(ib1)%index(i1);

                                  fg2 = fg_spx(ib2)%index(i2);
                                  ml2 = ml_spx(ib2)%index(i2);

                                  if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then
                                      alpha_old_x(it)%nnzblocks(ib1,ib2) = .true.;
                                      alpha_old_y(it)%nnzblocks(ib1,ib2) = .true.;
                                      alpha_new_x(it)%nnzblocks(ib1,ib2) = .true.;
                                      alpha_new_y(it)%nnzblocks(ib1,ib2) = .true.;
                                      alpha_tmp_x(it)%nnzblocks(ib1,ib2) = .true.;
                                      alpha_tmp_y(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if
                                  if( fg1/=fg2 .and. ( abs(ml1-ml2)==K_multipole .or. abs(ml1+ml2+1)==K_multipole ) ) then
                                      alpha_old_x(it)%nnzblocks(ib1,ib2) = .true.;
                                      alpha_old_y(it)%nnzblocks(ib1,ib2) = .true.;
                                      alpha_new_x(it)%nnzblocks(ib1,ib2) = .true.;
                                      alpha_new_y(it)%nnzblocks(ib1,ib2) = .true.;
                                      alpha_tmp_x(it)%nnzblocks(ib1,ib2) = .true.;
                                      alpha_tmp_y(it)%nnzblocks(ib1,ib2) = .true.;
                                  end if

                              end do
                          end do
                      end do
                  end do


                  allocate( alpha_old_x(it)%blocks(1:N_blocks,1:N_blocks) );
                  allocate( alpha_old_y(it)%blocks(1:N_blocks,1:N_blocks) );
                  allocate( alpha_new_x(it)%blocks(1:N_blocks,1:N_blocks) );
                  allocate( alpha_new_y(it)%blocks(1:N_blocks,1:N_blocks) );
                  allocate( alpha_tmp_x(it)%blocks(1:N_blocks,1:N_blocks) );
                  allocate( alpha_tmp_y(it)%blocks(1:N_blocks,1:N_blocks) );
                  do ib2 = 1 , N_blocks
                      do ib1 = 1 , N_blocks

                          if( alpha_old_x(it)%nnzblocks(ib1,ib2) ) &
                              allocate( alpha_old_x(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                          if( alpha_old_y(it)%nnzblocks(ib1,ib2) ) &
                              allocate( alpha_old_y(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                          if( alpha_new_x(it)%nnzblocks(ib1,ib2) ) &
                              allocate( alpha_new_x(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                          if( alpha_new_y(it)%nnzblocks(ib1,ib2) ) &
                              allocate( alpha_new_y(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                          if( alpha_tmp_x(it)%nnzblocks(ib1,ib2) ) &
                              allocate( alpha_tmp_x(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );

                          if( alpha_tmp_y(it)%nnzblocks(ib1,ib2) ) &
                              allocate( alpha_tmp_y(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) );
                      end do
                  end do

              end do

              return;
          end subroutine alloc_KPMdata

          subroutine dealloc_KPMdata()
              implicit none;
              integer :: it, ib1, ib2;

              do it = 1 , 2

                  do ib2 = 1 , size(alpha_old_x(it)%nnzblocks,2)
                      do ib1 = 1 , size(alpha_old_x(it)%nnzblocks,1)
                          if( alpha_old_x(it)%nnzblocks(ib1,ib2) ) deallocate( alpha_old_x(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( alpha_old_x(it)%blocks    );
                  deallocate( alpha_old_x(it)%nnzblocks );

                  do ib2 = 1 , size(alpha_old_y(it)%nnzblocks,2)
                      do ib1 = 1 , size(alpha_old_y(it)%nnzblocks,1)
                          if( alpha_old_y(it)%nnzblocks(ib1,ib2) ) deallocate( alpha_old_y(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( alpha_old_y(it)%blocks    );
                  deallocate( alpha_old_y(it)%nnzblocks );

                  do ib2 = 1 , size(alpha_new_x(it)%nnzblocks,2)
                      do ib1 = 1 , size(alpha_new_x(it)%nnzblocks,1)
                          if( alpha_new_x(it)%nnzblocks(ib1,ib2) ) deallocate( alpha_new_x(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( alpha_new_x(it)%blocks    );
                  deallocate( alpha_new_x(it)%nnzblocks );

                  do ib2 = 1 , size(alpha_new_y(it)%nnzblocks,2)
                      do ib1 = 1 , size(alpha_new_y(it)%nnzblocks,1)
                          if( alpha_new_y(it)%nnzblocks(ib1,ib2) ) deallocate( alpha_new_y(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( alpha_new_y(it)%blocks    );
                  deallocate( alpha_new_y(it)%nnzblocks );

                  do ib2 = 1 , size(alpha_tmp_x(it)%nnzblocks,2)
                      do ib1 = 1 , size(alpha_tmp_x(it)%nnzblocks,1)
                          if( alpha_tmp_x(it)%nnzblocks(ib1,ib2) ) deallocate( alpha_tmp_x(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( alpha_tmp_x(it)%blocks    );
                  deallocate( alpha_tmp_x(it)%nnzblocks );

                  do ib2 = 1 , size(alpha_tmp_y(it)%nnzblocks,2)
                      do ib1 = 1 , size(alpha_tmp_y(it)%nnzblocks,1)
                          if( alpha_tmp_y(it)%nnzblocks(ib1,ib2) ) deallocate( alpha_tmp_y(it)%blocks(ib1,ib2)%mat );
                      end do
                  end do
                  deallocate( alpha_tmp_y(it)%blocks    );
                  deallocate( alpha_tmp_y(it)%nnzblocks );

              end do

              return;
          end subroutine dealloc_KPMdata

      end module KPMdata

      module clock_module
          implicit none;

              type , public :: clockClass
                  integer , private :: clock_start;
                  integer , private :: clock_end;
                  integer , private :: clock_rate;
              contains
                  procedure , public :: tic     => clock_tic;
                  procedure , public :: toc     => clock_toc;
                  procedure , public :: getTime => clock_getTime;
              end type clockClass


              private :: clock_tic;
              private :: clock_toc;
              private :: clock_getTime;

          contains

              subroutine clock_tic(this)
                  implicit none;
                  class(clockClass) , intent(inout) :: this;
                      call system_clock( this%clock_start , this%clock_rate );
                  return;
              end subroutine clock_tic

              subroutine clock_toc(this)
                  implicit none;
                  class(clockClass) , intent(inout) :: this;
                      call system_clock( this%clock_end , this%clock_rate );
                  return;
              end subroutine clock_toc

              double precision function clock_getTime(this) result(time)
                  implicit none;
                  class(clockClass) , intent(in) :: this;
                      time = real( this%clock_end - this%clock_start , kind=8 ) / this%clock_rate;
                  return;
              end function clock_getTime

      end module clock_module
