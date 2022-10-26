!======================================================================!

      subroutine init_mesons()

!======================================================================!
      use dataTypes;
      use fam_input;
      use ddpc1ddme2;
      use basis;
      use quadrature;
      use mesmat;
      implicit none;

      integer                                                   :: K;
      integer                                                   :: nz,  nr;
      integer                                                   :: nz1, nz2;
      integer                                                   :: nr1, nr2;
      integer                                                   :: i, i1, i2;
      integer                                                   :: ihl, ih, il;
      double precision                                          :: z, r;
      double precision                                          :: bp_mes, bz_mes;
      type(integer1Darray) , dimension(0:J_MAX+1)               :: nz_mes;
      type(integer1Darray) , dimension(0:J_MAX+1)               :: nr_mes;
      type(real2Darray)    , dimension(0:J_MAX+1)               :: H;
      type(real2Darray)    , dimension(0:J_MAX+1)               :: phiz;
      type(real2Darray)    , dimension(0:J_MAX+1)               :: phirK;
      double precision     , dimension(:,:)       , allocatable :: L;
      double precision     , parameter                          :: hbc = 197.328284d0;
      double precision     , parameter                          :: one = 1.d0;
      integer                                                   :: INFO;
      double precision     , external                           :: phi_nz;
      double precision     , external                           :: phi_nr_ml;



      if( LagrangianModelName /= 'DD-ME2' ) then
          return;
      end if


      ! Construction of (nz>=0,nr>=0) pairs such that nz + 2*nr + K <= n0b.
      do K = 0 , J_MAX+1

          allocate( nz_mes(K)%index(1:nP(K)) );
          allocate( nr_mes(K)%index(1:nP(K)) );

          i = 0;
          do nr = 0 , (n0b-K)/2
              do nz = 0 , n0b-K-2*nr
                  i = i + 1;

                  nz_mes(K)%index(i) = nz;
                  nr_mes(K)%index(i) = nr;

              end do
          end do
          call assert( i == nP(K) , 'nP(K) wrong!' );

      end do


      ! Construction of H^{K} matrices, but without the meson's mass squared on the diagonal.
      do K = 0 , J_MAX+1

          allocate( H(K)%mat( 1:nP(K) , 1:nP(K) ) );

          bp_mes = bp/sqrt(2.d0);
          bz_mes = bz/sqrt(2.d0);

          H(K)%mat(:,:) = 0.d0;
          do i2 = 1 , nP(K)
              do i1 = 1 , nP(K)

                  nz1 = nz_mes(K)%index(i1);
                  nr1 = nr_mes(K)%index(i1);

                  nz2 = nz_mes(K)%index(i2);
                  nr2 = nr_mes(K)%index(i2);

                  if( nz1==nz2   .and. nr1==nr2   ) H(K)%mat(i1,i2) = H(K)%mat(i1,i2) + ( + (nz1+0.5d0)/bz_mes**2 + (2*nr1+K+1)/bp_mes**2  );
                  if( nz1==nz2   .and. nr1==nr2+1 ) H(K)%mat(i1,i2) = H(K)%mat(i1,i2) + ( + sqrt( 1.d0 * nr1*(nr1+K)     ) / bp_mes**2     );
                  if( nz1==nz2   .and. nr2==nr1+1 ) H(K)%mat(i1,i2) = H(K)%mat(i1,i2) + ( + sqrt( 1.d0 * nr2*(nr2+K)     ) / bp_mes**2     );
                  if( nz1==nz2+2 .and. nr1==nr2   ) H(K)%mat(i1,i2) = H(K)%mat(i1,i2) + ( - sqrt( 1.d0 * (nz2+1)*(nz2+2) ) / (2*bz_mes**2) );
                  if( nz2==nz1+2 .and. nr1==nr2   ) H(K)%mat(i1,i2) = H(K)%mat(i1,i2) + ( - sqrt( 1.d0 * (nz1+1)*(nz1+2) ) / (2*bz_mes**2) );

              end do
          end do

      end do


      ! Calculation of phiz and phirK.
      do K = 0 , J_MAX+1

          allocate(  phiz(K)%mat( -NGH:+NGH , 0:maxval(nz_mes(K)%index(:)) ) );
          allocate( phirK(K)%mat(    1:NGL  , 0:maxval(nr_mes(K)%index(:)) ) );

          do nz = 0 , maxval(nz_mes(K)%index(:))
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;
                  z = sign(1,ih) * zb_fam(abs(ih));
                  phiz(K)%mat(ih,nz) = phi_nz( nz , bz/sqrt(2.d0) , z );
              end do
          end do

          do nr = 0 , maxval(nr_mes(K)%index(:))
              do il = 1 , NGL
                  r = rb_fam(il);
                  phirK(K)%mat(il,nr) = phi_nr_ml( nr , K , bp/sqrt(2.d0) , r );
              end do
          end do

      end do


      ! Calculation of Psig.
      do K = 0 , J_MAX+1

          ! Construct Phi^{K} matrix.
          do i = 1 , nP(K)
              nz = nz_mes(K)%index(i);
              nr = nr_mes(K)%index(i);

              ihl = 0;
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;
                      ihl = ihl + 1;
                      Psig(K)%mat(i,ihl) = phiz(K)%mat(ih,nz) * phirK(K)%mat(il,nr);
                  end do
              end do
              call assert( ihl==size(Psig(K)%mat,2) , 'Error in meson P matrix.' );

          end do
          call assert( nP(K)==size(Psig(K)%mat,1) , 'Error in meson P matrix.' );


          ! Cholesky factorization of H^{K}.
          allocate( L( 1:nP(K) , 1:nP(K) ) );
          do i2 = 1 , nP(K)
              do i1 = 1 , nP(K)
                  L(i1,i2) = H(K)%mat(i1,i2) + merge( (m_sigma/hbc)**2 , 0.d0 , i1==i2 );
              end do
          end do
          call dpotrf( 'L' , nP(K) , L , size(L,1) , INFO );
          call assert( INFO==0 , 'dpotrf() failed!' );
          ! Now the lower triangle of L contains the lower triangular matrix L^{K} such that H^{K} = L^{K} * (L^{K})^T.

          ! Now we want to find P^{K} = (L^{K})^{-1} * Phi^{K}.
          ! Phi^{K} is stored in Psig(K)%mat, while L^{K} is stored in L.
          ! We actually solve the linear system: L^{K} * P^{K} = Phi^{K}.
          call dtrsm( 'L' , 'L' , 'N' , 'N'                     , &
                      size(Psig(K)%mat,1) , size(Psig(K)%mat,2) , &
                      one                                       , &
                      L(1,1)              , size(L,1)           , &
                      Psig(K)%mat(1,1)    , size(Psig(K)%mat,1)   );
          deallocate(L);

          ! Now Psig(K)%mat contains P^{K} = (L^{K})^{-1} * Phi^{K}, i.e. just what we need.

      end do


      ! Calculation of Pome.
      do K = 0 , J_MAX+1

          ! Construct Phi^{K} matrix.
          do i = 1 , nP(K)
              nz = nz_mes(K)%index(i);
              nr = nr_mes(K)%index(i);

              ihl = 0;
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;
                      ihl = ihl + 1;
                      Pome(K)%mat(i,ihl) = phiz(K)%mat(ih,nz) * phirK(K)%mat(il,nr);
                  end do
              end do
              call assert( ihl==size(Pome(K)%mat,2) , 'Error in meson P matrix.' );

          end do
          call assert( nP(K)==size(Pome(K)%mat,1) , 'Error in meson P matrix.' );


          ! Cholesky factorization of H^{K}.
          allocate( L( 1:nP(K) , 1:nP(K) ) );
          do i2 = 1 , nP(K)
              do i1 = 1 , nP(K)
                  L(i1,i2) = H(K)%mat(i1,i2) + merge( (m_omega/hbc)**2 , 0.d0 , i1==i2 );
              end do
          end do
          call dpotrf( 'L' , nP(K) , L , size(L,1) , INFO );
          call assert( INFO==0 , 'dpotrf() failed!' );
          ! Now the lower triangle of L contains the lower triangular matrix L^{K} such that H^{K} = L^{K} * (L^{K})^T.

          ! Now we want to find P^{K} = (L^{K})^{-1} * Phi^{K}.
          ! Phi^{K} is stored in Pome(K)%mat, while L^{K} is stored in L.
          ! We actually solve the linear system: L^{K} * P^{K} = Phi^{K}.
          call dtrsm( 'L' , 'L' , 'N' , 'N'                     , &
                      size(Pome(K)%mat,1) , size(Pome(K)%mat,2) , &
                      one                                       , &
                      L(1,1)              , size(L,1)           , &
                      Pome(K)%mat(1,1)    , size(Pome(K)%mat,1)   );
          deallocate(L);

          ! Now Pome(K)%mat contains P^{K} = (L^{K})^{-1} * Phi^{K}, i.e. just what we need.

      end do


      ! Calculation of Prho.
      do K = 0 , J_MAX+1

          ! Construct Phi^{K} matrix.
          do i = 1 , nP(K)
              nz = nz_mes(K)%index(i);
              nr = nr_mes(K)%index(i);

              ihl = 0;
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;
                      ihl = ihl + 1;
                      Prho(K)%mat(i,ihl) = phiz(K)%mat(ih,nz) * phirK(K)%mat(il,nr);
                  end do
              end do
              call assert( ihl==size(Prho(K)%mat,2) , 'Error in meson P matrix.' );

          end do
          call assert( nP(K)==size(Prho(K)%mat,1) , 'Error in meson P matrix.' );


          ! Cholesky factorization of H^{K}.
          allocate( L( 1:nP(K) , 1:nP(K) ) );
          do i2 = 1 , nP(K)
              do i1 = 1 , nP(K)
                  L(i1,i2) = H(K)%mat(i1,i2) + merge( (m_rho/hbc)**2 , 0.d0 , i1==i2 );
              end do
          end do
          call dpotrf( 'L' , nP(K) , L , size(L,1) , INFO );
          call assert( INFO==0 , 'dpotrf() failed!' );
          ! Now the lower triangle of L contains the lower triangular matrix L^{K} such that H^{K} = L^{K} * (L^{K})^T.

          ! Now we want to find P^{K} = (L^{K})^{-1} * Phi^{K}.
          ! Phi^{K} is stored in Prho(K)%mat, while L^{K} is stored in L.
          ! We actually solve the linear system: L^{K} * P^{K} = Phi^{K}.
          call dtrsm( 'L' , 'L' , 'N' , 'N'                     , &
                      size(Prho(K)%mat,1) , size(Prho(K)%mat,2) , &
                      one                                       , &
                      L(1,1)              , size(L,1)           , &
                      Prho(K)%mat(1,1)    , size(Prho(K)%mat,1)   );
          deallocate(L);

          ! Now Prho(K)%mat contains P^{K} = (L^{K})^{-1} * Phi^{K}, i.e. just what we need.

      end do


      ! Deallocation of local memory.
      do K = 0 , J_MAX+1
          deallocate( nz_mes(K)%index );
          deallocate( nr_mes(K)%index );

          deallocate( H(K)%mat );

          deallocate(  phiz(K)%mat );
          deallocate( phirK(K)%mat );
      end do




      return;
      end subroutine init_mesons
