!======================================================================!

      subroutine init_basis()

!======================================================================!
      use fam_input;
      use simplex;
      use quadrature;
      use basis;
      use wbasis;
      use PHI;
      use module_SVD;
      implicit none;

      double precision                                :: hbarw0;
      double precision                                :: hbarc;
      double precision                                :: nucleon_mass;
      double precision ,                  parameter   :: pi = 4.d0*atan(1.d0);
      integer                                         :: ib, i, ih, il, ihl;
      integer                                         :: nz, nr, ml;
      double precision                                :: z, r;
      double precision                                :: alpha;
      double precision                                :: xgh(1:2*NGH);
      double precision                                :: wgh(1:2*NGH);
      double precision                                :: xgl(1:NGL);
      double precision                                :: wgl(1:NGL);
      double precision , dimension(:,:) , allocatable :: PHI_tmp, U_PHI_tmp, V_PHI_tmp;
      double precision ,                  external    ::   phi_nz;
      double precision ,                  external    :: d_phi_nz;
      double precision ,                  external    ::   phi_nr_ml;
      double precision ,                  external    :: d_phi_nr_ml;




          ! Calculating oscillator length parameters bz and bp.
          ! bz = sqrt(hbar/(m*w0)) * exp(        sqrt(5/(16*pi)) * beta0 ), for hbar*w0 = 41*A^(-1/3) [MeV].
          ! bp = sqrt(hbar/(m*w0)) * exp( -1/2 * sqrt(5/(16*pi)) * beta0 ), for hbar*w0 = 41*A^(-1/3) [MeV].
          hbarw0       = 41 * (nucleusZ+nucleusN)**(-1.d0/3.d0); ! [MeV].
          hbarc        = 197.328284d0;                           ! [MeV*fm].
          nucleon_mass = 939.d0                                  ! [MeV/c^2].
          bz = sqrt( hbarc**2/(nucleon_mass*hbarw0) ) * exp(          sqrt(5/(16*pi)) * beta0 );
          bp = sqrt( hbarc**2/(nucleon_mass*hbarw0) ) * exp( -0.5d0 * sqrt(5/(16*pi)) * beta0 );




          ! Calculating Gaussian quadrature mesh and weights.
          call gauher( xgh , wgh , 2*NGH );
          do ih = 1 , NGH
              zb_fam(ih) = bz * xgh(1+NGH-ih);
                  wz(ih) = bz * wgh(1+NGH-ih) * exp(xgh(1+NGH-ih)**2);
          end do

          alpha = 0.d0;
          call gaulag( xgl , wgl , NGL , alpha );
          do il = 1 , NGL
              rb_fam(il) = bp * sqrt(xgl(il));
                  wr(il) = 0.5d0 * bp**2 * wgl(il) * exp(xgl(il)) / xgl(il)**alpha;
          end do

          do il = 1 , NGL
              do ih = 1 , NGH
                  wzwr(ih,il) = wz(ih)*wr(il);
              end do
          end do




          ! Calculating phi_z.
          do ib = 1 , N_blocks
              do i = 1 , id_spx(ib)
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;

                      nz = nz_spx(ib)%index(i);
                      z  = sign(1,ih) * zb_fam(abs(ih));

                      phi_z(ib)%mat(ih,i) = phi_nz( nz , bz , z );

                  end do
              end do
          end do

          ! Calculating dphi_z.
          do ib = 1 , N_blocks
              do i = 1 , id_spx(ib)
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;

                      nz = nz_spx(ib)%index(i);
                      z  = sign(1,ih) * zb_fam(abs(ih));

                      dphi_z(ib)%mat(ih,i) = d_phi_nz( nz , bz , z );

                  end do
              end do
          end do

          ! Calculating phi_r.
          do ib = 1 , N_blocks
              do i = 1 , id_spx(ib)
                  do il = 1 , NGL

                      nr = nr_spx(ib)%index(i);
                      ml = ml_spx(ib)%index(i);
                      r  = rb_fam(il);

                      phi_r(ib)%mat(il,i) = phi_nr_ml( nr , abs(ml) , bp , r );

                  end do
              end do
          end do

          ! Calculating dphi_r.
          do ib = 1 , N_blocks
              do i = 1 , id_spx(ib)
                  do il = 1 , NGL

                      nr = nr_spx(ib)%index(i);
                      ml = ml_spx(ib)%index(i);
                      r  = rb_fam(il);

                      dphi_r(ib)%mat(il,i) = d_phi_nr_ml( nr , abs(ml) , bp , r );

                  end do
              end do
          end do

          ! Calculating wPhi.
          do ib = 1 , N_blocks
              do i = 1 , id_spx(ib)
                  do il = 1 , NGL
                      do ih = 1 , NGH

                          call assert( wzwr(ih,il) >= 0.d0 , 'Negative quadrature weights!' );

                          wPhi(ib)%arr(ih,il,i) = sqrt(wzwr(ih,il)) * phi_z(ib)%mat(ih,i) * phi_r(ib)%mat(il,i);

                      end do
                  end do
              end do
          end do




          ! Calculating U_PHI and V_PHI such that PHI = U_PHI * V_PHI.
          allocate( PHI_tmp( 1:N_total , 1:(2*NGH)*NGL ) );
          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;
                  ihl = ihl + 1;

                  do ib = 1 , N_blocks
                      do i = 1 , id_spx(ib)

                          PHI_tmp( ia_spx(ib)+i-1 , ihl ) = phi_z(ib)%mat(ih,i) * phi_r(ib)%mat(il,i);

                      end do
                  end do

              end do
          end do

          call calculateSVD( A=PHI_tmp , rank=size(U_PHI,2) , U=U_PHI_tmp , V=V_PHI_tmp );

          U_PHI(:,:) = U_PHI_tmp(:,:);
          V_PHI(:,:) = V_PHI_tmp(:,:);
          call assert( getSVDerror(PHI_tmp,U_PHI,V_PHI) < 1.d-10 , 'SVD error is too large!' );

#ifdef DEBUG
          write(*,'(a,e15.6,a)') 'SVD error || PHI - U_PHI*V_PHI ||_F / || PHI ||_F = ' , getSVDerror( PHI_tmp , U_PHI , V_PHI ) , '.';
          write(*,'(a,i0,a,i0,a)') 'Rank of PHI matrix is ' , size(U_PHI,2) , ', compared to full rank of ' , min(size(PHI_tmp,1),size(PHI_tmp,2)) , '.';
#endif
          deallocate(   PHI_tmp );
          deallocate( U_PHI_tmp ); ! Allocated by calculateSVD.
          deallocate( V_PHI_tmp ); ! Allocated by calculateSVD.




      return;
      end subroutine init_basis
