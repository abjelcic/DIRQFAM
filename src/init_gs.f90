!======================================================================!

      subroutine init_gs()

!======================================================================!
      use fields;
      use dataTypes;
      use fam_input;
      use simplex;
      use v_matrix;
      use ddpc1ddme2;
      use basis;
      use gs_dens;
      use gs_mesons;
      implicit none;

      integer                                 :: it, ib, il, ih, ihl, i, j;
      character                               :: fg1, fg2;
      integer                                 :: ml1, ml2;
      type(complexBlockMatrix) , dimension(2) :: rho1;
      double complex           , parameter    :: cone  = cmplx( 1.d0 , 0.d0 , kind=8 );
      double complex           , parameter    :: czero = cmplx( 0.d0 , 0.d0 , kind=8 );
      double precision         , parameter    :: pi    = 4.d0*atan(1.d0);
      double precision         , parameter    :: hbc   = 197.328284d0;




      ! Density matrix rho = [ rho1 , 0 ; 0 , rho1^T ] = [ v*v' , 0 ; 0 , (v*v')^T ].
      do it = 1 , 2
          allocate( rho1(it)%blocks(1:N_blocks,1:N_blocks) );
          do ib = 1 , N_blocks
              allocate( rho1(it)%blocks(ib,ib)%mat(1:id_spx(ib),1:id_spx(ib)) );

              call zgemm( 'n' , 'c' , id_spx(ib) , id_spx(ib) , id_spx(ib)                    , &
                          cone                                                                , &
                          v(it)%blocks(ib,ib)%mat(1,1)    , size(v(it)%blocks(ib,ib)%mat,1)   , &
                          v(it)%blocks(ib,ib)%mat(1,1)    , size(v(it)%blocks(ib,ib)%mat,1)   , &
                          czero                                                               , &
                          rho1(it)%blocks(ib,ib)%mat(1,1) , size(rho1(it)%blocks(ib,ib)%mat,1)  );
          end do
      end do




      ! Calculation of ground state isoscalar-scalar and isoscalar-vector densities.
      rhov_GS( -NGH:+NGH , 1:NGL , 1:2 ) = 0.d0;
      rhos_GS( -NGH:+NGH , 1:NGL       ) = 0.d0;
      do it = 1 , 2
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;

                  do ib = 1 , N_blocks
                      do j = 1 , id_spx(ib)
                          do i = 1 , id_spx(ib)

                              fg1 = fg_spx(ib)%index(i);
                              ml1 = ml_spx(ib)%index(i);

                              fg2 = fg_spx(ib)%index(j);
                              ml2 = ml_spx(ib)%index(j);

                              if( all([fg1,fg2]==['f','f']) .and. abs(ml1-ml2)==0 ) then
                                  rhov_GS(ih,il,it) = rhov_GS(ih,il,it) + 2.d0*real(rho1(it)%blocks(ib,ib)%mat(i,j)) * phi_z(ib)%mat(ih,i)*phi_z(ib)%mat(ih,j) * phi_r(ib)%mat(il,i)*phi_r(ib)%mat(il,j) / (2.d0*pi);
                                  rhos_GS(ih,il   ) = rhos_GS(ih,il   ) + 2.d0*real(rho1(it)%blocks(ib,ib)%mat(i,j)) * phi_z(ib)%mat(ih,i)*phi_z(ib)%mat(ih,j) * phi_r(ib)%mat(il,i)*phi_r(ib)%mat(il,j) / (2.d0*pi);
                              end if
                              if( all([fg1,fg2]==['g','g']) .and. abs(ml1-ml2)==0 ) then
                                  rhov_GS(ih,il,it) = rhov_GS(ih,il,it) + 2.d0*real(rho1(it)%blocks(ib,ib)%mat(i,j)) * phi_z(ib)%mat(ih,i)*phi_z(ib)%mat(ih,j) * phi_r(ib)%mat(il,i)*phi_r(ib)%mat(il,j) / (2.d0*pi);
                                  rhos_GS(ih,il   ) = rhos_GS(ih,il   ) - 2.d0*real(rho1(it)%blocks(ib,ib)%mat(i,j)) * phi_z(ib)%mat(ih,i)*phi_z(ib)%mat(ih,j) * phi_r(ib)%mat(il,i)*phi_r(ib)%mat(il,j) / (2.d0*pi);
                              end if

                          end do
                      end do
                  end do

              end do
          end do
      end do




      ! Setting the ground state meson fields.
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih == 0 ) cycle;

              ihl = 1+abs(ih) + il*(NGH+1);

              sigma_GS(ih,il) = - g0_sigma * phi(ihl,1) / (m_sigma/hbc)**2;
              omega_GS(ih,il) = + g0_omega * phi(ihl,2) / (m_omega/hbc)**2;
                rho_GS(ih,il) = + g0_rho   * phi(ihl,4) / (m_rho  /hbc)**2;

          end do
      end do




      ! Deallocate rho1.
      do it = 1 , 2
          do ib = 1 , N_blocks
              deallocate( rho1(it)%blocks(ib,ib)%mat );
          end do
          deallocate( rho1(it)%blocks );
      end do




      return;
      end subroutine init_gs
