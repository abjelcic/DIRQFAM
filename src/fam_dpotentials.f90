!======================================================================!

      subroutine fam_dpotentials()

!======================================================================!
      use fam_input;
      use ddpc1ddme2;
      use gs_dens;
      use gs_mesons;
      use ddens;
      use dcurr;
      use dlaplace;
      use dmesons;
      use dcoulomb;
      use dpotentials;
      implicit none;

      double precision , parameter :: hbc = 197.328284d0;
      integer                      :: it, ih, il;
      double precision             :: rho0_s;
      double precision             :: rho0_v;
      double precision             :: rho0_tv;
      double precision , external  :: tau3;
      double complex               :: dSig_s ( -NGH:+NGH , 1:NGL       );
      double complex               :: dSig0  ( -NGH:+NGH , 1:NGL , 1:2 );
      double complex               :: dSig0_R( -NGH:+NGH , 1:NGL       );
      double complex               :: dSig_z ( -NGH:+NGH , 1:NGL , 1:2 );
      double complex               :: dSig_r ( -NGH:+NGH , 1:NGL , 1:2 );
      double complex               :: dSig_p ( -NGH:+NGH , 1:NGL , 1:2 );




      call assert( LagrangianModelName=='DD-PC1' .or. LagrangianModelName=='DD-ME2' , 'Only DD-PC1 and DD-ME2 supported.' );

      dSig_s (:,:  ) = cmplx( 0.d0 , 0.d0 , kind=8 );
      dSig0  (:,:,:) = cmplx( 0.d0 , 0.d0 , kind=8 );
      dSig0_R(:,:  ) = cmplx( 0.d0 , 0.d0 , kind=8 );
      dSig_z (:,:,:) = cmplx( 0.d0 , 0.d0 , kind=8 );
      dSig_r (:,:,:) = cmplx( 0.d0 , 0.d0 , kind=8 );
      dSig_p (:,:,:) = cmplx( 0.d0 , 0.d0 , kind=8 );




      if( LagrangianModelName == 'DD-PC1' ) then

          ! Calculation of dSig_s, see Eq. (41).
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;

                  rho0_s =        rhos_GS(ih,il);
                  rho0_v = sum([( rhov_GS(ih,il,it) , it=1,2 )]);

                  dSig_s(ih,il) = + ( dalpha_s(rho0_v)*rho0_s ) * sum([( drho_v(ih,il,it) , it=1,2 )]) &
                                  + (  alpha_s(rho0_v)        ) *        drho_s(ih,il)                 &
                                  + (  delta_s                ) *       ldrho_s(ih,il)                 ;
              end do
          end do

          ! Calculation of dSig0, see Eq. (42).
          do it = 1 , 2
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;

                      rho0_s  =                 rhos_GS(ih,il);
                      rho0_v  = sum([(          rhov_GS(ih,il,it) , it=1,2 )]);
                      rho0_tv = sum([( tau3(it)*rhov_GS(ih,il,it) , it=1,2 )]);

                      dSig0(ih,il,it) = + ( dalpha_v(rho0_v)*rho0_v + alpha_v(rho0_v) + tau3(it)*dalpha_tv(rho0_v)*rho0_tv ) * sum([(          drho_v(ih,il,it) , it=1,2 )]) &
                                        + ( tau3(it)*alpha_tv(rho0_v)                                                      ) * sum([( tau3(it)*drho_v(ih,il,it) , it=1,2 )]) ;
                  end do
              end do
          end do

          ! Calculation of dSig0_R, see Eq. (43).
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;

                  rho0_s  =                 rhos_GS(ih,il);
                  rho0_v  = sum([(          rhov_GS(ih,il,it) , it=1,2 )]);
                  rho0_tv = sum([( tau3(it)*rhov_GS(ih,il,it) , it=1,2 )]);

                  dSig0_R(ih,il) = + 0.5d0 * ( ddalpha_s(rho0_v)*rho0_s**2 + ddalpha_v(rho0_v)*rho0_v**2 + ddalpha_tv(rho0_v)*rho0_tv**2 ) * sum([(          drho_v(ih,il,it) , it=1,2 )]) &
                                   +         ( dalpha_s(rho0_v)*rho0_s                                                                   ) *                 drho_s(ih,il)                 &
                                   +         ( dalpha_v(rho0_v)*rho0_v                                                                   ) * sum([(          drho_v(ih,il,it) , it=1,2 )]) &
                                   +         ( dalpha_tv(rho0_v)*rho0_tv                                                                 ) * sum([( tau3(it)*drho_v(ih,il,it) , it=1,2 )]) ;
              end do
          end do

          ! Calculation of dSig_z, see Eq. (44).
          do it = 1 , 2
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;

                      rho0_v = sum([( rhov_GS(ih,il,it) , it=1,2 )]);

                      dSig_z(ih,il,it) = + ( alpha_v(rho0_v)           ) * sum([(          dj_z(ih,il,it) , it=1,2 )]) &
                                         + ( tau3(it)*alpha_tv(rho0_v) ) * sum([( tau3(it)*dj_z(ih,il,it) , it=1,2 )]) ;
                  end do
              end do
          end do

          ! Calculation of dSig_r, see Eq. (44).
          do it = 1 , 2
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;

                      rho0_v = sum([( rhov_GS(ih,il,it) , it=1,2 )]);

                      dSig_r(ih,il,it) = + ( alpha_v(rho0_v)           ) * sum([(          dj_r(ih,il,it) , it=1,2 )]) &
                                         + ( tau3(it)*alpha_tv(rho0_v) ) * sum([( tau3(it)*dj_r(ih,il,it) , it=1,2 )]) ;
                  end do
              end do
          end do

          ! Calculation of dSig_p, see Eq. (44).
          do it = 1 , 2
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;

                      rho0_v = sum([( rhov_GS(ih,il,it) , it=1,2 )]);

                      dSig_p(ih,il,it) = + ( alpha_v(rho0_v)           ) * sum([(          dj_p(ih,il,it) , it=1,2 )]) &
                                         + ( tau3(it)*alpha_tv(rho0_v) ) * sum([( tau3(it)*dj_p(ih,il,it) , it=1,2 )]) ;
                  end do
              end do
          end do

      end if




      if( LagrangianModelName == 'DD-ME2' ) then
          call fam_dmesons();

          ! Calculation of dSig_s.
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;

                  rho0_v = sum([( rhov_GS(ih,il,it) , it=1,2 )]);

                  dSig_s(ih,il) = + (  g_sigma(rho0_v)                   ) *        dsigma(ih,il)                 &
                                  + ( dg_sigma(rho0_v) * sigma_GS(ih,il) ) * sum([( drho_v(ih,il,it) , it=1,2 )]) ;
              end do
          end do

          ! Calculation of dSig0.
          do it = 1 , 2
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;

                      rho0_v = sum([( rhov_GS(ih,il,it) , it=1,2 )]);

                      dSig0(ih,il,it) = + (             g_omega(rho0_v)                   ) *        domega_0(ih,il)               &
                                        + (            dg_omega(rho0_v) * omega_GS(ih,il) ) * sum([( drho_v(ih,il,it) , it=1,2 )]) &
                                        + ( tau3(it) *  g_rho(rho0_v)                     ) *        drho_0(ih,il)                 &
                                        + ( tau3(it) * dg_rho(rho0_v)   * rho_GS(ih,il)   ) * sum([( drho_v(ih,il,it) , it=1,2 )]) ;
                  end do
              end do
          end do

          ! Calculation of dSig0_R.
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;

                  rho0_s  =                 rhos_GS(ih,il);
                  rho0_v  = sum([(          rhov_GS(ih,il,it) , it=1,2 )]);
                  rho0_tv = sum([( tau3(it)*rhov_GS(ih,il,it) , it=1,2 )]);

                  dSig0_R(ih,il) = + ( + ddg_sigma(rho0_v) * rho0_s  * sigma_GS(ih,il)                                                     &
                                       + ddg_omega(rho0_v) * rho0_v  * omega_GS(ih,il)                                                     &
                                       +   ddg_rho(rho0_v) * rho0_tv *   rho_GS(ih,il)   ) * sum([(          drho_v(ih,il,it) , it=1,2 )]) &
                                   + ( dg_sigma(rho0_v) * sigma_GS(ih,il)                ) *                 drho_s(ih,il)                 &
                                   + ( dg_omega(rho0_v) * omega_GS(ih,il)                ) * sum([(          drho_v(ih,il,it) , it=1,2 )]) &
                                   + (   dg_rho(rho0_v) *   rho_GS(ih,il)                ) * sum([( tau3(it)*drho_v(ih,il,it) , it=1,2 )]) &
                                   + ( dg_sigma(rho0_v) * rho0_s                         ) *                 dsigma(ih,il)                 &
                                   + ( dg_omega(rho0_v) * rho0_v                         ) *               domega_0(ih,il)                 &
                                   + (   dg_rho(rho0_v) * rho0_tv                        ) *                 drho_0(ih,il)                 ;
              end do
          end do

          ! Calculation of dSig_z.
          do it = 1 , 2
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;

                      rho0_v = sum([( rhov_GS(ih,il,it) , it=1,2 )]);

                      dSig_z(ih,il,it) = + ( g_omega(rho0_v)        ) * domega_z(ih,il) &
                                         + ( tau3(it)*g_rho(rho0_v) ) *   drho_z(ih,il) ;
                  end do
              end do
          end do

          ! Calculation of dSig_r.
          do it = 1 , 2
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;

                      rho0_v = sum([( rhov_GS(ih,il,it) , it=1,2 )]);

                      dSig_r(ih,il,it) = + ( g_omega(rho0_v)        ) * domega_r(ih,il) &
                                         + ( tau3(it)*g_rho(rho0_v) ) *   drho_r(ih,il) ;
                  end do
              end do
          end do

          ! Calculation of dSig_p.
          do it = 1 , 2
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;

                      rho0_v = sum([( rhov_GS(ih,il,it) , it=1,2 )]);

                      dSig_p(ih,il,it) = + ( g_omega(rho0_v)        ) * domega_p(ih,il) &
                                         + ( tau3(it)*g_rho(rho0_v) ) *   drho_p(ih,il) ;
                  end do
              end do
          end do

      end if




      ! Calculation of the induced potentials.
      do it = 1 , 2
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;

                  dVplusS (ih,il,it) = hbc * ( ( dSig0(ih,il,it) + dSig0_R(ih,il) ) + dSig_s(ih,il) );
                  dVminusS(ih,il,it) = hbc * ( ( dSig0(ih,il,it) + dSig0_R(ih,il) ) - dSig_s(ih,il) );
                  dSigma_z(ih,il,it) = hbc * dSig_z(ih,il,it);
                  dSigma_r(ih,il,it) = hbc * dSig_r(ih,il,it);
                  dSigma_p(ih,il,it) = hbc * dSig_p(ih,il,it);

              end do
          end do
      end do

      ! Induced Coulomb interaction.
      call fam_dcoulomb();
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih == 0 ) cycle;

              dVplusS (ih,il,2) = dVplusS (ih,il,2) + dVCoulomb_0(ih,il);
              dVminusS(ih,il,2) = dVminusS(ih,il,2) + dVCoulomb_0(ih,il);
              dSigma_z(ih,il,2) = dSigma_z(ih,il,2) + dVCoulomb_z(ih,il);
              dSigma_r(ih,il,2) = dSigma_r(ih,il,2) + dVCoulomb_r(ih,il);
              dSigma_p(ih,il,2) = dSigma_p(ih,il,2) + dVCoulomb_p(ih,il);

          end do
      end do




      return;
      end subroutine fam_dpotentials
