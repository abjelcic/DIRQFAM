!======================================================================!

      subroutine fam_dmesons()

!======================================================================!
      use fam_input;
      use ddpc1ddme2;
      use gs_dens;
      use mesmat;
      use ddens;
      use dcurr;
      use dmesons;
      implicit none;

      integer                     :: it, ih, il;
      double precision            :: rho0_s;
      double precision            :: rho0_v;
      double precision            :: rho0_tv;
      double precision , external :: tau3;
      double complex              :: dPhi0( -NGH:+NGH , 1:NGL );
      double complex              :: dPhi1( -NGH:+NGH , 1:NGL );
      double complex              :: dPhi2( -NGH:+NGH , 1:NGL );
      double complex              :: dPhi3( -NGH:+NGH , 1:NGL );
      double complex              :: dS0  ( -NGH:+NGH , 1:NGL );
      double complex              :: dS1  ( -NGH:+NGH , 1:NGL );
      double complex              :: dS2  ( -NGH:+NGH , 1:NGL );
      double complex              :: dS3  ( -NGH:+NGH , 1:NGL );




      ! Calculating dsigma.
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih == 0 ) cycle;

              rho0_s =        rhos_GS(ih,il);
              rho0_v = sum([( rhov_GS(ih,il,it) , it=1,2 )]);

              dS0(ih,il) = - (  g_sigma(rho0_v)          ) *        drho_s(ih,il)                 &
                           - ( dg_sigma(rho0_v) * rho0_s ) * sum([( drho_v(ih,il,it) , it=1,2 )]) ;

          end do
      end do
      call solveKleinGordon( dS0 , Psig(abs(K_multipole)) , dPhi0 );
      dsigma(:,:) = dPhi0(:,:);




      ! Calculating domega_0.
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih == 0 ) cycle;

              rho0_v = sum([( rhov_GS(ih,il,it) , it=1,2 )]);

              dS0(ih,il) = + (  g_omega(rho0_v)          ) * sum([( drho_v(ih,il,it) , it=1,2 )]) &
                           + ( dg_omega(rho0_v) * rho0_v ) * sum([( drho_v(ih,il,it) , it=1,2 )]) ;

          end do
      end do
      call solveKleinGordon( dS0 , Pome(abs(K_multipole)) , dPhi0 );
      domega_0(:,:) = dPhi0(:,:);




      ! Calculating drho_0.
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih == 0 ) cycle;

              rho0_v  = sum([(            rhov_GS(ih,il,it) , it=1,2 )]);
              rho0_tv = sum([( tau3(it) * rhov_GS(ih,il,it) , it=1,2 )]);


              dS0(ih,il) = + (  g_rho(rho0_v)           ) * sum([( tau3(it) * drho_v(ih,il,it) , it=1,2 )]) &
                           + ( dg_rho(rho0_v) * rho0_tv ) * sum([(            drho_v(ih,il,it) , it=1,2 )]) ;

          end do
      end do
      call solveKleinGordon( dS0 , Prho(abs(K_multipole)) , dPhi0 );
      drho_0(:,:) = dPhi0(:,:);




      ! Calculating domega_r, domega_p and domega_z.
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih == 0 ) cycle;

              rho0_v = sum([( rhov_GS(ih,il,it) , it=1,2 )]);

              dS1(ih,il) = + ( g_omega(rho0_v) ) * sum([( dj_1(ih,il,it) , it=1,2 )]);
              dS2(ih,il) = + ( g_omega(rho0_v) ) * sum([( dj_2(ih,il,it) , it=1,2 )]);
              dS3(ih,il) = + ( g_omega(rho0_v) ) * sum([( dj_3(ih,il,it) , it=1,2 )]);

          end do
      end do
      call solveKleinGordon( dS1 , Pome(abs(K_multipole-1)) , dPhi1 );
      call solveKleinGordon( dS2 , Pome(abs(K_multipole+1)) , dPhi2 );
      call solveKleinGordon( dS3 , Pome(abs(K_multipole+0)) , dPhi3 );
      domega_r(:,:) = ( + dPhi1(:,:) + dPhi2(:,:) ) / sqrt(2.d0);
      domega_p(:,:) = ( - dPhi1(:,:) + dPhi2(:,:) ) / sqrt(2.d0);
      domega_z(:,:) = ( + dPhi3(:,:)              );




      ! Calculating drho_r, drho_p and drho_z.
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih == 0 ) cycle;

              rho0_v = sum([( rhov_GS(ih,il,it) , it=1,2 )]);

              dS1(ih,il) = + ( g_rho(rho0_v) ) * sum([( tau3(it) * dj_1(ih,il,it) , it=1,2 )]);
              dS2(ih,il) = + ( g_rho(rho0_v) ) * sum([( tau3(it) * dj_2(ih,il,it) , it=1,2 )]);
              dS3(ih,il) = + ( g_rho(rho0_v) ) * sum([( tau3(it) * dj_3(ih,il,it) , it=1,2 )]);

          end do
      end do
      call solveKleinGordon( dS1 , Prho(abs(K_multipole-1)) , dPhi1 );
      call solveKleinGordon( dS2 , Prho(abs(K_multipole+1)) , dPhi2 );
      call solveKleinGordon( dS3 , Prho(abs(K_multipole+0)) , dPhi3 );
      drho_r(:,:) = ( + dPhi1(:,:) + dPhi2(:,:) ) / sqrt(2.d0);
      drho_p(:,:) = ( - dPhi1(:,:) + dPhi2(:,:) ) / sqrt(2.d0);
      drho_z(:,:) = ( + dPhi3(:,:)              );




      return;
      end subroutine fam_dmesons






!======================================================================!

      subroutine solveKleinGordon( dS , P , dPhi );

!======================================================================!
      use dataTypes;
      use fam_input;
      use quadrature;
      implicit none;

      double complex    , dimension( -NGH:+NGH , 1:NGL ) , intent(in)    :: dS;
      type(real2Darray) ,                                  intent(in)    :: P;
      double complex    , dimension( -NGH:+NGH , 1:NGL ) , intent(inout) :: dPhi;

      integer                      :: ihl, ih, il;
      double precision , parameter :: one  = 1.d0;
      double precision , parameter :: zero = 0.d0;
      double precision             ::    wdS( 1:size(P%mat,2) , 1:2 );
      double precision             ::   PwdS( 1:size(P%mat,1) , 1:2 );
      double precision             :: PtPwdS( 1:size(P%mat,2) , 1:2 );


      ihl = 0;
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih == 0 ) cycle;
              ihl = ihl + 1;

              wdS(ihl,1) = wzwr(abs(ih),il) * real( dS(ih,il) );
              wdS(ihl,2) = wzwr(abs(ih),il) * imag( dS(ih,il) );

          end do
      end do
      call assert( ihl == size(P%mat,2) , 'P size error in solveKleinGordon.' );

      ! PwdS(:,1:2) = P%mat * wdS(:,1:2).
      call dgemm( 'n' , 'n' , size(P%mat,1) , 2 , size(P%mat,2) , &
                  one                                           , &
                  P%mat(1,1)  , size(P%mat,1)                   , &
                  wdS(1,1)    , size(wdS,1)                     , &
                  zero                                          , &
                  PwdS(1,1)   , size(PwdS,1)                      );

      ! PtPwdS(:,1:2) = [(P%mat)^T*(P%mat)] * wdS(:,1:2).
      call dgemm( 't' , 'n' , size(P%mat,2) , 2 , size(P%mat,1) , &
                  one                                           , &
                  P%mat(1,1)  , size(P%mat,1)                   , &
                  PwdS(1,1)   , size(PwdS,1)                    , &
                  zero                                          , &
                  PtPwdS(1,1) , size(PtPwdS,1)                    );

      ihl = 0;
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih == 0 ) cycle;
              ihl = ihl + 1;

              dPhi(ih,il) = cmplx( PtPwdS(ihl,1) , PtPwdS(ihl,2) , kind=8 );

          end do
      end do
      call assert( ihl == size(PtPwdS,1) , 'PtPwdS size error in solveKleinGordon.' );


      return;
      end subroutine solveKleinGordon
