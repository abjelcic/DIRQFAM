!======================================================================!

      subroutine fam_dcoulomb()

!======================================================================!
      use fam_input;
      use fam_green;
      use dlaplace;
      use dcoulomb;
      implicit none;

      integer                    :: il, ih, il1, ih1;
      double complex             :: dV0, dV1, dV2, dV3;
      double complex , parameter :: czero = cmplx( 0.d0 , 0.d0 , kind=8 );


      dVCoulomb_0( -NGH:+NGH , 1:NGL ) = czero;
      dVCoulomb_r( -NGH:+NGH , 1:NGL ) = czero;
      dVCoulomb_p( -NGH:+NGH , 1:NGL ) = czero;
      dVCoulomb_z( -NGH:+NGH , 1:NGL ) = czero;

      if( include_coulomb == 0 ) then
          return;
      end if

      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih == 0 ) cycle;

              dV0 = czero;
              dV1 = czero;
              dV2 = czero;
              dV3 = czero;
              do il1 = 1 , NGL
                  do ih1 = -NGH , +NGH
                      if( ih1 == 0 ) cycle;

                      dV0 = dV0 + G3(ih1,il1,ih,il) * ldrho_vprot(ih1,il1);
                      dV1 = dV1 + G1(ih1,il1,ih,il) *   ldj_1prot(ih1,il1);
                      dV2 = dV2 + G2(ih1,il1,ih,il) *   ldj_2prot(ih1,il1);
                      dV3 = dV3 + G3(ih1,il1,ih,il) *   ldj_3prot(ih1,il1);

                  end do
              end do

              dVCoulomb_0(ih,il) = dV0;
              dVCoulomb_r(ih,il) = ( + dV1 + dV2 )/sqrt(2.d0);
              dVCoulomb_p(ih,il) = ( - dV1 + dV2 )/sqrt(2.d0);
              dVCoulomb_z(ih,il) = dV3;

          end do
      end do


      return;
      end subroutine fam_dcoulomb
