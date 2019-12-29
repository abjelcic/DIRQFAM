c======================================================================c

      subroutine fam_dcoulomb( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /fam/ omega_start, omega_end, delta_omega, omega_print,
     &             omega, gamma_smear,
     &             i_calculation_type, i_coulomb, i_pairing,
     &             J_multipole, K_multipole, ISO;

      common /fam_green/ G1( -NGH:NGH , 1:NGL , -NGH:NGH , 1:NGL ),
     &                   G2( -NGH:NGH , 1:NGL , -NGH:NGH , 1:NGL ),
     &                   G3( -NGH:NGH , 1:NGL , -NGH:NGH , 1:NGL );
                        !G (     z'   ,  r'   ,     z    ,  r    )

      COMPLEX*16 ldrho_vp, ldrho_s, ldj_1p, ldj_2p, ldj_3p;
      common /laplace/ ldrho_vp( -NGH:NGH , 1:NGL ),
     &                 ldrho_s ( -NGH:NGH , 1:NGL ),
     &                 ldj_1p  ( -NGH:NGH , 1:NGL ),
     &                 ldj_2p  ( -NGH:NGH , 1:NGL ),
     &                 ldj_3p  ( -NGH:NGH , 1:NGL );

      COMPLEX*16 dVCou_0, dVCou_r, dVCou_p, dVCou_z;
      common /fam_coul/ dVCou_0( -NGH:NGH , 1:NGL ),
     &                  dVCou_r( -NGH:NGH , 1:NGL ),
     &                  dVCou_p( -NGH:NGH , 1:NGL ),
     &                  dVCou_z( -NGH:NGH , 1:NGL );



      COMPLEX*16 acc0, acc1, acc2, acc3;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_dcoulomb() ************************';
      write(6,*) '';
      endif


      dVCou_0 = COMPLEX( 0.D0 , 0.D0 );
      dVCou_r = COMPLEX( 0.D0 , 0.D0 );
      dVCou_p = COMPLEX( 0.D0 , 0.D0 );
      dVCou_z = COMPLEX( 0.D0 , 0.D0 );


      if( i_coulomb .eq. 0 ) then
          return;
      endif


      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              acc0 = COMPLEX( 0.D0 , 0.D0 );
              acc1 = COMPLEX( 0.D0 , 0.D0 );
              acc2 = COMPLEX( 0.D0 , 0.D0 );
              acc3 = COMPLEX( 0.D0 , 0.D0 );
              do il1 = 1 , NGL
                  do ih1 = -NGH , +NGH
                      if( ih1 .eq. 0 ) CYCLE;
                      acc0 = acc0 + G3(ih1,il1,ih,il)*ldrho_vp(ih1,il1);
                      acc1 = acc1 + G1(ih1,il1,ih,il)*  ldj_1p(ih1,il1);
                      acc2 = acc2 + G2(ih1,il1,ih,il)*  ldj_2p(ih1,il1);
                      acc3 = acc3 + G3(ih1,il1,ih,il)*  ldj_3p(ih1,il1);
                  enddo
              enddo

              dVCou_0(ih,il) = acc0;
              dVCou_r(ih,il) = ( + acc1 + acc2 ) / DSQRT(2.D0);
              dVCou_p(ih,il) = ( - acc1 + acc2 ) / DSQRT(2.D0);
              dVCou_z(ih,il) = acc3;

          enddo
      enddo



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_dcoulomb() **************************';
      write(6,*) '';
      endif

      return;
      end;
