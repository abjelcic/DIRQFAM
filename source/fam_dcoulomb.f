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

      COMPLEX*16 drho_v, drho_s, ldrho_v, ldrho_s;
      common /ind_dens/ drho_v ( -NGH:NGH , 1:NGL , 2 ),
     &                  drho_s ( -NGH:NGH , 1:NGL     ),
     &                  ldrho_v( -NGH:NGH , 1:NGL , 2 ),
     &                  ldrho_s( -NGH:NGH , 1:NGL     );

      COMPLEX*16 dV_Cou;
      common /fam_coul/ dV_Cou( -NGH:NGH , 1:NGL ),
     &                       G( -NGH:NGH , 1:NGL , -NGH:NGH , 1:NGL );



      COMPLEX*16 z;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_dcoulomb() *************************';
      write(6,*) '';
      endif


      dV_Cou = COMPLEX( 0.D0 , 0.D0 );


      if( i_coulomb .eq. 0 ) then
          return;
      endif

      do il = 1 , NGL
          do ih = -NGH , NGH
              if( ih .eq. 0 ) CYCLE;

              z = COMPLEX( 0.D0 , 0.D0 );
              do il1 = 1 , NGL
                  do ih1 = -NGH , NGH
                      if( ih1 .eq. 0 ) CYCLE;
                      z = z + G(ih1,il1,ih,il) * ldrho_v(ih1,il1,2);
                  enddo
              enddo

              dV_Cou(ih,il) = z;

          enddo
      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_dcoulomb() ***************************';
      write(6,*) '';
      endif

      return;
      end;
