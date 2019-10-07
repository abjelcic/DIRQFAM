c======================================================================c

      subroutine init_coulomb( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /fam/ omega_start, omega_end, delta_omega, omega_print,
     &             omega, gamma_smear,
     &             i_calculation_type, i_coulomb, i_pairing,
     &             J_multipole, K_multipole, ISO;

      common /quadrature/ zb_fam( 1:NGH ), wz( 1:NGH ),
     &                    rb_fam( 1:NGL ), wr( 1:NGL ),
     &                    wzwr( 1:NGH , 1:NGL );

      COMPLEX*16 dV_Cou;
      common /fam_coul/ dV_Cou( -NGH:NGH , 1:NGL ),
     &                       G( -NGH:NGH , 1:NGL , -NGH:NGH , 1:NGL );
                            !G( z' , r' , z , r )



      REAL*8 I_K;
      hbc   = 197.328284D0;
      alpha = 1.D0 / 137.03602D0;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN init_coulomb() *************************';
      write(6,*) '';
      endif


      if( i_coulomb .eq. 0 ) then
          return;
      endif






      G = 0.D0;
      do il1 = 1 , NGL
          r1 = rb_fam(il1);

          do ih1 = -NGH , +NGH
              if( ih1 .eq. 0 ) CYCLE;
              z1 = DBLE(isign(1,ih1)) * zb_fam(abs(ih1));

              do il2 = 1 , NGL
                  r2 = rb_fam(il2);

                  do ih2 = -NGH , +NGH
                      if( ih2 .eq. 0 ) CYCLE;
                      z2 = DBLE(isign(1,ih2)) * zb_fam(abs(ih2));

                      a = 4.D0*r1*r2/( (r1+r2)**2.D0 + (z1-z2)**2.D0 );

                      fac =  wzwr(abs(ih2),il2) * (2.D0*hbc*alpha)
     &                     * DSQRT( (r1+r2)**2.D0 + (z1-z2)**2.D0 );

                      G(ih2,il2,ih1,il1) = fac * I_K( a , K_multipole );

                  enddo

              enddo

          enddo

      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END init_coulomb() ***************************';
      write(6,*) '';
      endif

      return;
      end;
