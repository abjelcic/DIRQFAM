c======================================================================c

      subroutine fam_dmesons( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /fam/ omega_start, omega_end, delta_omega, omega_print,
     &             omega, gamma_smear,
     &             i_calculation_type, i_coulomb, i_pairing,
     &             J_multipole, K_multipole, ISO;

      REAL*8 m_sig, m_ome, m_rho;
      common /DDPC1_DDME2/  a_s , b_s , c_s , d_s ,
     &                      a_v , b_v , c_v , d_v ,
     &                      a_tv, b_tv, c_tv, d_tv,
     &                      del_s,
     &
     &                      a_sig, b_sig, c_sig, d_sig, g0_sig, m_sig,
     &                      a_ome, b_ome, c_ome, d_ome, g0_ome, m_ome,
     &                      a_rho,                      g0_rho, m_rho,
     &
     &                      rho_sat;

      common /gs_dens/ rhov_GS( -NGH:NGH , 1:NGL , 2 ),
     &                 rhos_GS( -NGH:NGH , 1:NGL     );

      common /mesmat/ Psig( NHMAX , NCOORD , 0:J_MAX+1 ),
     &                Pome( NHMAX , NCOORD , 0:J_MAX+1 ),
     &                Prho( NHMAX , NCOORD , 0:J_MAX+1 ),
     &                nP( 0:J_MAX+1 );

      COMPLEX*16 drho_v, drho_s;
      common /ind_dens/ drho_v( -NGH:NGH , 1:NGL , 2 ),
     &                  drho_s( -NGH:NGH , 1:NGL     );

      COMPLEX*16 dj_r, dj_p, dj_z, dj_1, dj_2, dj_3;
      common /ind_curr/ dj_r( -NGH:NGH , 1:NGL , 2 ),
     &                  dj_p( -NGH:NGH , 1:NGL , 2 ),
     &                  dj_z( -NGH:NGH , 1:NGL , 2 ),
     &                  dj_1( -NGH:NGH , 1:NGL , 2 ),
     &                  dj_2( -NGH:NGH , 1:NGL , 2 ),
     &                  dj_3( -NGH:NGH , 1:NGL , 2 );

      COMPLEX*16 dsig  ,
     &           dome_0, dome_r, dome_p, dome_z,
     &           drho_0, drho_r, drho_p, drho_z;
      common /ind_mesons/ dsig  ( -NGH:NGH , 1:NGL ),
     &                    dome_0( -NGH:NGH , 1:NGL ),
     &                    dome_r( -NGH:NGH , 1:NGL ),
     &                    dome_p( -NGH:NGH , 1:NGL ),
     &                    dome_z( -NGH:NGH , 1:NGL ),
     &                    drho_0( -NGH:NGH , 1:NGL ),
     &                    drho_r( -NGH:NGH , 1:NGL ),
     &                    drho_p( -NGH:NGH , 1:NGL ),
     &                    drho_z( -NGH:NGH , 1:NGL );



      COMPLEX*16 s ( -NGH:NGH , 1:NGL );
      COMPLEX*16 s1( -NGH:NGH , 1:NGL );
      COMPLEX*16 s2( -NGH:NGH , 1:NGL );
      COMPLEX*16 s3( -NGH:NGH , 1:NGL );
      COMPLEX*16 u1( -NGH:NGH , 1:NGL );
      COMPLEX*16 u2( -NGH:NGH , 1:NGL );
      COMPLEX*16 u3( -NGH:NGH , 1:NGL );

      f0(x,a,b,c,d) =   a       * (1+b*(x+d)**2)   / (1+c*(x+d)**2);
      f1(x,a,b,c,d) = 2*a*(b-c) * (x+d)            / (1+c*(x+d)**2)**2;
      f2(x,a,b,c,d) = 2*a*(b-c) * (1-3*c*(x+d)**2) / (1+c*(x+d)**2)**3;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_dmesons() *************************';
      write(6,*) '';
      endif






c-----Calculating dsig
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              rs = rhos_GS(ih,il);
              rv = rhov_GS(ih,il,1) + rhov_GS(ih,il,2);
              x  = rv / rho_sat;

              g  = g0_sig * f0(x,a_sig,b_sig,c_sig,d_sig);
              dg = g0_sig * f1(x,a_sig,b_sig,c_sig,d_sig) / rho_sat;

              s(ih,il) = - g     * (+drho_s(ih,il)                  )
     &                   - dg*rs * (+drho_v(ih,il,1)+drho_v(ih,il,2));

          enddo
      enddo
      K = K_multipole;
      call solve_kg( s , dsig , nP(K) , Psig(1,1,K) , NHMAX );






c-----Calculating dome_0
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              rv = rhov_GS(ih,il,1) + rhov_GS(ih,il,2);
              x  = rv / rho_sat;

              g  = g0_ome * f0(x,a_ome,b_ome,c_ome,d_ome);
              dg = g0_ome * f1(x,a_ome,b_ome,c_ome,d_ome) / rho_sat;

              s(ih,il) = + g     * (+drho_v(ih,il,1)+drho_v(ih,il,2))
     &                   + dg*rv * (+drho_v(ih,il,1)+drho_v(ih,il,2));

          enddo
      enddo
      K = K_multipole;
      call solve_kg( s , dome_0 , nP(K) , Pome(1,1,K) , NHMAX );






c-----Calculating drho_0
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              rv  = + rhov_GS(ih,il,1) + rhov_GS(ih,il,2);
              rtv = - rhov_GS(ih,il,1) + rhov_GS(ih,il,2);
              x   = rv / rho_sat;

              g   = g0_rho *            DEXP(-a_rho*(x-1.D0));
              dg  = g0_rho * (-a_rho) * DEXP(-a_rho*(x-1.D0)) /rho_sat;

              s(ih,il) = + g      * (-drho_v(ih,il,1)+drho_v(ih,il,2))
     &                   + dg*rtv * (+drho_v(ih,il,1)+drho_v(ih,il,2));

          enddo
      enddo
      K = K_multipole;
      call solve_kg( s , drho_0 , nP(K) , Prho(1,1,K) , NHMAX );






c-----Calculating (dome_r, dome_p, dome_z)
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              rv = rhov_GS(ih,il,1) + rhov_GS(ih,il,2);
              x  = rv / rho_sat;

              g  = g0_ome * f0(x,a_ome,b_ome,c_ome,d_ome);

              s1(ih,il) = + g * (+dj_1(ih,il,1)+dj_1(ih,il,2));
              s2(ih,il) = + g * (+dj_2(ih,il,1)+dj_2(ih,il,2));
              s3(ih,il) = + g * (+dj_3(ih,il,1)+dj_3(ih,il,2));

          enddo
      enddo
      K = iabs( K_multipole - 1 );
      call solve_kg( s1 , u1 , nP(K) , Pome(1,1,K) , NHMAX );
      K = iabs( K_multipole + 1 );
      call solve_kg( s2 , u2 , nP(K) , Pome(1,1,K) , NHMAX );
      K = iabs( K_multipole + 0 );
      call solve_kg( s3 , u3 , nP(K) , Pome(1,1,K) , NHMAX );
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              dome_r(ih,il) = ( + u1(ih,il) + u2(ih,il) ) / DSQRT(2.D0);
              dome_p(ih,il) = ( - u1(ih,il) + u2(ih,il) ) / DSQRT(2.D0);
              dome_z(ih,il) =   + u3(ih,il);

          enddo
      enddo






c-----Calculating (drho_r, drho_p, drho_z)
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              rv = rhov_GS(ih,il,1) + rhov_GS(ih,il,2);
              x  = rv / rho_sat;

              g  = g0_rho * DEXP(-a_rho*(x-1.D0));

              s1(ih,il) = + g * (-dj_1(ih,il,1)+dj_1(ih,il,2));
              s2(ih,il) = + g * (-dj_2(ih,il,1)+dj_2(ih,il,2));
              s3(ih,il) = + g * (-dj_3(ih,il,1)+dj_3(ih,il,2));

          enddo
      enddo
      K = iabs( K_multipole - 1 );
      call solve_kg( s1 , u1 , nP(K) , Prho(1,1,K) , NHMAX );
      K = iabs( K_multipole + 1 );
      call solve_kg( s2 , u2 , nP(K) , Prho(1,1,K) , NHMAX );
      K = iabs( K_multipole + 0 );
      call solve_kg( s3 , u3 , nP(K) , Prho(1,1,K) , NHMAX );
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              drho_r(ih,il) = ( + u1(ih,il) + u2(ih,il) ) / DSQRT(2.D0);
              drho_p(ih,il) = ( - u1(ih,il) + u2(ih,il) ) / DSQRT(2.D0);
              drho_z(ih,il) =   + u3(ih,il);

          enddo
      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_dmesons() ***************************';
      write(6,*) '';
      endif

      return;
      end;






c======================================================================c

      subroutine solve_kg( src , u   ,
     &                     n   , P   , LDP );

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'

      common /quadrature/ zb_fam( 1:NGH ), wz( 1:NGH ),
     &                    rb_fam( 1:NGL ), wr( 1:NGL ),
     &                    wzwr( 1:NGH , 1:NGL );

      COMPLEX*16 src( -NGH:NGH , 1:NGL );
      COMPLEX*16 u  ( -NGH:NGH , 1:NGL );
      INTEGER*4  n;
      INTEGER*4  LDP;
      REAL*8     P( LDP , NCOORD );
      REAL*8    ws(       NCOORD );
      REAL*8   Pws( LDP          );
      REAL*8 PtPws(       NCOORD );



      u = COMPLEX( 0.D0 , 0.D0 );



c-----Real part
      ihl = 0;
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;
              ihl = ihl + 1;
              ws(ihl) = wzwr(iabs(ih),il) * DREAL( src(ih,il) );
          enddo
      enddo
      call dgemv(  'N' ,
     &               n , NCOORD ,
     &            1.D0 ,
     &               P ,    LDP ,
     &              ws ,      1 ,
     &            0.D0 ,
     &             Pws ,      1   );
      call dgemv(  'T' ,
     &               n , NCOORD ,
     &            1.D0 ,
     &               P ,    LDP ,
     &             Pws ,      1 ,
     &            0.D0 ,
     &           PtPws ,      1   );
      ihl = 0;
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;
              ihl = ihl + 1;
              u(ih,il) = u(ih,il) + COMPLEX( PtPws(ihl) , 0.D0 );
          enddo
      enddo



c-----Imaginary part
      ihl = 0;
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;
              ihl = ihl + 1;
              ws(ihl) = wzwr(iabs(ih),il) * DIMAG( src(ih,il) );
          enddo
      enddo
      call dgemv(  'N' ,
     &               n , NCOORD ,
     &            1.D0 ,
     &               P ,    LDP ,
     &              ws ,      1 ,
     &            0.D0 ,
     &             Pws ,      1   );
      call dgemv(  'T' ,
     &               n , NCOORD ,
     &            1.D0 ,
     &               P ,    LDP ,
     &             Pws ,      1 ,
     &            0.D0 ,
     &           PtPws ,      1   );
      ihl = 0;
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;
              ihl = ihl + 1;
              u(ih,il) = u(ih,il) + COMPLEX( 0.D0 , PtPws(ihl) );
          enddo
      enddo



      return;
      end;
