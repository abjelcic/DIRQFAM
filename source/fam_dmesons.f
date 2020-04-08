c======================================================================c

      subroutine fam_dmesons( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE ddpc1ddme2;
      USE gs_dens;
      USE mesmat;
      USE ddens;
      USE dcurr;
      USE dmesons;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      DOUBLE COMPLEX s ( -NGH:NGH , 1:NGL );
      DOUBLE COMPLEX s1( -NGH:NGH , 1:NGL );
      DOUBLE COMPLEX s2( -NGH:NGH , 1:NGL );
      DOUBLE COMPLEX s3( -NGH:NGH , 1:NGL );
      DOUBLE COMPLEX u1( -NGH:NGH , 1:NGL );
      DOUBLE COMPLEX u2( -NGH:NGH , 1:NGL );
      DOUBLE COMPLEX u3( -NGH:NGH , 1:NGL );

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

      USE dirqfampar;
      USE quadrature;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)

      DOUBLE COMPLEX src( -NGH:NGH , 1:NGL );
      DOUBLE COMPLEX u  ( -NGH:NGH , 1:NGL );
      INTEGER n;
      INTEGER LDP;
      DOUBLE PRECISION     P( LDP , NCOORD );
      DOUBLE PRECISION    ws(       NCOORD );
      DOUBLE PRECISION   Pws( LDP          );
      DOUBLE PRECISION PtPws(       NCOORD );



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
