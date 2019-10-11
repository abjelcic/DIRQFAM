c======================================================================c

      subroutine fam_dmesons( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

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

      common /gs_dens/ rhov_GS ( -NGH:NGH , 1:NGL , 2 ),
     &                 rhos_GS ( -NGH:NGH , 1:NGL     ),
     &                 rhovK_GS( -NGH:NGH , 1:NGL , 2 );

      COMPLEX*16 drho_v, drho_s, ldrho_v, ldrho_s;
      common /ind_dens/ drho_v ( -NGH:NGH , 1:NGL , 2 ),
     &                  drho_s ( -NGH:NGH , 1:NGL     ),
     &                  ldrho_v( -NGH:NGH , 1:NGL , 2 ),
     &                  ldrho_s( -NGH:NGH , 1:NGL     );

      COMPLEX*16 dj_r, dj_p, dj_z;
      common /ind_curr/ dj_r( -NGH:NGH , 1:NGL , 2 ),
     &                  dj_p( -NGH:NGH , 1:NGL , 2 ),
     &                  dj_z( -NGH:NGH , 1:NGL , 2 );

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



      COMPLEX*16 s( -NGH:NGH , 1:NGL );

      f0(x,a,b,c,d) =   a       * (1+b*(x+d)**2)   / (1+c*(x+d)**2);
      f1(x,a,b,c,d) = 2*a*(b-c) * (x+d)            / (1+c*(x+d)**2)**2;
      f2(x,a,b,c,d) = 2*a*(b-c) * (1-3*c*(x+d)**2) / (1+c*(x+d)**2)**3;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_dmesons() **************************';
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

              s(ih,il) = - g     * (drho_s(ih,il)                  )
     &                   - dg*rs * (drho_v(ih,il,1)+drho_v(ih,il,2));

          enddo
      enddo
      call KG_solve( 'sigma' , s , dsig );






c-----Calculating dome_0
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              rv = rhov_GS(ih,il,1) + rhov_GS(ih,il,2);
              x  = rv / rho_sat;

              g  = g0_ome * f0(x,a_ome,b_ome,c_ome,d_ome);
              dg = g0_ome * f1(x,a_ome,b_ome,c_ome,d_ome) / rho_sat;

              s(ih,il) = + g     * (drho_v(ih,il,1)+drho_v(ih,il,2))
     &                   + dg*rv * (drho_v(ih,il,1)+drho_v(ih,il,2));

          enddo
      enddo
      call KG_solve( 'omega0' , s , dome_0 );






c-----Calculating dome_r
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              rv = rhovK_GS(ih,il,1) + rhovK_GS(ih,il,2);
              x  = rv / rho_sat;

              g  = g0_ome * f0(x,a_ome,b_ome,c_ome,d_ome);

              s(ih,il) = + g * (dj_r(ih,il,1)+dj_r(ih,il,2));

          enddo
      enddo
      call KG_solve( 'omega' , s , dome_r );






c-----Calculating dome_p
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              rv = rhovK_GS(ih,il,1) + rhovK_GS(ih,il,2);
              x  = rv / rho_sat;

              g  = g0_ome * f0(x,a_ome,b_ome,c_ome,d_ome);

              s(ih,il) = + g * (dj_p(ih,il,1)+dj_p(ih,il,2));

          enddo
      enddo
      call KG_solve( 'omega' , s , dome_p );






c-----Calculating dome_z
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              rv = rhovK_GS(ih,il,1) + rhovK_GS(ih,il,2);
              x  = rv / rho_sat;

              g  = g0_ome * f0(x,a_ome,b_ome,c_ome,d_ome);

              s(ih,il) = + g * (dj_z(ih,il,1)+dj_z(ih,il,2));

          enddo
      enddo
      call KG_solve( 'omega' , s , dome_z );






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
      call KG_solve( 'rho0' , s , drho_0 );






c-----Calculating drho_r
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              rv = rhovK_GS(ih,il,1) + rhovK_GS(ih,il,2);
              x  = rv / rho_sat;

              g  = g0_rho * DEXP(-a_rho*(x-1.D0));

              s(ih,il) = + g * (-dj_r(ih,il,1)+dj_r(ih,il,2));

          enddo
      enddo
      call KG_solve( 'rho' , s , drho_r );






c-----Calculating drho_p
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              rv = rhovK_GS(ih,il,1) + rhovK_GS(ih,il,2);
              x  = rv / rho_sat;

              g  = g0_rho * DEXP(-a_rho*(x-1.D0));

              s(ih,il) = + g * (-dj_p(ih,il,1)+dj_p(ih,il,2));

          enddo
      enddo
      call KG_solve( 'rho' , s , drho_p );






c-----Calculating drho_z
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              rv = rhovK_GS(ih,il,1) + rhovK_GS(ih,il,2);
              x  = rv / rho_sat;

              g  = g0_rho * DEXP(-a_rho*(x-1.D0));

              s(ih,il) = + g * (-dj_z(ih,il,1)+dj_z(ih,il,2));

          enddo
      enddo
      call KG_solve( 'rho' , s , drho_z );






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_dmesons() ****************************';
      write(6,*) '';
      endif

      return;
      end;






c======================================================================c

      subroutine KG_solve( meson , src , u )

c======================================================================c
      ! Solves the Klein-Gordon equation in QHO basis expansion:
      !    ( -Delta + (mc^2/hbarc)^2 ) u(z,r,phi) = s(z,r,phi).
      !
      ! Source s(z,r,phi) is assumed to have s(z,r)*cos/sin(K*phi) form.
      !
      ! One can show that the solution u(z,r,phi) will follow the same
      ! angular form, i.e. u(z,r,phi) = u(z,r)*cos/sin(K*phi).
      !
      !
      !  [in]: meson is the name of the meson field being calculated
      !        src   is s(z,r) source term given on Gaussian mesh
      !
      ! [out]: u     is u(z,r) part of the solution on Gaussian mesh

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'

      common /mesmat/ Qsig( NHSIZE , NCOORD ), QKsig( NHSIZE , NCOORD ),
     &                Qome( NHSIZE , NCOORD ), QKome( NHSIZE , NCOORD ),
     &                Qrho( NHSIZE , NCOORD ), QKrho( NHSIZE , NCOORD ),
     &                   P( NHSIZE , NCOORD );

      CHARACTER( LEN = * ) meson;
      COMPLEX*16 src( -NGH:NGH , 1:NGL );
      COMPLEX*16 u  ( -NGH:NGH , 1:NGL );



      REAL*8    s( NCOORD );
      REAL*8   Qs( NHSIZE );
      REAL*8 PtQs( NCOORD );



      u = COMPLEX( 0.D0 , 0.D0 );






c-----Real part
      ihl = 0;
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;
              ihl = ihl + 1;
              s(ihl) = DREAL( src(ih,il) );
          enddo
      enddo

      select case( meson )
          case( 'sigma'  )
              call dgemv(  'N'  , NHSIZE , NCOORD ,
     &                    1.D0  ,
     &                    Qsig  , NHSIZE ,
     &                       s  ,      1 ,
     &                    0.D0  ,
     &                      Qs  ,      1           );
          case( 'omega0' )
              call dgemv(  'N'  , NHSIZE , NCOORD ,
     &                    1.D0  ,
     &                    Qome  , NHSIZE ,
     &                       s  ,      1 ,
     &                    0.D0  ,
     &                      Qs  ,      1           );
          case( 'omega'  )
              call dgemv(  'N'  , NHSIZE , NCOORD ,
     &                    1.D0  ,
     &                    QKome , NHSIZE ,
     &                       s  ,      1 ,
     &                    0.D0  ,
     &                      Qs  ,      1           );
          case( 'rho0'   )
              call dgemv(  'N'  , NHSIZE , NCOORD ,
     &                    1.D0  ,
     &                    Qrho  , NHSIZE ,
     &                       s  ,      1 ,
     &                    0.D0  ,
     &                      Qs  ,      1           );
          case( 'rho'    )
              call dgemv(  'N'  , NHSIZE , NCOORD ,
     &                    1.D0  ,
     &                    QKrho , NHSIZE ,
     &                       s  ,      1 ,
     &                    0.D0  ,
     &                      Qs  ,      1           );
          case default
              stop 'Error: Unknown meson!';
      end select

      call dgemv(  'T'  , NHSIZE , NCOORD ,
     &            1.D0  ,
     &               P  , NHSIZE ,
     &              Qs  ,      1 ,
     &            0.D0  ,
     &            PtQs  ,      1           );

      ihl = 0;
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;
              ihl = ihl + 1;
              u(ih,il) = u(ih,il) + COMPLEX( PtQs(ihl) , 0.D0 );
          enddo
      enddo






c-----Imaginary part
      ihl = 0;
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;
              ihl = ihl + 1;
              s(ihl) = DIMAG( src(ih,il) );
          enddo
      enddo

      select case( meson )
          case( 'sigma'  )
              call dgemv(  'N'  , NHSIZE , NCOORD ,
     &                    1.D0  ,
     &                    Qsig  , NHSIZE ,
     &                       s  ,      1 ,
     &                    0.D0  ,
     &                      Qs  ,      1           );
          case( 'omega0' )
              call dgemv(  'N'  , NHSIZE , NCOORD ,
     &                    1.D0  ,
     &                    Qome  , NHSIZE ,
     &                       s  ,      1 ,
     &                    0.D0  ,
     &                      Qs  ,      1           );
          case( 'omega'  )
              call dgemv(  'N'  , NHSIZE , NCOORD ,
     &                    1.D0  ,
     &                    QKome , NHSIZE ,
     &                       s  ,      1 ,
     &                    0.D0  ,
     &                      Qs  ,      1           );
          case( 'rho0'   )
              call dgemv(  'N'  , NHSIZE , NCOORD ,
     &                    1.D0  ,
     &                    Qrho  , NHSIZE ,
     &                       s  ,      1 ,
     &                    0.D0  ,
     &                      Qs  ,      1           );
          case( 'rho'    )
              call dgemv(  'N'  , NHSIZE , NCOORD ,
     &                    1.D0  ,
     &                    QKrho , NHSIZE ,
     &                       s  ,      1 ,
     &                    0.D0  ,
     &                      Qs  ,      1           );
          case default
              stop 'Error: Unknown meson!';
      end select

      call dgemv(  'T'  , NHSIZE , NCOORD ,
     &            1.D0  ,
     &               P  , NHSIZE ,
     &              Qs  ,      1 ,
     &            0.D0  ,
     &            PtQs  ,      1           );

      ihl = 0;
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;
              ihl = ihl + 1;
              u(ih,il) = u(ih,il) + COMPLEX( 0.D0 , PtQs(ihl) );
          enddo
      enddo






      return;
      end;
