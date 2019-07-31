c======================================================================c

      subroutine fam_dpotentials( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /dens/ ro(MG,4), dro(MG,4);
      common /gaussh/ xh(0:NGH), wh(0:NGH), zb(0:NGH);
      common /gaussl/ xl(0:NGL), wl(0:NGL), sxl(0:NGL), rb(0:NGL);

      common /DDPC1/  a_s , b_s , c_s , d_s ,
     &                a_v , b_v , c_v , d_v ,
     &                a_tv, b_tv, c_tv, d_tv,
     &                del_s,
     &                rho_sat;

      COMPLEX*16 drho_v, drho_s, ldrho_v, ldrho_s;
      common /ind_dens/ drho_v ( -NGH:NGH , 1:NGL , 2 ),
     &                  drho_s ( -NGH:NGH , 1:NGL     ),
     &                  ldrho_v( -NGH:NGH , 1:NGL , 2 ),
     &                  ldrho_s( -NGH:NGH , 1:NGL     );

      COMPLEX*16 dj_r, dj_p, dj_z;
      common /ind_curr/ dj_r( -NGH:NGH , 1:NGL , 2 ),
     &                  dj_p( -NGH:NGH , 1:NGL , 2 ),
     &                  dj_z( -NGH:NGH , 1:NGL , 2 );

      COMPLEX*16 dV_Cou;
      common /fam_coul/ dV_Cou( -NGH:NGH , 1:NGL ),
     &                       G( -NGH:NGH , 1:NGL , -NGH:NGH , 1:NGL );

      COMPLEX*16 dVpS, dVmS, dSig_z, dSig_r, dSig_p;
      common /fam_pot/ dVpS  ( -NGH:NGH , 1:NGL , 2 ),
     &                 dVmS  ( -NGH:NGH , 1:NGL , 2 ),
     &                 dSig_z( -NGH:NGH , 1:NGL , 2 ),
     &                 dSig_r( -NGH:NGH , 1:NGL , 2 ),
     &                 dSig_p( -NGH:NGH , 1:NGL , 2 );



      COMPLEX*16 Sig0  ( -NGH:NGH , 1:NGL , 2 );
      COMPLEX*16 Sig_s ( -NGH:NGH , 1:NGL     );
      COMPLEX*16 Sig0_R( -NGH:NGH , 1:NGL     );
      COMPLEX*16 Sig_z ( -NGH:NGH , 1:NGL , 2 );
      COMPLEX*16 Sig_r ( -NGH:NGH , 1:NGL , 2 );
      COMPLEX*16 Sig_p ( -NGH:NGH , 1:NGL , 2 );

      REAL*8 a0_s( -NGH:NGH , 1:NGL );
      REAL*8 a1_s( -NGH:NGH , 1:NGL );
      REAL*8 a2_s( -NGH:NGH , 1:NGL );

      REAL*8 a0_v( -NGH:NGH , 1:NGL );
      REAL*8 a1_v( -NGH:NGH , 1:NGL );
      REAL*8 a2_v( -NGH:NGH , 1:NGL );

      REAL*8 a0_tv( -NGH:NGH , 1:NGL );
      REAL*8 a1_tv( -NGH:NGH , 1:NGL );
      REAL*8 a2_tv( -NGH:NGH , 1:NGL );

      COMPLEX*16 z, w, rhos, rhov, rhotv, lrhos;
      REAL*8 hbc;

      a0(x,a,b,c,d) = a + ( b + c*x )*DEXP(-d*x);
      a1(x,a,b,c,d) = ( c - d*(b+c*x) )*DEXP(-d*x);
      a2(x,a,b,c,d) = ( -2*c*d + d*d*(b+c*x) )*DEXP(-d*x);

      hbc = 197.328284D0;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_dpotentials() **********************';
      write(6,*) '';
      endif






c-----Calculation of density-dependent coupling "constants"
      do il = 1 , NGL
          do ih = -NGH , NGH
              if( ih .eq. 0 ) CYCLE;

              rs  = ro( 1+abs(ih) + il*(NGH+1) , 1 ) / rho_sat;
              rv  = ro( 1+abs(ih) + il*(NGH+1) , 2 ) / rho_sat;
              rtv = ro( 1+abs(ih) + il*(NGH+1) , 4 ) / rho_sat;

              a0_s(ih,il)  = a0(rv, a_s , b_s , c_s , d_s );
              a1_s(ih,il)  = a1(rv, a_s , b_s , c_s , d_s ) * rs;
              a2_s(ih,il)  = a2(rv, a_s , b_s , c_s , d_s ) * rs*rs;

              a0_v(ih,il)  = a0(rv, a_v , b_v , c_v , d_v );
              a1_v(ih,il)  = a1(rv, a_v , b_v , c_v , d_v ) * rv;
              a2_v(ih,il)  = a2(rv, a_v , b_v , c_v , d_v ) * rv*rv;

              a0_tv(ih,il) = a0(rv, a_tv, b_tv, c_tv, d_tv);
              a1_tv(ih,il) = a1(rv, a_tv, b_tv, c_tv, d_tv) * rtv;
              a2_tv(ih,il) = a2(rv, a_tv, b_tv, c_tv, d_tv) * rtv*rtv;

          enddo
      enddo






c-----Calculation of Sig0
      do it = 1 , 2
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;

                  rhov  = + drho_v(ih,il,1) + drho_v(ih,il,2);
                  rhotv = - drho_v(ih,il,1) + drho_v(ih,il,2);

                  z = ( a0_v(ih,il) + a1_v(ih,il) )*rhov;
                  w = a0_tv(ih,il)*rhotv + a1_tv(ih,il)*rhov;

                  if( it .eq. 1 ) then
                      w = - w;
                  endif

                  Sig0(ih,il,it) = hbc * ( z + w );

              enddo
          enddo
      enddo






c-----Calculation of Sig0_R
      do il = 1 , NGL
          do ih = -NGH , NGH
              if( ih .eq. 0 ) CYCLE;

              rhos  = + drho_s(ih,il);
              rhov  = + drho_v(ih,il,1) + drho_v(ih,il,2);
              rhotv = - drho_v(ih,il,1) + drho_v(ih,il,2);

              z = 0.5D0*( a2_s(ih,il)+a2_v(ih,il)+a2_tv(ih,il) ) * rhov;
              z = z + a1_s(ih,il) *rhos;
              z = z + a1_v(ih,il) *rhov;
              z = z + a1_tv(ih,il)*rhotv;

              Sig0_R(ih,il) = hbc * z;

          enddo
      enddo






c-----Calculation of Sig_s
      do il = 1 , NGL
          do ih = -NGH , NGH
              if( ih .eq. 0 ) CYCLE;

              rhos  = +  drho_s(ih,il);
              rhov  = +  drho_v(ih,il,1) +  drho_v(ih,il,2);
              lrhos = + ldrho_s(ih,il);

              z = a1_s(ih,il)*rhov + a0_s(ih,il)*rhos + del_s*lrhos;

              Sig_s(ih,il) = hbc * z;

          enddo
      enddo






c-----Calculation of Sig_z
      do it = 1 , 2
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;

                  z = a0_v(ih,il) *( + dj_z(ih,il,1) + dj_z(ih,il,2) );
                  w = a0_tv(ih,il)*( - dj_z(ih,il,1) + dj_z(ih,il,2) );

                  if( it .eq. 1 ) then
                      w = - w;
                  endif

                  Sig_z(ih,il,it) = hbc * ( z + w );

              enddo
          enddo
      enddo






c-----Calculation of Sig_r
      do it = 1 , 2
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;

                  z = a0_v(ih,il) *( + dj_r(ih,il,1) + dj_r(ih,il,2) );
                  w = a0_tv(ih,il)*( - dj_r(ih,il,1) + dj_r(ih,il,2) );

                  if( it .eq. 1 ) then
                      w = - w;
                  endif

                  Sig_r(ih,il,it) = hbc * ( z + w );

              enddo
          enddo
      enddo






c-----Calculation of Sig_p
      do it = 1 , 2
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;

                  z = a0_v(ih,il) *( + dj_p(ih,il,1) + dj_p(ih,il,2) );
                  w = a0_tv(ih,il)*( - dj_p(ih,il,1) + dj_p(ih,il,2) );

                  if( it .eq. 1 ) then
                      w = - w;
                  endif

                  Sig_p(ih,il,it) = hbc * ( z + w );

              enddo
          enddo
      enddo






c-----Calculation of the induced potentials
      call fam_dcoulomb( .false. );

      do it = 1 , 2
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;

                  z = Sig0(ih,il,it) + Sig0_R(ih,il);
                  w = Sig_s(ih,il);

                  dVpS  (ih,il,it) = z + w;
                  dVmS  (ih,il,it) = z - w;
                  dSig_z(ih,il,it) = Sig_z(ih,il,it);
                  dSig_r(ih,il,it) = Sig_r(ih,il,it);
                  dSig_p(ih,il,it) = Sig_p(ih,il,it);

                  if( it .eq. 2 ) then
                      dVpS(ih,il,it) = dVpS(ih,il,it) + dV_Cou(ih,il);
                      dVmS(ih,il,it) = dVmS(ih,il,it) + dV_Cou(ih,il);
                  endif

              enddo
          enddo
      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_dpotentials() ************************';
      write(6,*) '';
      endif

      return;
      end;
