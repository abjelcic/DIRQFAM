c======================================================================c

      subroutine fam_dpotentials( lpr )

c======================================================================c

      USE dirqfampar;
      USE ddpc1ddme2;
      USE gs_dens;
      USE gs_mesons;
      USE ddens;
      USE dcurr;
      USE dlaplace;
      USE dmesons;
      USE dcoul;
      USE dpot;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;

      CHARACTER*10 parname;
      common /partyp/ parname;



      DOUBLE COMPLEX Sig0  ( -NGH:NGH , 1:NGL , 2 );
      DOUBLE COMPLEX Sig0_R( -NGH:NGH , 1:NGL     );
      DOUBLE COMPLEX Sig_s ( -NGH:NGH , 1:NGL     );
      DOUBLE COMPLEX Sig_z ( -NGH:NGH , 1:NGL , 2 );
      DOUBLE COMPLEX Sig_r ( -NGH:NGH , 1:NGL , 2 );
      DOUBLE COMPLEX Sig_p ( -NGH:NGH , 1:NGL , 2 );

      DOUBLE PRECISION a0_s ( -NGH:NGH , 1:NGL );
      DOUBLE PRECISION a1_s ( -NGH:NGH , 1:NGL );
      DOUBLE PRECISION a2_s ( -NGH:NGH , 1:NGL );

      DOUBLE PRECISION a0_v ( -NGH:NGH , 1:NGL );
      DOUBLE PRECISION a1_v ( -NGH:NGH , 1:NGL );
      DOUBLE PRECISION a2_v ( -NGH:NGH , 1:NGL );

      DOUBLE PRECISION a0_tv( -NGH:NGH , 1:NGL );
      DOUBLE PRECISION a1_tv( -NGH:NGH , 1:NGL );
      DOUBLE PRECISION a2_tv( -NGH:NGH , 1:NGL );

      DOUBLE COMPLEX z, w, drhos, drhov, drhotv, ldrhos;

      a0(x,a,b,c,d) =              a + (b+c*x)*DEXP(-d*x);
      a1(x,a,b,c,d) =        ( c - d*(b+c*x) )*DEXP(-d*x);
      a2(x,a,b,c,d) = ( -2*c*d + d*d*(b+c*x) )*DEXP(-d*x);

      f0(x,a,b,c,d) =   a       * (1+b*(x+d)**2)   / (1+c*(x+d)**2);
      f1(x,a,b,c,d) = 2*a*(b-c) * (x+d)            / (1+c*(x+d)**2)**2;
      f2(x,a,b,c,d) = 2*a*(b-c) * (1-3*c*(x+d)**2) / (1+c*(x+d)**2)**3;

      hbc = 197.328284D0;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_dpotentials() *********************';
      write(6,*) '';
      endif






      if( parname .eq. 'DD-PC1' ) then



c---------Calculation of density-dependent coupling "constants"
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  rs  =  +rhos_GS(ih,il)                     / rho_sat;
                  rv  = (+rhov_GS(ih,il,1)+rhov_GS(ih,il,2)) / rho_sat;
                  rtv = (-rhov_GS(ih,il,1)+rhov_GS(ih,il,2)) / rho_sat;

                  a0_s(ih,il)  = a0(rv, a_s , b_s , c_s , d_s );
                  a1_s(ih,il)  = a1(rv, a_s , b_s , c_s , d_s )*rs;
                  a2_s(ih,il)  = a2(rv, a_s , b_s , c_s , d_s )*rs*rs;

                  a0_v(ih,il)  = a0(rv, a_v , b_v , c_v , d_v );
                  a1_v(ih,il)  = a1(rv, a_v , b_v , c_v , d_v )*rv;
                  a2_v(ih,il)  = a2(rv, a_v , b_v , c_v , d_v )*rv*rv;

                  a0_tv(ih,il) = a0(rv, a_tv, b_tv, c_tv, d_tv);
                  a1_tv(ih,il) = a1(rv, a_tv, b_tv, c_tv, d_tv)*rtv;
                  a2_tv(ih,il) = a2(rv, a_tv, b_tv, c_tv, d_tv)*rtv*rtv;

              enddo
          enddo



c---------Calculation of Sig0
          do it = 1 , 2
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih .eq. 0 ) CYCLE;

                      drhov  = + drho_v(ih,il,1) + drho_v(ih,il,2);
                      drhotv = - drho_v(ih,il,1) + drho_v(ih,il,2);

                      z = ( a0_v(ih,il) + a1_v(ih,il) )*drhov;
                      w = a0_tv(ih,il)*drhotv + a1_tv(ih,il)*drhov;

                      if( it .eq. 1 ) then
                          w = - w;
                      endif

                      Sig0(ih,il,it) = hbc * ( z + w );

                  enddo
              enddo
          enddo



c---------Calculation of Sig0_R
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  drhos  = + drho_s(ih,il);
                  drhov  = + drho_v(ih,il,1) + drho_v(ih,il,2);
                  drhotv = - drho_v(ih,il,1) + drho_v(ih,il,2);

                  z = (a2_s(ih,il)+a2_v(ih,il)+a2_tv(ih,il))/2.D0*drhov;
                  z = z + a1_s (ih,il) * drhos;
                  z = z + a1_v (ih,il) * drhov;
                  z = z + a1_tv(ih,il) * drhotv;

                  Sig0_R(ih,il) = hbc * z;

              enddo
          enddo



c---------Calculation of Sig_s
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  drhos  = +  drho_s(ih,il);
                  drhov  = +  drho_v(ih,il,1) +  drho_v(ih,il,2);
                  ldrhos = + ldrho_s(ih,il);

                  z = + a1_s(ih,il)*drhov
     &                + a0_s(ih,il)*drhos
     &                + del_s*ldrhos;

                  Sig_s(ih,il) = hbc * z;

              enddo
          enddo



c---------Calculation of Sig_z
          do it = 1 , 2
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih .eq. 0 ) CYCLE;

                      z = a0_v (ih,il) * (+dj_z(ih,il,1)+dj_z(ih,il,2));
                      w = a0_tv(ih,il) * (-dj_z(ih,il,1)+dj_z(ih,il,2));

                      if( it .eq. 1 ) then
                          w = - w;
                      endif

                      Sig_z(ih,il,it) = hbc * ( z + w );

                  enddo
              enddo
          enddo



c---------Calculation of Sig_r
          do it = 1 , 2
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih .eq. 0 ) CYCLE;

                      z = a0_v (ih,il) * (+dj_r(ih,il,1)+dj_r(ih,il,2));
                      w = a0_tv(ih,il) * (-dj_r(ih,il,1)+dj_r(ih,il,2));

                      if( it .eq. 1 ) then
                          w = - w;
                      endif

                      Sig_r(ih,il,it) = hbc * ( z + w );

                  enddo
              enddo
          enddo



c---------Calculation of Sig_p
          do it = 1 , 2
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih .eq. 0 ) CYCLE;

                      z = a0_v (ih,il) * (+dj_p(ih,il,1)+dj_p(ih,il,2));
                      w = a0_tv(ih,il) * (-dj_p(ih,il,1)+dj_p(ih,il,2));

                      if( it .eq. 1 ) then
                          w = - w;
                      endif

                      Sig_p(ih,il,it) = hbc * ( z + w );

                  enddo
              enddo
          enddo



      endif






      if( parname .eq. 'DD-ME2' ) then

          call fam_dmesons( .false. );



c---------Calculation of Sig0
          do it = 1 , 2
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih .eq. 0 ) CYCLE;

                      rv = rhov_GS(ih,il,1) + rhov_GS(ih,il,2);
                      x  = rv/rho_sat;

                       g_ome = g0_ome * f0(x,a_ome,b_ome,c_ome,d_ome);
                      dg_ome = g0_ome * f1(x,a_ome,b_ome,c_ome,d_ome)
     &                                  / rho_sat;
                       g_rho = g0_rho * DEXP(-a_rho*(x-1.D0));
                      dg_rho = g0_rho * DEXP(-a_rho*(x-1.D0))
     &                                * (-a_rho/rho_sat);

                      drhov = + drho_v(ih,il,1) + drho_v(ih,il,2);

                      z = +  g_ome * dome_0(ih,il)
     &                    + dg_ome * ome_GS(ih,il) * drhov;

                      w = +  g_rho * drho_0(ih,il)
     &                    + dg_rho * rho_GS(ih,il) * drhov;


                      if( it .eq. 1 ) then
                          w = - w;
                      endif

                      Sig0(ih,il,it) = hbc * ( z + w );

                  enddo
              enddo
          enddo



c---------Calculation of Sig0_R
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  rs  = + rhos_GS(ih,il);
                  rv  = + rhov_GS(ih,il,1) + rhov_GS(ih,il,2);
                  rtv = - rhov_GS(ih,il,1) + rhov_GS(ih,il,2);
                  x   = rv/rho_sat;

                    g_sig = g0_sig * f0(x,a_sig,b_sig,c_sig,d_sig);
                   dg_sig = g0_sig * f1(x,a_sig,b_sig,c_sig,d_sig)
     &                               / rho_sat;
                  ddg_sig = g0_sig * f2(x,a_sig,b_sig,c_sig,d_sig)
     &                               / rho_sat**2.D0;
                    g_ome = g0_ome * f0(x,a_ome,b_ome,c_ome,d_ome);
                   dg_ome = g0_ome * f1(x,a_ome,b_ome,c_ome,d_ome)
     &                               / rho_sat;
                  ddg_ome = g0_ome * f2(x,a_ome,b_ome,c_ome,d_ome)
     &                               / rho_sat**2.D0;
                    g_rho = g0_rho * DEXP(-a_rho*(x-1.D0));
                   dg_rho = g0_rho * DEXP(-a_rho*(x-1.D0))
     &                                 * (-a_rho/rho_sat);
                  ddg_rho = g0_rho * DEXP(-a_rho*(x-1.D0))
     &                                 * (-a_rho/rho_sat)**2.D0;

                  drhos  = + drho_s(ih,il);
                  drhov  = + drho_v(ih,il,1) + drho_v(ih,il,2);
                  drhotv = - drho_v(ih,il,1) + drho_v(ih,il,2);


                  z = + ( + ddg_sig * sig_GS(ih,il) * rs
     &                    + ddg_ome * ome_GS(ih,il) * rv
     &                    + ddg_rho * rho_GS(ih,il) * rtv ) * drhov
     &
     &                + ( dg_sig *   dsig(ih,il) * rs     )
     &                + ( dg_ome * dome_0(ih,il) * rv     )
     &                + ( dg_rho * drho_0(ih,il) * rtv    )
     &
     &                + ( dg_sig * sig_GS(ih,il) * drhos  )
     &                + ( dg_ome * ome_GS(ih,il) * drhov  )
     &                + ( dg_rho * rho_GS(ih,il) * drhotv );


                  Sig0_R(ih,il) =  hbc * z;

              enddo
          enddo



c---------Calculation of Sig_s
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  rv = rhov_GS(ih,il,1) + rhov_GS(ih,il,2);
                  x  = rv/rho_sat;

                   g_sig = g0_sig * f0(x,a_sig,b_sig,c_sig,d_sig);
                  dg_sig = g0_sig * f1(x,a_sig,b_sig,c_sig,d_sig)
     &                              / rho_sat;

                  drhov = + drho_v(ih,il,1) + drho_v(ih,il,2);

                  z = +  g_sig *   dsig(ih,il)
     &                + dg_sig * sig_GS(ih,il) * drhov;


                  Sig_s(ih,il) = hbc * z;

              enddo
          enddo



c---------Calculation of Sig_z
          do it = 1 , 2
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih .eq. 0 ) CYCLE;

                      rv = rhov_GS(ih,il,1) + rhov_GS(ih,il,2);
                      x  = rv/rho_sat;

                      g_ome = g0_ome * f0(x,a_ome,b_ome,c_ome,d_ome);
                      g_rho = g0_rho * DEXP(-a_rho*(x-1.D0));

                      z = g_ome * dome_z(ih,il);
                      w = g_rho * drho_z(ih,il);

                      if( it .eq. 1 ) then
                          w = - w;
                      endif

                      Sig_z(ih,il,it) = hbc * ( z + w );

                  enddo
              enddo
          enddo



c---------Calculation of Sig_r
          do it = 1 , 2
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih .eq. 0 ) CYCLE;

                      rv = rhov_GS(ih,il,1) + rhov_GS(ih,il,2);
                      x  = rv/rho_sat;

                      g_ome = g0_ome * f0(x,a_ome,b_ome,c_ome,d_ome);
                      g_rho = g0_rho * DEXP(-a_rho*(x-1.D0));

                      z = g_ome * dome_r(ih,il);
                      w = g_rho * drho_r(ih,il);

                      if( it .eq. 1 ) then
                          w = - w;
                      endif

                      Sig_r(ih,il,it) = hbc * ( z + w );

                  enddo
              enddo
          enddo



c---------Calculation of Sig_p
          do it = 1 , 2
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih .eq. 0 ) CYCLE;

                      rv = rhov_GS(ih,il,1) + rhov_GS(ih,il,2);
                      x  = rv/rho_sat;

                      g_ome = g0_ome * f0(x,a_ome,b_ome,c_ome,d_ome);
                      g_rho = g0_rho * DEXP(-a_rho*(x-1.D0));

                      z = g_ome * dome_p(ih,il);
                      w = g_rho * drho_p(ih,il);

                      if( it .eq. 1 ) then
                          w = - w;
                      endif

                      Sig_p(ih,il,it) = hbc * ( z + w );

                  enddo
              enddo
          enddo



      endif






c-----Calculation of the induced potentials
      do it = 1 , 2
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  z = Sig0(ih,il,it) + Sig0_R(ih,il);
                  w = Sig_s(ih,il);

                  dVpS  (ih,il,it) = z + w;
                  dVmS  (ih,il,it) = z - w;
                  dSig_z(ih,il,it) = Sig_z(ih,il,it);
                  dSig_r(ih,il,it) = Sig_r(ih,il,it);
                  dSig_p(ih,il,it) = Sig_p(ih,il,it);

              enddo
          enddo
      enddo

      call fam_dcoulomb( .false. );
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              dVpS  (ih,il,2) = dVpS  (ih,il,2) + dVCou_0(ih,il);
              dVmS  (ih,il,2) = dVmS  (ih,il,2) + dVCou_0(ih,il);
              dSig_z(ih,il,2) = dSig_z(ih,il,2) + dVCou_z(ih,il);
              dSig_r(ih,il,2) = dSig_r(ih,il,2) + dVCou_r(ih,il);
              dSig_p(ih,il,2) = dSig_p(ih,il,2) + dVCou_p(ih,il);

          enddo
      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_dpotentials() ***********************';
      write(6,*) '';
      endif

      return;
      end;
