c======================================================================c

      subroutine init_gs( lpr )

c======================================================================c

      IMPLICIT REAL*8    (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /fields/ phi(MG,4);

      CHARACTER fg_spx;
      common /simplex/ N_total         , N_blocks        ,
     &                 ia_spx(NBX)     , id_spx(NBX)     ,
     &                 nf_size(NBX)    , ng_size(NBX)    ,
     &                 nz_spx(NBSX,NBX), nr_spx(NBSX,NBX),
     &                 ml_spx(NBSX,NBX), fg_spx(NBSX,NBX);

      COMPLEX*16 v;
      common /v_matrix/ v( NBSX , NBSX , NBX , 2 );

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

      common /basis/ phi_z( -NGH:NGH , NTX ),
     &              dphi_z( -NGH:NGH , NTX ),
     &               phi_r(    1:NGL , NTX ),
     &              dphi_r(    1:NGL , NTX );

      common /gs_dens/ rhov_GS( -NGH:NGH , 1:NGL , 2 ),
     &                 rhos_GS( -NGH:NGH , 1:NGL     );

      common /gs_mesons/ sig_GS( -NGH:NGH , 1:NGL ),
     &                   ome_GS( -NGH:NGH , 1:NGL ),
     &                   rho_GS( -NGH:NGH , 1:NGL );



      CHARACTER fg1, fg2;
      COMPLEX*16 dens_mat( NBSX , NBSX , NBX , 2 );
      hbc = 197.328284D0;
      pi  = 3.14159265358979324D0;



      if(lpr) then
      write(6,*) ''
      write(6,*) '****** BEGIN init_gs() *****************************';
      write(6,*) '';
      endif






c-----Calculation of density matrix
      do it = 1 , 2
          do ib = 1 , N_blocks
              call zgemm( 'n'        ,   'c'        ,
     &                    id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                    COMPLEX( 1.D0 , 0.D0 )    ,
     &                    v(1,1,ib,it)              , NBSX       ,
     &                    v(1,1,ib,it)              , NBSX       ,
     &                    COMPLEX( 0.D0 , 0.D0 )    ,
     &                    dens_mat(1,1,ib,it)       , NBSX         );
          enddo
      enddo






c-----Calculation of Ground-State densities
      rhov_GS = 0.D0;
      rhos_GS = 0.D0;
      do it = 1 , 2
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  do ib = 1 , N_blocks
                      do i = 1 , id_spx(ib)
                          fg1 = fg_spx(i,ib);
                          ml1 = ml_spx(i,ib);
                          ii  = i-1+ia_spx(ib);

                          do j = 1 , id_spx(ib)
                              fg2 = fg_spx(j,ib);
                              ml2 = ml_spx(j,ib);
                              jj  = j-1+ia_spx(ib);

                              if( fg1.ne.fg2 ) CYCLE;
                              if( ml1.ne.ml2 ) CYCLE;

                              x = 2.D0 * DREAL(dens_mat(i,j,ib,it))
     &                                 * phi_z(ih,ii) * phi_z(ih,jj)
     &                                 * phi_r(il,ii) * phi_r(il,jj)
     &                                 / (2.D0*pi);

                              y = x;
                              if( fg1.eq.'g' .and. fg2.eq.'g' ) y = -y;

                              rhov_GS(ih,il,it) = rhov_GS(ih,il,it) + x;
                              rhos_GS(ih,il   ) = rhos_GS(ih,il   ) + y;

                          enddo
                      enddo
                  enddo

              enddo
          enddo
      enddo






c-----Setting Ground-State meson fields
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              ii = 1 + abs(ih) + il*(NGH+1);

              sig_GS(ih,il) = - g0_sig * phi(ii,1) / (m_sig/hbc)**2.D0;
              ome_GS(ih,il) = + g0_ome * phi(ii,2) / (m_ome/hbc)**2.D0;
              rho_GS(ih,il) = + g0_rho * phi(ii,4) / (m_rho/hbc)**2.D0;

          enddo
      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END init_gs() *******************************';
      write(6,*) '';
      endif

      return;
      end;
