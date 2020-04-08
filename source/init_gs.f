c======================================================================c

      subroutine init_gs( lpr )

c======================================================================c

      USE dirqfampar;
      USE fields; !from ground state code
      USE simplex;
      USE v_matrix;
      USE ddpc1ddme2;
      USE basis;
      USE gs_dens;
      USE gs_mesons;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      CHARACTER fg1, fg2;
      DOUBLE COMPLEX dens_mat( NBSX , NBSX , NBX , 2 );
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
