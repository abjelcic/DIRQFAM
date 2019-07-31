c======================================================================c

      subroutine check_gs_dens( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /herpol/ qh(0:NZX,0:NGH), qh1(0:NZX,0:NGH);
      common /lagpol/ ql(0:2*NRX,0:MLX,0:NGL), ql1(0:2*NRX,0:MLX,0:NGL);
      common /gaucor/ ww(MG);
      common /gaussh/ xh(0:NGH), wh(0:NGH), zb(0:NGH);
      common /gaussl/ xl(0:NGL), wl(0:NGL), sxl(0:NGL), rb(0:NGL);
      common /dens/   ro(MG,4), dro(MG,4);

      CHARACTER fg_spx;
      common /simplex/ N_total         , N_blocks        ,
     &                 ia_spx(NBX)     , id_spx(NBX)     ,
     &                 nf_size(NBX)    , ng_size(NBX)    ,
     &                 nz_spx(NBSX,NBX), nr_spx(NBSX,NBX),
     &                 ml_spx(NBSX,NBX), fg_spx(NBSX,NBX);

      COMPLEX*16 v;
      common /v_matrix/ v( NBSX , NBSX , NBX , 2 );



      CHARACTER fg1, fg2;
      COMPLEX*16 density_matrix( NBSX , NBSX , NBX , 2 );
      REAL*8 max_error;
      REAL*8 RHO_v(MG,2);
      REAL*8 RHO_s(MG,2);



      ihl(ih,il) = 1+ih + il*(NGH+1);
      RHO_v = 0.D0;
      RHO_s = 0.D0;



      if(lpr) then
      write(6,*) ''
      write(6,*) '****** BEGIN check_gs_dens() ************************';
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
     &                    density_matrix(1,1,ib,it) , NBSX         );
          enddo
      enddo






c-----Calculation of Ground-State densities
      do it = 1 , 2
          do il = 0 , NGL
              do ih = 0 , NGH
                  ii = ihl(ih,il);
                  do ib = 1 , N_blocks
                      do i = 1 , id_spx(ib)
                          fg1 = fg_spx(i,ib);
                          nz1 = nz_spx(i,ib);
                          nr1 = nr_spx(i,ib);
                          ml1 = ml_spx(i,ib);
                          qhql1 = qh(nz1,ih) * ql(nr1,abs(ml1),il);
                          do j = 1 , id_spx(ib)
                              fg2 = fg_spx(j,ib);
                              nz2 = nz_spx(j,ib);
                              nr2 = nr_spx(j,ib);
                              ml2 = ml_spx(j,ib);
                              qhql2 = qh(nz2,ih) * ql(nr2,abs(ml2),il);


                              if( ml1.ne.ml2 .or. fg1.ne.fg2 ) CYCLE;


                              if( fg1.eq.'f' .and. fg2.eq.'f' ) then
                                  x = DREAL(density_matrix(i,j,ib,it));
                                  x = x * 2.D0*qhql1*qhql2;

                                  RHO_v(ii,it) = RHO_v(ii,it) + x
                                  RHO_s(ii,it) = RHO_s(ii,it) + x
                              endif

                              if( fg1.eq.'g' .and. fg2.eq.'g' ) then
                                  x = DREAL(density_matrix(i,j,ib,it));
                                  x = x * 2.D0*qhql1*qhql2;

                                  RHO_v(ii,it) = RHO_v(ii,it) + x
                                  RHO_s(ii,it) = RHO_s(ii,it) - x
                              endif

                          enddo
                      enddo
                  enddo

                  RHO_v(ii,it) = RHO_v(ii,it) / ww(ii);
                  RHO_s(ii,it) = RHO_s(ii,it) / ww(ii);
              enddo
          enddo
      enddo






c-----Prints densities
      if( lpr ) then

          write(6,*) '';
          write(6,*) 'Comparison of the isoscalar-vector densities';
          write(6,*) '';

          write(6,*) '         Coordinate [fm]          ',
     &               '        rho_new   ',
     &               '        rho_old   ',
     &               '        |error|  ',
     &               '     rel.error';

          write(6,*) '';

          do il = 0 , NGL
              do ih = 0 , NGH
                  rho_new   = RHO_v(ihl(ih,il),1) + RHO_v(ihl(ih,il),2);
                  rho_old   = ro(ihl(ih,il),2);
                  abs_error = DABS( rho_new - rho_old );
                  rel_error = abs_error / DABS(rho_old);

                  write(6,100)   '( r(', il, ') =', rb(il), ',  ',
     &                           '  z(', ih, ') =', zb(ih), ' )  ',
     &                           rho_new,
     &                           rho_old,
     &                           abs_error,
     &                           rel_error;

              enddo
          enddo


          write(6,*) '';
          write(6,*) 'Comparison of the isoscalar-scalar densities';
          write(6,*) '';

          write(6,*) '         Coordinate [fm]          ',
     &               '        rho_new   ',
     &               '        rho_old   ',
     &               '        |error|  ',
     &               '     rel.error';

          write(6,*) '';

          do il = 0 , NGL
              do ih = 0 , NGH
                  rho_new   = RHO_s(ihl(ih,il),1) + RHO_s(ihl(ih,il),2);
                  rho_old   = ro(ihl(ih,il),1);
                  abs_error = DABS( rho_new - rho_old );
                  rel_error = abs_error / DABS(rho_old);

                  write(6,100)   '( r(', il, ') =', rb(il), ',  ',
     &                           '  z(', ih, ') =', zb(ih), ' )  ',
     &                           rho_new,
     &                           rho_old,
     &                           abs_error,
     &                           rel_error;

              enddo
          enddo

      endif






c-----Compares densities
      max_error = -1.D0;
      do il = 0 , NGL
          do ih = 0 , NGH

              rho_new   = RHO_v(ihl(ih,il),1) + RHO_v(ihl(ih,il),2);
              rho_old   = ro(ihl(ih,il),2);
              rel_error = DABS(rho_new-rho_old)/DABS(rho_old);
              max_error = MAX( max_error , rel_error );

              rho_new   = RHO_s(ihl(ih,il),1) + RHO_s(ihl(ih,il),2);
              rho_old   = ro(ihl(ih,il),1);
              rel_error = DABS(rho_new-rho_old)/DABS(rho_old);
              max_error = MAX( max_error , rel_error );

          enddo
      enddo
      if( max_error .gt. 1.D-5 ) then
          write(6,*) 'max_error = ', max_error;
          stop 'Error: Ground-State densities wrong!';
      endif
      if( lpr ) then
          write(6,*) 'max_error = ', max_error;
      endif






  100 format( a, i2, a, 1f6.2, a,
     &        a, i2, a, 1f6.2, a,
     &        1f18.13,
     &        1f18.13,
     &        1e15.5,
     &        1e15.5            );

      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END check_gs_dens() **************************';
      write(6,*) '';
      endif

      return;
      end;
