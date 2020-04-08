c======================================================================c

      subroutine check_gs_dens( lpr )

c======================================================================c

      USE dirqfampar;
      USE simplex;
      USE v_matrix;
      USE herpol;
      USE lagpol;
      USE gaucor;
      USE gaussh;
      USE gaussl;
      USE dens;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      CHARACTER fg1, fg2;
      DOUBLE COMPLEX dens_mat( NBSX , NBSX , NBX , 2 );
      DOUBLE PRECISION RHO_v(MG,2);
      DOUBLE PRECISION RHO_s(MG,2);



      if(lpr) then
      write(6,*) ''
      write(6,*) '****** BEGIN check_gs_dens() ***********************';
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
      RHO_v = 0.D0;
      RHO_s = 0.D0;
      do it = 1 , 2
          do il = 0 , NGL
              do ih = 0 , NGH
                  ii = 1+ih + il*(NGH+1);

                  do ib = 1 , N_blocks
                      do i = 1 , id_spx(ib)
                          fg1   = fg_spx(i,ib);
                          nz1   = nz_spx(i,ib);
                          nr1   = nr_spx(i,ib);
                          ml1   = ml_spx(i,ib);
                          qhql1 = qh(nz1,ih) * ql(nr1,abs(ml1),il);
                          do j = 1 , id_spx(ib)
                              fg2   = fg_spx(j,ib);
                              nz2   = nz_spx(j,ib);
                              nr2   = nr_spx(j,ib);
                              ml2   = ml_spx(j,ib);
                              qhql2 = qh(nz2,ih) * ql(nr2,abs(ml2),il);


                              if( ml1.ne.ml2 .or. fg1.ne.fg2 ) CYCLE;


                              if( fg1.eq.'f' .and. fg2.eq.'f' ) then
                                  x = 2.D0 * DREAL(dens_mat(i,j,ib,it))
     &                                     * qhql1*qhql2;

                                  RHO_v(ii,it) = RHO_v(ii,it) + x
                                  RHO_s(ii,it) = RHO_s(ii,it) + x
                              endif

                              if( fg1.eq.'g' .and. fg2.eq.'g' ) then
                                  x = 2.D0 * DREAL(dens_mat(i,j,ib,it))
     &                                     * qhql1*qhql2;

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






c-----Compares densities
      err1 = 0.D0; eucl1 = 0.D0;
      err2 = 0.D0; eucl2 = 0.D0;
      do il = 0 , NGL
          do ih = 0 , NGH
              ii = 1+ih + il*(NGH+1);

              rho_new   = RHO_v(ii,1) + RHO_v(ii,2);
              rho_old   = ro(ii,2);
              err1  = err1  + ( rho_new - rho_old )**2.D0;
              eucl1 = eucl1 + (           rho_old )**2.D0;

              rho_new   = RHO_s(ii,1) + RHO_s(ii,2);
              rho_old   = ro(ii,1);
              err2  = err2  + ( rho_new - rho_old )**2.D0;
              eucl2 = eucl2 + (           rho_old )**2.D0;

          enddo
      enddo
      err1  = DSQRT(err1);
      err2  = DSQRT(err2);
      eucl1 = DSQRT(eucl1);
      eucl2 = DSQRT(eucl2);
      call assert( err1/eucl1 .lt. 1.D-8 , 'GS densities problem' );
      call assert( err2/eucl2 .lt. 1.D-8 , 'GS densities problem' );






c-----Prints densities
      if( lpr ) then

          write(6,*) 'Comparison of the isoscalar-vector densities';
          write(6,*) '';
          write(6,*) '         Coordinate [fm]         ',
     &               '         rho_new  ',
     &               '         rho_old  ',
     &               '         |error|  ',
     &               '     rel.error';

          write(6,*) '';

          do il = 0 , NGL
              do ih = 0 , NGH
                  ii = 1+ih + il*(NGH+1);

                  rho_new   = RHO_v(ii,1) + RHO_v(ii,2);
                  rho_old   = ro(ii,2);
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


          write(6,*) 'Comparison of the isoscalar-scalar densities';
          write(6,*) '';
          write(6,*) '         Coordinate [fm]         ',
     &               '         rho_new  ',
     &               '         rho_old  ',
     &               '         |error|  ',
     &               '     rel.error';

          write(6,*) '';

          do il = 0 , NGL
              do ih = 0 , NGH
                  ii = 1+ih + il*(NGH+1);

                  rho_new   = RHO_s(ii,1) + RHO_s(ii,2);
                  rho_old   = ro(ii,1);
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






  100 format( a, i2, a, 1f6.2, a,
     &        a, i2, a, 1f6.2, a,
     &        1f18.13,
     &        1f18.13,
     &        1e15.5,
     &        1e15.5            );


      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END check_gs_dens() *************************';
      write(6,*) '';
      endif

      return;
      end;
