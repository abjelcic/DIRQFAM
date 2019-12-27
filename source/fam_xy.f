c======================================================================c

      subroutine fam_xy( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /fam/ omega_start, omega_end, delta_omega, omega_print,
     &             omega, gamma_smear,
     &             i_calculation_type, i_coulomb, i_pairing,
     &             J_multipole, K_multipole, ISO;

      CHARACTER fg_spx;
      common /simplex/ N_total         , N_blocks        ,
     &                 ia_spx(NBX)     , id_spx(NBX)     ,
     &                 nf_size(NBX)    , ng_size(NBX)    ,
     &                 nz_spx(NBSX,NBX), nr_spx(NBSX,NBX),
     &                 ml_spx(NBSX,NBX), fg_spx(NBSX,NBX);

      common /fam_energies/ E_fam( NTX , 2 );

      COMPLEX*16 f20, f02;
      common /f02_f20_matrix/ f20( NTX , NTX , 2 ),
     &                        f02( NTX , NTX , 2 );

      COMPLEX*16 h20, h02;
      common /h20h02/ h20( NTX , NTX , 2 ),
     &                h02( NTX , NTX , 2 );

      COMPLEX*16 x_fam, y_fam;
      common /xy_fam/ x_fam( NTX , NTX , 2 ),
     &                y_fam( NTX , NTX , 2 );



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_xy() ******************************';
      write(6,*) '';
      endif


      if( i_calculation_type .eq. 0 ) then
          h20 = COMPLEX( 0.D0 , 0.D0 );
          h02 = COMPLEX( 0.D0 , 0.D0 );
      endif


c-----Calculation of xmn = -( f20mn + h20mn )/( Em + En - hw - igamma )
c-----Calculation of ymn = -( f02mn + h02mn )/( Em + En + hw + igamma )
      do it = 1 , 2
          do j = 1 , N_total
              E_mu = E_fam(j,it);
              do i = 1 , N_total
                  E_nu = E_fam(i,it);

                  x_fam(i,j,it) = - ( f20(i,j,it) + h20(i,j,it) ) /
     &                   COMPLEX( E_mu + E_nu - omega , - gamma_smear );

                  y_fam(i,j,it) = - ( f02(i,j,it) + h02(i,j,it) ) /
     &                   COMPLEX( E_mu + E_nu + omega , + gamma_smear );

              enddo
          enddo
      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_xy() ********************************';
      write(6,*) '';
      endif

      return;
      end;
