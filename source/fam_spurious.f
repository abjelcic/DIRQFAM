c======================================================================c

      subroutine fam_spurious( lpr )

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

      COMPLEX*16 x_fam, y_fam;
      common /xy_fam/ x_fam( NTX , NTX , 2 ),
     &                y_fam( NTX , NTX , 2 );

      COMPLEX*16 r20, p20, RcmPcm_commutator, lamR, lamP;
      common /spurious/ r20( NTX , NTX , 2 ),
     &                  p20( NTX , NTX , 2 ),
     &                  RcmPcm_commutator(2),
     &                  lamR, lamP;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_spurious() ************************';
      write(6,*) '';
      endif


      if( J_multipole.ne.1 .and. J_multipole.ne.3 ) then
          return;
      endif
      if( K_multipole.ne.0 .and. K_multipole.ne.1 ) then
          return;
      endif


      lamR = COMPLEX( 0.D0 , 0.D0 );
      lamP = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2
          do j = 1 , N_total
              do i = 1 , N_total
                  lamR = lamR + ( - x_fam(i,j,it) + y_fam(i,j,it) )
     &                        * DCONJG( p20(i,j,it) );
                  lamP = lamP + ( + x_fam(i,j,it) + y_fam(i,j,it) )
     &                        * DCONJG( r20(i,j,it) );
              enddo
          enddo
      enddo
      lamR = lamR / ( RcmPcm_commutator(1)+RcmPcm_commutator(2) );
      lamP = lamP / ( RcmPcm_commutator(1)+RcmPcm_commutator(2) );


      do it = 1 , 2
          do j = 1 , N_total
              do i = 1 , N_total
                  x_fam(i,j,it) = x_fam(i,j,it) - lamR*r20(i,j,it)
     &                                          - lamP*p20(i,j,it);
                  y_fam(i,j,it) = y_fam(i,j,it) + lamR*r20(i,j,it)
     &                                          - lamP*p20(i,j,it);
              enddo
          enddo
      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_spurious() **************************';
      write(6,*) '';
      endif

      return;
      end;
