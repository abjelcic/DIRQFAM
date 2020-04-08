c======================================================================c

      subroutine fam_xy( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE simplex;
      USE fam_energies;
      USE f02f20matrix;
      USE h02h20matrix;
      USE xyfam;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



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
