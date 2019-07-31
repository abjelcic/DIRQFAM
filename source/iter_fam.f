c======================================================================c

      subroutine iter_fam( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /fam/ omega_start, omega_end, delta_omega, omega_print,
     &             omega, gamma_smear,
     &             i_calculation_type, i_coulomb, i_pairing,
     &             J_multipole, K_multipole, ISO;

      common /fam_iter/ error, tol, iter, iter_max;



      CHARACTER*27 text1;
      CHARACTER*25 text2;
      CHARACTER*21 text3;

      tol      = 1.D-5;  ! Broyden error tolerance
      iter     = 0;      ! Index of current iteration
      iter_max = 999;    ! Maximum number of iterations
      error    = 1.D+10; ! Error in current iteration

      text1    = 'ITERATION INTERRUPTED AFTER';
      text2    = 'ITERATION CONVERGED AFTER';
      text3    = ' STEPS, FINAL error =';



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN iter_fam() ****************************';
      write(6,*) '';
      endif






      write(6,*) '';
      write(6,'(2a,1f6.3,a)') 'START ITERATION ',
     &                        '---------------------> omega = ', omega ,
     &                        '[MeV/hbar]';






      ! Initializing first Broyden vector
      call fam_broyden( .false. , .true. );

c-----Main FAM iteration for fixed frequency
      do iter = 1 , iter_max

          write(6,'(i3,a,1f10.5,2a)') iter, '.Iteration, error = ',
     &                                error,'                    ',
     &                                      '         |';
          call flush(6);



          call fam_drhodkappa ( .false. );
          call fam_ddensdcurr ( .false. );
          call fam_dpotentials( .false. );
          call fam_dh1        ( .false. );
          call fam_ddelta     ( .false. );
          call fam_broyden    ( .false. , .false. );



          if( iter.gt.3 .and. error.lt.tol ) then
              write(6,'(a,i4,a,e13.5)') text2, iter, text3, error;
              write(6,*) '';
              EXIT;
          endif
          if( iter .eq. iter_max ) then
              write(6,'(a,i4,a,e11.5)') text1, iter, text3, error;
              EXIT;
          endif

      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END iter_fam() ******************************';
      write(6,*) '';
      endif

      return;
      end;
