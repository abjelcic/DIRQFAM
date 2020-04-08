c======================================================================c

      subroutine iter_fam( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE fam_iter;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      CHARACTER*28 text1;
      CHARACTER*26 text2;
      CHARACTER*21 text3;
      CHARACTER*13 format1;
      CHARACTER*18 format2;
#ifdef DEBUG
      format1   = '(i3,a,1f14.9)';
      format2   = '(2x,a,1e14.9,1x,a)';
#else
      format1   = '(i3,a,1f10.5)';
      format2   = '(6x,a,1e14.9,1x,a)';
#endif

      tol      = tol;    ! In order to change tolerance, see start_fam.f
      iter     = 0;      ! Index of current iteration
      iter_max = 999;    ! Maximum number of iterations
      error    = 1.D+99; ! Error in current iteration

      text1    = 'ITERATIONS INTERRUPTED AFTER';
      text2    = 'ITERATIONS CONVERGED AFTER';
      text3    = ' STEPS, FINAL error =';



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN iter_fam() ****************************';
      write(6,*) '';
      endif






      write(6,*) '';
      write(6,'(2a,1f6.3,a)') 'START ITERATIONS ',
     &                        '---------------------> omega = ', omega ,
     &                        '[MeV/hbar]';






      ! Initializing first Broyden's vector
      call fam_broyden( .false. , .true. );

c-----Main QFAM iteration for fixed energy
      do iter = 1 , iter_max

          write(6,format1,advance='no')
     &    iter , '.Iteration, error = ' , error;
          call flush(6);


          call fam_h20h02     ( .false. );
          call fam_xy         ( .false. );
          call fam_drho       ( .false. );
          call fam_dkappa     ( .false. );
          call fam_ddensdcurr ( .false. );
          call fam_dpotentials( .false. );
          call fam_dh1        ( .false. );
          call fam_ddelta     ( .false. );
          call fam_broyden    ( .false. , .false. );

          Sn = fam_strength( .false. , 1 );
          Sp = fam_strength( .false. , 2 );


          write(6,format2,advance='yes')
     &    'S(f,w) = ' , Sn + Sp , '|';
          call flush(6);



          if( iter.gt.3 .and. error.lt.tol ) then
              write(6,'(a,i4,a,e13.5,/)') text2, iter, text3, error;
              EXIT;
          endif
          if( iter .eq. iter_max ) then
              write(6,'(a,i4,a,e11.3,/)') text1, iter, text3, error;
              EXIT;
          endif

      enddo






      if( J_multipole.eq.1 .or. J_multipole.eq.3 ) then
      if( K_multipole.eq.0 .or. K_multipole.eq.1 ) then

          ! Elimination of the spurious Nambu-Goldstone mode
          call fam_spurious( .false. );

          ! Update drho matrix since we only need
          ! drho(omega) to calculate S(f,omega)
          call fam_drho( .false. );

          ! Densities and currents also have
          ! to be updated before printing
          if( i_calculation_type .eq. 2 ) then
              call fam_ddensdcurr( .false. );
          endif

      endif
      endif






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END iter_fam() ******************************';
      write(6,*) '';
      endif

      return;
      end;
