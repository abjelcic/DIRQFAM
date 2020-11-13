c======================================================================c

      subroutine start_fam( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam_iter;
      USE fam;
      USE dh;
      USE dDelta;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN start_fam() ***************************';
      write(6,*) '';
      endif






#ifdef DEBUG
      tol = 1.D-9;  ! Self-consistency tolerance
#else
      tol = 1.D-5;  ! Self-consistency tolerance
#endif






c-----Printing strength.out header
      open( tape_strength , file   = './output/QFAM_output/strength.out'
     &                    , status = 'unknown' );

      write( tape_strength , '(a)' , advance = 'no' )
     &      'S(f,omega) = -1/pi * Im[Tr[hermconj(f)*drho(omega)]] ';

      ! Unit of measurement of S(f,omega)
      select case( J_multipole )
          case( 0 )
              write(tape_strength,*) '[fm^4/MeV]';
          case( 1 )
              select case( ISO )
                  case( 0 )
                      write(tape_strength,*) '[fm^2/MeV]';
                  case( 1 )
                      write(tape_strength,*) '[e^2fm^2/MeV]';
                  case default
                      stop 'Error: Wrong ISO!';
              end select
          case( 2 )
              write(tape_strength,*) '[fm^4/MeV]';
          case( 3 )
              write(tape_strength,*) '[fm^6/MeV]';
          case default
              stop 'Error: J_multipole > 3 not implemented!';
      end select

      call print_header( tape_strength );

      write(tape_strength,'(/,a,21x,a,/)')'omega[MeV/hbar]',
     &                                    'S(f,omega)';

      if( i_calculation_type .eq. 0 ) then
          write(6,'(/,17x,a,/)') 'Free response';
          write(6,'(a,21x,a,/)') 'omega[MeV/hbar]' , 'S(f,omega)';
      endif

      call flush(tape_strength);
      call flush(6);






c-----Free response, energy sweep
      if( i_calculation_type .eq. 0 ) then

          ! Free response is defined as a response
          ! which does not take into account induced
          ! self-consistent Hamiltonian (H02 = H20 = 0).

          omega = omega_start;
          do while( omega .le. omega_end + 1.D-4 )

              call fam_xy      ( .false. );
              call fam_spurious( .false. );
              call fam_drho    ( .false. );

              Sn = fam_strength( .false. , 1 );
              Sp = fam_strength( .false. , 2 );

              write(            6,'(1f15.5,1f31.10)') omega , Sp+Sn;
              write(tape_strength,'(1f15.5,1f31.10)') omega , Sp+Sn;
              call flush(            6);
              call flush(tape_strength);

              omega = omega + delta_omega;
          enddo

      endif






c-----Fully self-consistent response, energy sweep
      if( i_calculation_type .eq. 1 ) then

          ! Initial guess for self-consistent solution
          ! vector for initial sweep energy is zero.
          ! Otherwise, we use self-consistent solution
          ! of previous energy as initial guess
          ! for the following energy.
          dh_1      = COMPLEX( 0.D0 , 0.D0 );
          dh_2      = COMPLEX( 0.D0 , 0.D0 );
          dDelta_pl = COMPLEX( 0.D0 , 0.D0 );
          dDelta_mi = COMPLEX( 0.D0 , 0.D0 );

          ! One can easily parallelize the following energy
          ! loop for systematic, large scale calculations
          omega = omega_start;
          do while( omega .le. omega_end + 1.D-4 )

              call iter_fam( .false. );

              Sn = fam_strength( .false. , 1 );
              Sp = fam_strength( .false. , 2 );

              write(tape_strength,'(1f15.5,1f31.10)') omega , Sp+Sn;
              call flush(tape_strength);

              omega = omega + delta_omega;
          enddo

      endif






c-----Calculation for given energy, printing vector density
      if( i_calculation_type .eq. 2 ) then

          dh_1      = COMPLEX( 0.D0 , 0.D0 );
          dh_2      = COMPLEX( 0.D0 , 0.D0 );
          dDelta_pl = COMPLEX( 0.D0 , 0.D0 );
          dDelta_mi = COMPLEX( 0.D0 , 0.D0 );

          omega = omega_print;
          call iter_fam( .false. );

          Sn = fam_strength( .false. , 1 );
          Sp = fam_strength( .false. , 2 );

          write(tape_strength,'(1f15.5,1f31.10)') omega , Sp+Sn;
          call flush(tape_strength);

          call print_dens( .false. );

          ! Calculates and prints the localization function
          !call       nuclocfunc( .false. );
          !call print_nuclocfunc( .false. );

      endif






      close(tape_strength);






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END start_fam() *****************************';
      write(6,*) '';
      endif

      return;
      end;
