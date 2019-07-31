c======================================================================c

      subroutine start_fam( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /fam/ omega_start, omega_end, delta_omega, omega_print,
     &             omega, gamma_smear,
     &             i_calculation_type, i_coulomb, i_pairing,
     &             J_multipole, K_multipole, ISO;

      COMPLEX*16 dh_1, dh_2;
      common /delta_h/ dh_1( NTX , NTX , 2 ),
     &                 dh_2( NTX , NTX , 2 );

      COMPLEX*16 dDelta_pl, dDelta_mi;
      common /dDelta/ dDelta_pl( NTX , NTX , 2 ),
     &                dDelta_mi( NTX , NTX , 2 );



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN start_fam() ***************************';
      write(6,*) '';
      endif






c-----Printing strength.out
      open( 100 , file   = './output/QFAM_output/strength.out' ,
     &            status = 'unknown' );

      write(100,'(2a)',advance='no') 'S(f,omega) = -1/pi * ',
     &                               'Im[Tr[hermconj(f)*drho(omega)]] ';

      select case( J_multipole )
          case( 0 )
              write(100,*) '[fm^4/MeV]';
          case( 1 )
              select case( ISO )
                  case( 0 )
                      write(100,*) '[fm^2/MeV]';
                  case( 1 )
                      write(100,*) '[e^2fm^2/MeV]';
                  case default
                      stop 'Error: Wrong ISO!';
              end select
          case( 2 )
              write(100,*) '[fm^4/MeV]';
          case( 3 )
              write(100,*) '[fm^6/MeV]';
          case default
              stop 'Error: J_multipole > 3 not implemented!';
      end select

      call print_header( 100 );

      write(100,*) '';
      write(100,'(a,a)') 'omega[MeV/hbar]        ',
     &                   '             S(f,omega)';
      write(100,*) '';


      if( i_calculation_type .eq. 0 ) then
          write(6,*) '';
          write(6,*) '                ',
     &               'Hartree response';
          write(6,*) '';
          write(6,'(a,a)') 'omega[MeV/hbar]        ',
     &                     '             S(f,omega)';
          write(6,*) '';
      endif






c-----Hartree response, frequency sweep
      if( i_calculation_type .eq. 0 ) then

          ! Hartree response is defined as a response
          ! which doesn't take into account induced
          ! self-consistent Hamiltonian H02,H20.
          ! Therefore, we set dh and dDelta to zero.
          dh_1      = COMPLEX( 0.D0 , 0.D0 );
          dh_2      = COMPLEX( 0.D0 , 0.D0 );
          dDelta_pl = COMPLEX( 0.D0 , 0.D0 );
          dDelta_mi = COMPLEX( 0.D0 , 0.D0 );

          omega = omega_start;
          do while( omega .le. omega_end + 1.D-4 )

              call fam_drhodkappa( .false. );

              Sn = fam_strength( .false. , 1 );
              Sp = fam_strength( .false. , 2 );

              write(  6,'(1f15.5,1f31.10)') omega , Sp+Sn;
              call flush(6);

              write(100,'(1f15.5,1f31.10)') omega , Sp+Sn;
              call flush(100);

              omega = omega + delta_omega;
          enddo

      endif






c-----Fully self-consistent response, frequency sweep
      if( i_calculation_type .eq. 1 ) then

          ! Initial guess of self-consistent Broyden
          ! vector for initial frequency is zero.
          ! Otherwise, we use self-consistent solution
          ! of previous frequency as initial guess
          ! for the following frequency
          dh_1      = COMPLEX( 0.D0 , 0.D0 );
          dh_2      = COMPLEX( 0.D0 , 0.D0 );
          dDelta_pl = COMPLEX( 0.D0 , 0.D0 );
          dDelta_mi = COMPLEX( 0.D0 , 0.D0 );

          ! One can easily parallelize this frequency
          ! loop for systematic, large scale calculations
          omega = omega_start;
          do while( omega .le. omega_end + 1.D-4 )

              call iter_fam( .false. );

              Sn = fam_strength( .false. , 1 );
              Sp = fam_strength( .false. , 2 );

              write(100,'(1f15.5,1f31.10)') omega , Sp+Sn;
              call flush(100);

              omega = omega + delta_omega;
          enddo

      endif






c-----Calculation for given frequency, printing vector density
      if( i_calculation_type .eq. 2 ) then

          omega = omega_print;

          dh_1      = COMPLEX( 0.D0 , 0.D0 );
          dh_2      = COMPLEX( 0.D0 , 0.D0 );
          dDelta_pl = COMPLEX( 0.D0 , 0.D0 );
          dDelta_mi = COMPLEX( 0.D0 , 0.D0 );

          call iter_fam( .false. );

          Sn = fam_strength( .false. , 1 );
          Sp = fam_strength( .false. , 2 );

          write(100,'(1f15.5,1f31.10)') omega , Sp+Sn;
          call flush(100);

          !Printing time dependent vector density
          call print_dens( .false. );

      endif






      close(100);






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END start_fam() *****************************';
      write(6,*) '';
      endif

      return;
      end;
