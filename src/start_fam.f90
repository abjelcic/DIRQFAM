!======================================================================!

      subroutine start_fam()

!======================================================================!
      use fam_input;
      use dh02dh20matrix;
      implicit none;

      double precision , parameter :: pi = 4.d0*atan(1.d0);
      double complex   , parameter :: II = cmplx( 0.d0 , 1.d0 , kind=8 );
      integer                      :: tape_strength;
      double precision             :: omega;
      integer                      :: it;
      double complex   , external  :: fam_strength;
      double complex               :: S;
      double precision             :: dBdw;
      double precision , external  :: iter_fam;
      double precision             :: relResError;
      double complex               :: contour_integral;
      integer                      :: ipoint;
      double precision             :: phi_n;
      double complex               :: omega_gamma;
      double precision             :: w_n;
      double precision             :: Omega_b;
      integer                      :: N_it;




      ! Free response, energy sweep.
      if( calculation_type == 0 ) then

          write(6,'(/,x,a)') 'Calculating free response, i.e. the residual interaction is ignored.';

          open( newunit=tape_strength , file='./output/QFAM_output/strength.out' , action='write' );

          write(tape_strength,'(a,a,a)') 'dB/dw(f,omega) = -1/pi * Im[S(f,omega)].';

          call print_header( tape_strength );

          write(            6,'(/,x,a,x,a,x,a,x,a,x,a,x,a,/)') 'omega'          , [ character(len=15) :: '[MeV]'                                        ] , &
                                                               'dB/dw(f,omega)' , [ character(len=15) :: strengthUnitOfMeasurement(J_multipole,Isospin) ] , &
                                                               'S(f,omega)'     , [ character(len=15) :: strengthUnitOfMeasurement(J_multipole,Isospin) ];
          write(tape_strength,'(/,x,a,x,a,x,a,x,a,x,a,x,a,/)') 'omega'          , [ character(len=15) :: '[MeV]'                                        ] , &
                                                               'dB/dw(f,omega)' , [ character(len=15) :: strengthUnitOfMeasurement(J_multipole,Isospin) ] , &
                                                               'S(f,omega)'     , [ character(len=15) :: strengthUnitOfMeasurement(J_multipole,Isospin) ];

          call flush(tape_strength);
          call flush(            6);


          omega = omega_start;
          do while( omega <= omega_end+1.d-6 )

              do it = 1 , 2
                  call setToZeroBlockMatrix( dh20(it) );
                  call setToZeroBlockMatrix( dh02(it) );
              end do

              call fam_xy( omega + gamma_smear*II );
              call fam_spurious();

              S    = fam_strength();
              dBdw = -1/pi * imag(S);

              write(            6,'(e16.9,6x,e16.9,15x,e16.9,a,e16.9,a)') omega , dBdw , real(S),' + ',imag(S),' i';
              write(tape_strength,'(e16.9,6x,e16.9,15x,e16.9,a,e16.9,a)') omega , dBdw , real(S),' + ',imag(S),' i';
              call flush(            6);
              call flush(tape_strength);

              omega = omega + delta_omega;
          end do


          close(tape_strength);

      end if




      ! Fully self-consistent response, energy sweep.
      if( calculation_type == 1 ) then

          open( newunit=tape_strength , file='./output/QFAM_output/strength.out' , action='write' );

          write(tape_strength,'(a,a,a)') 'dB/dw(f,omega) = -1/pi * Im[S(f,omega)].';

          call print_header( tape_strength );

          write(tape_strength,'(/,x,a,x,a,x,a,x,a,x,a,x,a,x,a,a,/)') 'omega'             , [ character(len=15) :: '[MeV]'                                        ] , &
                                                                     'dB/dw(f,omega)'    , [ character(len=15) :: strengthUnitOfMeasurement(J_multipole,Isospin) ] , &
                                                                     'S(f,omega)'        , [ character(len=15) :: strengthUnitOfMeasurement(J_multipole,Isospin) ] , &
                                                                     repeat(' ',15)      , 'rel. residual error';

          call flush(tape_strength);


          omega = omega_start;
          do while( omega <= omega_end+1.d-6 )

              relResError = iter_fam( omega+gamma_smear*II , NoArnoldiVectors , selfConsistencyTolerance );

              S    = fam_strength();
              dBdw = -1/pi * imag(S);

              write(tape_strength,'(e16.9,6x,e16.9,15x,e16.9,a,e16.9,a,5x,e10.3)',advance='no') omega , dBdw , real(S),' + ',imag(S),' i' , relResError;
              write(tape_strength,'(a)',advance='yes') merge( ' (possible low precision)' , repeat(' ',25) , relResError > selfConsistencyTolerance );
              call flush(tape_strength);

              omega = omega + delta_omega;
          end do


          close(tape_strength);

      end if




      ! Calculation and printing for given energy.
      if( calculation_type == 2 ) then

          open( newunit=tape_strength , file='./output/QFAM_output/strength.out' , action='write' );

          write(tape_strength,'(a,a,a)') 'dB/dw(f,omega) = -1/pi * Im[S(f,omega)].';

          call print_header( tape_strength );

          write(tape_strength,'(/,x,a,x,a,x,a,x,a,x,a,x,a,x,a,a,/)') 'omega'             , [ character(len=15) :: '[MeV]'                                        ] , &
                                                                     'dB/dw(f,omega)'    , [ character(len=15) :: strengthUnitOfMeasurement(J_multipole,Isospin) ] , &
                                                                     'S(f,omega)'        , [ character(len=15) :: strengthUnitOfMeasurement(J_multipole,Isospin) ] , &
                                                                     repeat(' ',15)      , 'rel. residual error';

          call flush(tape_strength);


          relResError = iter_fam( omega_print+gamma_smear*II , NoArnoldiVectors , selfConsistencyTolerance );
          call calculate_dh11();
          call calculate_nuclocfunc();

          S    = fam_strength();
          dBdw = -1/pi * imag(S);

          write(tape_strength,'(e16.9,6x,e16.9,15x,e16.9,a,e16.9,a,5x,e10.3)',advance='no') omega_print , dBdw , real(S),' + ',imag(S),' i' , relResError;
          write(tape_strength,'(a)',advance='yes') merge( ' (possible low precision)' , repeat(' ',25) , relResError > selfConsistencyTolerance );
          call flush(tape_strength);


          call print_basis     (             );
          call print_qpEnergies(             );
          call print_uv        (             );
          call print_densities ( omega_print );
          call print_currents  ( omega_print );
          call print_nuclocfunc( omega_print );
          call print_xy        ( omega_print );
          call print_dh20dh02  ( omega_print );
          call print_dh11      ( omega_print );


          close(tape_strength);

      end if




      ! Contour integration method.
      if( calculation_type == 3 ) then

          open( newunit=tape_strength , file='./output/QFAM_output/strength.out' , action='write' );

          write(tape_strength,'(a,a,a)') 'dB/dw(f,omega) = -1/pi * Im[S(f,omega)].';

          call print_header( tape_strength );

          write(tape_strength,'(/,x,a,x,a,x,a,x,a,x,a,x,a,x,a,a,/)') 'omega           +  gamma*i'   , [ character(len=15) :: '[MeV]'                                        ] , &
                                                                     'dB/dw(f,omega)'               , [ character(len=15) :: strengthUnitOfMeasurement(J_multipole,Isospin) ] , &
                                                                     'S(f,omega)'                   , [ character(len=15) :: strengthUnitOfMeasurement(J_multipole,Isospin) ] , &
                                                                     repeat(' ',15)                 , 'rel. residual error';

          call flush(tape_strength);



          ! We need even number of integration points for Simpson rule.
          NoContourPoints = NoContourPoints + merge( 1 , 0 , mod(NoContourPoints,2)==1 );

          contour_integral = cmplx( 0.d0 , 0.d0 , kind=8 );
          do ipoint = 0 , NoContourPoints

              phi_n = (2*pi*ipoint)/NoContourPoints;

              omega_gamma = omega_center + omega_radius*exp(II*phi_n);
              relResError = iter_fam( omega_gamma , NoArnoldiVectors , selfConsistencyTolerance );

              S    = fam_strength();
              dBdw = -1/pi * imag(S);

              write(tape_strength,'(e16.9,a,e16.9,a,6x,e16.9,15x,e16.9,a,e16.9,a,5x,e10.3)',advance='no') real(omega_gamma),' + ',imag(omega_gamma),' i' , dBdw , real(S),' + ',imag(S),' i' , relResError;
              write(tape_strength,'(a)',advance='yes') merge( ' (possible low precision)' , repeat(' ',25) , relResError > selfConsistencyTolerance );
              call flush(tape_strength);

              ! Simpson rule: integral_[a,b]f(x)dx = h/3 * sum_{i=0}^n w_i*f(x_i)
              ! x_i = a + i*h, i=0,1,2,...,n, and h=(b-a)/n, n even,
              ! where: w_0=1, w_n=1, and others w_i=4 if i is odd, or w_i=2 if i is even.
              w_n = merge( 4.d0 , 2.d0 , mod(ipoint,2)==1                       );
              w_n = merge( 1.d0 , w_n  , ipoint==0 .or. ipoint==NoContourPoints );

              contour_integral = contour_integral + (2*pi/NoContourPoints)/3 * w_n * ( omega_radius/(2*pi) * S * exp(II*phi_n) );

          end do

          write(tape_strength,'(/,a,f14.7,a,f14.7,a)') 'Contour integration of strength function along a circle with center omega0 = ' , omega_center , ', and radius R = ' , omega_radius , '.';
          write(tape_strength,'(a,i0,a)') 'Simpson''s integration rule with ' , NoContourPoints , ' points is used.';
          write(tape_strength,'(a,e16.9,a,e16.9,a)') '1/(2*pi*i) * integral_{C(omega0,R)} S(f,omega) domega = ' , real(contour_integral),' + ',imag(contour_integral),'i.';
          call flush(tape_strength);


          close(tape_strength);

      end if




      ! Kernel Polynomial Method.
      if( calculation_type == 4 ) then
          Omega_b = +4500.d0;
          N_it    = 1000;
          call kpm( Omega_b , N_it );
      end if




      return;

      contains

          pure function strengthUnitOfMeasurement( J , Isospin ) result(ans)
              implicit none;
              integer        , intent(in)  :: J;
              integer        , intent(in)  :: Isospin;
              character(:)   , allocatable :: ans;
              character(1024)              :: tmp;

                  write(tmp,'(a,i0,a)') '[fm^' , 2*J , '/MeV]';

                  if( Isospin==1 .and. J==1 ) &
                      write(tmp,'(a)') '[e^2fm^2/MeV]';

                  if( J == 0 ) &
                      write(tmp,'(a)') '[fm^4/MeV]';

                  ans = trim(adjustl(tmp));

              return;
          end function strengthUnitOfMeasurement

      end subroutine start_fam
