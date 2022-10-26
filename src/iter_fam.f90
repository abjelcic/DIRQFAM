!======================================================================!

      double precision function iter_fam( omega_gamma , iter_max , tolerance ) result(relResError)

!======================================================================!
      use fam_input;
      use clock_module;
      implicit none;
      double complex   , intent(in) :: omega_gamma;
      integer          , intent(in) :: iter_max;
      double precision , intent(in) :: tolerance;
      integer                       :: iter;
      double precision , external   :: fam_gmres;
      type(clockClass)              :: clock;


      write(6,'(/,a,f7.4,a,f7.4,a)') 'Start QFAM solver for omega = ' , real(omega_gamma),' + ',imag(omega_gamma),'i.';


      ! This costs +1 QFAM iterations.
      relResError = fam_gmres( 0 , omega_gamma );

      do iter = 1 , iter_max

          write(6,'(i3,a,f12.9,a)') iter , '.Iteration, rel. residual error = ' , relResError , '.';

          ! call clock%tic();
          ! This costs +1 QFAM iterations.
              call fam_dh20dh02();
              call fam_xy( omega_gamma );
              call fam_drhodkappa();
              call fam_ddensdcurr();
              call fam_dpotentials();
              call fam_dh();
              call fam_ddelta();
          ! call clock%toc();
          ! write(6,'(a,f12.7,a)') 'Elapsed time per QFAM iteration: ' , clock%getTime() , ' s.';

          relResError = fam_gmres( iter , omega_gamma );

          if( relResError < tolerance ) then
              call fam_gmres_getSolution( iter );
              write(6,'(a,i0,a,/)') 'QFAM solver converged after ' , iter+2 , ' QFAM iterations.';
              exit;
          end if

          if( iter == iter_max ) then
              call fam_gmres_getSolution( iter );
              write(6,'(a,i0,a,/)') 'QFAM solver interrupted after ' , iter+2 , ' QFAM iterations.';
              exit;
          end if

      end do

      ! This costs +1 QFAM iterations.
      call fam_dh20dh02();
      call fam_xy( omega_gamma );
      call fam_spurious(); ! Translational spurious mode removal modifies only x and y.
      call fam_drhodkappa();
      call fam_ddensdcurr();
      call fam_dpotentials();
      call fam_dh();
      call fam_ddelta();
      call fam_dh20dh02();


      return;
      end function iter_fam
