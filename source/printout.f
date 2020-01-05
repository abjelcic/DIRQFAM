c======================================================================c

      subroutine print_dens( lpr )

c======================================================================c

      IMPLICIT REAL*8    (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n)
      include 'dirqfam.par'
      LOGICAL lpr;

      INTEGER*4 tape_strength, tape_rhov;
      common /out_tapes/ tape_strength, tape_rhov;

      common /quadrature/ zb_fam( 1:NGH ), wz( 1:NGH ),
     &                    rb_fam( 1:NGL ), wr( 1:NGL ),
     &                    wzwr( 1:NGH , 1:NGL );

      common /gs_dens/ rhov_GS( -NGH:NGH , 1:NGL , 2 ),
     &                 rhos_GS( -NGH:NGH , 1:NGL     );

      common /fam/ omega_start, omega_end, delta_omega, omega_print,
     &             omega, gamma_smear,
     &             i_calculation_type, i_coulomb, i_pairing,
     &             J_multipole, K_multipole, ISO;

      COMPLEX*16 drho_v, drho_s;
      common /ind_dens/ drho_v( -NGH:NGH , 1:NGL , 2 ),
     &                  drho_s( -NGH:NGH , 1:NGL     );



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN print_dens() **************************';
      write(6,*) '';
      endif






      open( tape_rhov , file   = './output/QFAM_output/rhov.out' ,
     &                  status = 'unknown'                        );

      write(tape_rhov,'(a,a)')'rho_v(r,z,phi,t) = rho0_v(r,z) + 2*eta*',
     &                        'Re[ exp(-i*omega*t) * drho_v(r,z,phi) ]';

      write(tape_rhov,'(a,a,i1,a)') 'drho_v(r,z,phi) = drho_v(r,z)',
     &                              ' * cos(',K_multipole,'*phi)  ';

      call print_header( tape_rhov );

      write(tape_rhov,'(a,1f7.3,a)') 'omega = ', omega, ' [MeV/hbar]';




      write(tape_rhov,'(/,/,4x,a,6x,a,7x,a,11x,a,/)')
     & 'r[fm]' , 'z[fm]' , 'rho0_v[fm^-3]' , 'drho_v(r,z)[fm^-3]';

      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

            write(tape_rhov,'(f10.5,f12.5,E18.7E3,E17.7E3,a,E16.7E3,a)')
     &      rb_fam(il), DBLE(isign(1,ih))*zb_fam(abs(ih)),
     &      rhov_GS(ih,il,1)+rhov_GS(ih,il,2),
     &      DREAL( drho_v(ih,il,1)+drho_v(ih,il,2) ),
     &      '  +',
     &      DIMAG( drho_v(ih,il,1)+drho_v(ih,il,2) ),
     &      ' i';
          enddo
      enddo

      close(tape_rhov);






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END print_dens() ****************************';
      write(6,*) '';
      endif

      return;
      end;






c======================================================================c

      subroutine print_header( tape )

c======================================================================c

      IMPLICIT REAL*8    (a-h,o-z)
      IMPLICIT INTEGER*4 (i-n)
      include 'dirqfam.par'
      INTEGER*4 tape;

      CHARACTER parname*10;
      common /partyp/ parname;

      common /basnnn/ n0f, n0b;

      CHARACTER*2 nucnam;
      common /nucnuc/ amas, nneu, npro, nmas, nucnam;

      common /defbas/ beta0, q, bp, bz;

      common /fam/ omega_start, omega_end, delta_omega, omega_print,
     &             omega, gamma_smear,
     &             i_calculation_type, i_coulomb, i_pairing,
     &             J_multipole, K_multipole, ISO;

      common /fam_iter/ error, tol, iter, iter_max;






      if( i_calculation_type .eq. 0 ) then
          write(tape,'(a)') 'Free response';
      else
          write(tape,'(2a,E7.1)') 'Fully self-consistent response ',
     &                            'with tolerance = ', tol;
      endif

      write(tape,'(i3,a,a,i2,a,a,f7.3,a,a)') nmas, nucnam, ', ',
     &                                       n0f, ' shells, ',
     &                                       'beta0 = ', beta0, ', ',
     &                                       parname;

      if( i_coulomb .eq. 1 ) then
          write(tape,'(a)')'Electromagnetism included in calculation';
      else
          write(tape,'(a)')'Electromagnetism excluded from calculation';
      endif

      if( i_pairing .eq. 1 ) then
          write(tape,'(a)')'Pairing included in calculation';
      else
          write(tape,'(a)')'Pairing excluded from calculation';
      endif

      if( ISO .eq. 0 ) then
          write(tape,'(a)')'Isoscalar excitation';
      else
          write(tape,'(a)')'Isovector excitation';
      endif

      write(tape,'(a,i1)') 'J = ', J_multipole;
      write(tape,'(a,i1)') 'K = ', K_multipole;
      write(tape,'(a,1f7.3,a)') 'gamma_smear = ', gamma_smear, ' [MeV]';
      write(tape,'(a,i2,a,i2)') 'NGH = ', NGH, ', NGL = ', NGL;






      return;
      end;
