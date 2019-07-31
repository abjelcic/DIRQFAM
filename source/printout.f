c======================================================================c

      subroutine print_dens( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /gaussh/ xh(0:NGH), wh(0:NGH), zb(0:NGH);
      common /gaussl/ xl(0:NGL), wl(0:NGL), sxl(0:NGL), rb(0:NGL);
      common /dens/ ro(MG,4), dro(MG,4);

      common /fam/ omega_start, omega_end, delta_omega, omega_print,
     &             omega, gamma_smear,
     &             i_calculation_type, i_coulomb, i_pairing,
     &             J_multipole, K_multipole, ISO;

      COMPLEX*16 drho_v;
      COMPLEX*16 drho_s;
      COMPLEX*16 ldrho_v;
      COMPLEX*16 ldrho_s;
      common /ind_dens/ drho_v ( -NGH:NGH , 1:NGL , 2 ),
     &                  drho_s ( -NGH:NGH , 1:NGL     ),
     &                  ldrho_v( -NGH:NGH , 1:NGL , 2 ),
     &                  ldrho_s( -NGH:NGH , 1:NGL     );



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN print_dens() ***************************';
      write(6,*) '';
      endif






      open( 300 , file = './output/QFAM_output/rhov.out' ,
     &      status = 'unknown' );

      write(300,'(a,a)') 'rho_v(r,z,phi,t) = rho0_v(r,z) + 2*eta*',
     &                   'Re[ exp(-i*omega*t) * drho_v(r,z,phi) ]';

      write(300,'(a,a,i1,a)')
     &          'drho_v(r,z,phi) = drho_v(r,z)',
     &          ' * cos(',K_multipole,'*phi)  ';

      call print_header( 300 );

      write(300,'(a,1f7.3,a)') 'omega = ', omega,
     &                         ' [MeV/hbar]';




      write(300,*) '';
      write(300,*) '';
      write(300,'(4x,a,6x,a,7x,a,11x,a)')
     &                  'r[fm]',
     &                  'z[fm]',
     &                  'rho0_v[fm^-3]',
     &                  'drho_v(r,z)[fm^-3]';
      write(300,*) '';

      do il = 1 , NGL
          do ih = -NGH , NGH
              if( ih .eq. 0 ) CYCLE;

              write(300,'(f10.5,f12.5,E18.7E3,E17.7E3,a,E16.7E3,a)')
     &                  rb(il), DBLE(isign(1,ih))*zb(abs(ih)),
     &                  ro( 1+abs(ih) + il*(NGH+1) , 2 ),
     &                  DREAL( drho_v(ih,il,1)+drho_v(ih,il,2) ),
     &                  '  +',
     &                  DIMAG( drho_v(ih,il,1)+drho_v(ih,il,2) ),
     &                  ' i';
          enddo
      enddo


      close(300);






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END print_dens() *****************************';
      write(6,*) '';
      endif

      return;
      end;






c======================================================================c

      subroutine print_header( tape )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      INTEGER tape;

      CHARACTER parname*10;
      common /partyp/ parname;

      common /basnnn/ n0f, n0b;

      CHARACTER*2 nucnam;
      common /nucnuc/ amas, nneu, npro, nmas, nucnam;

      common /fam/ omega_start, omega_end, delta_omega, omega_print,
     &             omega, gamma_smear,
     &             i_calculation_type, i_coulomb, i_pairing,
     &             J_multipole, K_multipole, ISO;






      if( i_calculation_type .eq. 0 ) then
          write(tape,'(a)') 'Hartree response';
      else
          write(tape,'(a)') 'Fully self-consistent response';
      endif

      write(tape,'(i3,a,a,i2,a,a)') nmas, nucnam, ', ',
     &                              n0f, ' shells, ',
     &                              parname;

      if( i_coulomb .eq. 1 ) then
          write(tape,'(a)') 'Electromagnetism included in calculation';
      else
          write(tape,'(a)') 'Electromagnetism excluded from calculation';
      endif

      if( i_pairing .eq. 1 ) then
          write(tape,'(a)') 'Pairing included in calculation';
      else
          write(tape,'(a)') 'Pairing excluded from calculation';
      endif

      if( ISO .eq. 0 ) then
          write(tape,'(a)') 'Isoscalar excitation';
      else
          write(tape,'(a)') 'Isovector excitation';
      endif

      write(tape,'(a,i1)') 'J = ', J_multipole;
      write(tape,'(a,i1)') 'K = ', K_multipole;
      write(tape,'(a,1f7.3,a)') 'gamma_smear = ', gamma_smear, ' [MeV]';
      write(tape,'(a,i2,a,i2)') 'NGH = ', NGH, ', NGL = ', NGL;






      return;
      end;
