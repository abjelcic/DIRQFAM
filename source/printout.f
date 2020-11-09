c======================================================================c

      subroutine print_dens( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE quadrature;
      USE gs_dens;
      USE ddens;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



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

      subroutine print_nuclocfunc( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE quadrature;
      USE nuclearlocfunc;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN print_nuclocfunc() ********************';
      write(6,*) '';
      endif






      do it = 1 , 2

          itape = tape_nuclocfunc;

          if( it .eq. 1 ) then
              open( itape , file='./output/QFAM_output/neutlocfunc.out',
     &                      status='unknown'                          );
          else
              open( itape , file='./output/QFAM_output/protlocfunc.out',
     &                      status='unknown'                          );
          endif

          write(itape,'(a,a)') 'C(r,z,phi,t) = C(r,z) + 2*eta*',
     &                         'Re[ exp(-i*omega*t) * dC(r,z,phi) ]';

          write(itape,'(a,a,i1,a)') 'dC(r,z,phi) = dC(r,z)',
     &                              ' * cos(',K_multipole,'*phi)  ';

          call print_header( itape );

          write(itape,'(a,1f7.3,a)') 'omega = ', omega, ' [MeV/hbar]';




          write(itape,'(/,/,4x,a,6x,a,9x,a,19x,a,/)')
     &     'r[fm]' , 'z[fm]' , 'C0(r,z)' , 'dC(r,z)';

          do il = 1 , NGL
            do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              write(itape,'(f10.5,f12.5,E18.7E3,E17.7E3,a,E16.7E3,a)')
     &        rb_fam(il), DBLE(isign(1,ih))*zb_fam(abs(ih)),
     &        C0(ih,il,it),
     &        DREAL( dC(ih,il,it) ),
     &        '  +',
     &        DIMAG( dC(ih,il,it) ),
     &        ' i';
            enddo
          enddo

          close(itape);

      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END print_nuclocfunc() **********************';
      write(6,*) '';
      endif

      return;
      end;






c======================================================================c

      subroutine print_header( tape )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE fam_iter;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      INTEGER tape;

      CHARACTER*10 parname;
      CHARACTER*2  nucnam;
      common /partyp/ parname;
      common /basnnn/ n0f, n0b;
      common /nucnuc/ amas, nneu, npro, nmas, nucnam;
      common /defbas/ beta0, q, bp, bz;






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

      call flush(tape);






      return;
      end;
