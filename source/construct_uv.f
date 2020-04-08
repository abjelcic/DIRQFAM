c======================================================================c

      subroutine construct_v( lpr , it )

c======================================================================c

      USE dirqfampar;
      USE simplex;
      USE v_matrix;
      USE v_energy;
      USE blokap;
      USE bloosc;
      USE quaosc;
      USE blodir;
      USE waveuv;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;
      INTEGER it;



      CHARACTER fg;
      DOUBLE COMPLEX z;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN construct_v() *************************';
      endif






      do ib = 1 , nb

          jj = 1; !Current column of matrix v( : , jj , ib , it )
                  !being under construction

          nf  = id(ib,1);
          ng  = id(ib,2);
          i0f = ia(ib,1);
          kap = kb(ib);
          nh  = nf + ng;


c---------Particles (positive E in Fermi sea)
          klp = ka(ib,it);
          do k = 1 , nf

              klp = klp + 1;
              E_fam_v(jj,ib,it) = equ(klp,it);

              do n = nh+1 , nh+nh

                  if( n .le. nh+nf ) then
                      fg = 'f';
                  else
                      fg = 'g';
                  endif


                  nzz = nz( i0f + n-nh );
                  nrr = nr( i0f + n-nh );
                  mll = ml( i0f + n-nh );
                  mx  = 2*( iabs(kap) - mll ) - 1;

                  z = COMPLEX( -fguv(n,klp,it) , 0.D0 );


                  if( mx.eq.+1 .and. fg.eq.'g' ) then
                      z = z * COMPLEX( 0.D0 , 1.D0 );
                  endif
                  if( mx.eq.-1 .and. fg.eq.'f' ) then
                      z = z * COMPLEX( 0.D0 , 1.D0 );
                  endif

                  ii = index_of_vector( fg , nzz , nrr , mx*mll , ib );
                  v(ii,jj,ib,it) = z;

              enddo

              jj = jj + 1;
          enddo


c---------Antiparticles (negative E in Dirac sea)
          kla = ka(ib,it+2);
          do k = 1 , ng

              kla = kla + 1;
              E_fam_v(jj,ib,it) = equ(kla,it+2);

              do n = nh+1 , nh+nh

                  if( n .le. nh+nf ) then
                      fg = 'f';
                  else
                      fg = 'g';
                  endif


                  nzz = nz( i0f + n-nh );
                  nrr = nr( i0f + n-nh );
                  mll = ml( i0f + n-nh );
                  mx  = 2*( iabs(kap) - mll ) - 1;

                  z = COMPLEX( -fguv(n,kla,it+2) , 0.D0 );


                  if( mx.eq.+1 .and. fg.eq.'g' ) then
                      z = z * COMPLEX( 0.D0 , 1.D0 );
                  endif
                  if( mx.eq.-1 .and. fg.eq.'f' ) then
                      z = z * COMPLEX( 0.D0 , 1.D0 );
                  endif

                  ii = index_of_vector( fg , nzz , nrr , mx*mll , ib );
                  v(ii,jj,ib,it) = z;

              enddo

              jj = jj + 1;
          enddo

          call assert( jj-1 .eq. id_spx(ib) , 'v matrix construction' );

      enddo






      if( lpr ) then

          write(6,*) '================================================';
          write(6,*) '=================== MATRIX v ===================';
          write(6,*) '================================================';
          write(6,*) '';
          write(6,*) '';
          write(6,*) '----------------------------------';
          if( it .eq. 1 ) write(6,*) '------------ Neutrons ----------';
          if( it .eq. 2 ) write(6,*) '------------ Protons -----------';
          write(6,*) '----------------------------------';
          write(6,*) '';
          write(6,*) '';


          do ib = 1 , N_blocks

              write(6,'(a,i2)') 'Block number #', ib;
              write(6,*) '--------------------------------------------';

              do j = 1 , id_spx(ib)

                  write(6,111) 'Column ', j, ', ',
     &                         'E = ', E_fam_v(j,ib,it);

                  write(6,*) 'fg  nz  nr  ml';
                  do i = 1 , id_spx(ib)
                      write(6,222) '  ',
     &                             fg_spx(i,ib), nz_spx(i,ib),
     &                             nr_spx(i,ib), ml_spx(i,ib),
     &                             DREAL( v(i,j,ib,it) ),
     &                             '        +',
     &                             DIMAG( v(i,j,ib,it) ), ' i';
                  enddo
                  write(6,*) '';

              enddo

          enddo

  111 format(a,i4,a,a,1f18.10);
  222 format(a,a1,i4,i4,i4,1f20.10,a,1f20.10,a);


      endif






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END construct_v() ***************************';
      endif

      return;
      end;






c======================================================================c

      subroutine construct_u( lpr , it )

c======================================================================c

      USE dirqfampar;
      USE simplex;
      USE u_matrix;
      USE u_energy;
      USE blokap;
      USE bloosc;
      USE quaosc;
      USE blodir;
      USE waveuv;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;
      INTEGER it;



      CHARACTER fg;
      DOUBLE COMPLEX z;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN construct_u() *************************';
      endif






      do ib = 1 , nb

          jj = 1; !Current column of matrix u( : , jj , ib , it )
                  !being under construction

          nf  = id(ib,1);
          ng  = id(ib,2);
          i0f = ia(ib,1);
          kap = kb(ib);
          nh  = nf + ng;


c---------Particles (positive E in Fermi sea)
          klp = ka(ib,it);
          do k = 1 , nf

              klp = klp + 1;
              E_fam_u(jj,ib,it) = equ(klp,it);

              do n = 1 , nh

                  if( n .le. nf ) then
                      fg = 'f';
                  else
                      fg = 'g';
                  endif


                  nzz = nz( i0f + n );
                  nrr = nr( i0f + n );
                  mll = ml( i0f + n );
                  mx  = 2*( iabs(kap) - mll ) - 1;

                  z = COMPLEX( +fguv(n,klp,it) , 0.D0 );


                  if( mx.eq.+1 .and. fg.eq.'g' ) then
                      z = z * COMPLEX( 0.D0 , 1.D0 );
                  endif
                  if( mx.eq.-1 .and. fg.eq.'f' ) then
                      z = z * COMPLEX( 0.D0 , 1.D0 );
                  endif

                  ii = index_of_vector( fg , nzz , nrr , mx*mll , ib );
                  u(ii,jj,ib,it) = z;

              enddo

              jj = jj + 1;
          enddo


c---------Antiparticles (negative E in Dirac sea)
          kla = ka(ib,it+2);
          do k = 1 , ng

              kla = kla + 1;
              E_fam_u(jj,ib,it) = equ(kla,it+2);

              do n = 1 , nh

                  if( n .le. nf ) then
                      fg = 'f';
                  else
                      fg = 'g';
                  endif


                  nzz = nz( i0f + n );
                  nrr = nr( i0f + n );
                  mll = ml( i0f + n );
                  mx  = 2*( iabs(kap) - mll ) - 1;

                  z = COMPLEX( +fguv(n,kla,it+2) , 0.D0 );


                  if( mx.eq.+1 .and. fg.eq.'g' ) then
                      z = z * COMPLEX( 0.D0 , 1.D0 );
                  endif
                  if( mx.eq.-1 .and. fg.eq.'f' ) then
                      z = z * COMPLEX( 0.D0 , 1.D0 );
                  endif

                  ii = index_of_vector( fg , nzz , nrr , mx*mll , ib );
                  u(ii,jj,ib,it) = z;

              enddo

              jj = jj + 1;
          enddo

          call assert( jj-1 .eq. id_spx(ib) , 'u matrix construction' );

      enddo






      if( lpr ) then

          write(6,*) '================================================';
          write(6,*) '=================== MATRIX u ===================';
          write(6,*) '================================================';
          write(6,*) '';
          write(6,*) '';
          write(6,*) '----------------------------------';
          if( it .eq. 1 ) write(6,*) '------------ Neutrons ----------';
          if( it .eq. 2 ) write(6,*) '------------ Protons -----------';
          write(6,*) '----------------------------------';
          write(6,*) '';
          write(6,*) '';


          do ib = 1 , N_blocks

              write(6,'(a,i2)') 'Block number #', ib;
              write(6,*) '--------------------------------------------';

              do j = 1 , id_spx(ib)

                  write(6,111) 'Column ', j, ', ',
     &                         'E = ', E_fam_u(j,ib,it);

                  write(6,*) 'fg  nz  nr  ml';
                  do i = 1 , id_spx(ib)
                      write(6,222) '  ',
     &                             fg_spx(i,ib), nz_spx(i,ib),
     &                             nr_spx(i,ib), ml_spx(i,ib),
     &                             DREAL( u(i,j,ib,it) ),
     &                             '        +',
     &                             DIMAG( u(i,j,ib,it) ), ' i';
                  enddo
                  write(6,*) '';

              enddo

          enddo

  111 format(a,i4,a,a,1f18.10);
  222 format(a,a1,i4,i4,i4,1f20.10,a,1f20.10,a);


      endif






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END construct_u() ***************************';
      endif

      return;
      end;
