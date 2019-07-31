c======================================================================c

      subroutine init_pairing( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /mathco/ zero, one, two, half, third, pi;
      common /baspar/ hom, hb0, b0;
      common /defbas/ beta0, q, bp, bz;
      common /gfviv / iv ( -IGFV : IGFV ); ! iv(n)  = (-1)^n
      common /gfvfi / fi (   0   : IGFV ); ! fi(n)  = 1/n!
      common /gfvwf / wf (   0   : IGFV ); ! wf(n)  = sqrt(n!)

      common /fam/ omega_start, omega_end, delta_omega, omega_print,
     &             omega, gamma_smear,
     &             i_calculation_type, i_coulomb, i_pairing,
     &             J_multipole, K_multipole, ISO;

      common /TMR_param/ G_pairing, a_pairing;

      CHARACTER fg_spx;
      common /simplex/ N_total         , N_blocks        ,
     &                 ia_spx(NBX)     , id_spx(NBX)     ,
     &                 nf_size(NBX)    , ng_size(NBX)    ,
     &                 nz_spx(NBSX,NBX), nr_spx(NBSX,NBX),
     &                 ml_spx(NBSX,NBX), fg_spx(NBSX,NBX);

      common /W_pairing/ W( NWMAX );



      parameter( NZMAX = N0FX );
      parameter( NRMAX = N0FX/2 );
      parameter( MLMAX = N0FX );

      REAL*8 Vz( 0:NZMAX , 0:NZMAX , 0:2*NZMAX );
      REAL*8 Vr( 0:NRMAX , -MLMAX:MLMAX ,
     &           0:NRMAX , -MLMAX:MLMAX ,
     &           0:(2*NRMAX+MLMAX)       );

      CHARACTER fg1, fg2;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN init_pairing() *************************';
      write(6,*) '';
      endif


      if( i_pairing .eq. 0 ) then
          return;
      endif

      a2 = a_pairing*a_pairing;






      Vz = 0.D0;
      do nz1 = 0 , NZMAX
          do nz2 = 0 , NZMAX
              do NZZ = 0 , nz1+nz2
                  if( mod(nz1+nz2,2) .ne. mod(NZZ,2) ) CYCLE;

                  nz = nz1+nz2-NZZ;
                  nz_half = nz/2;

                  x = TalmiMoshinsky_1d( nz1 , nz2 , NZZ , nz );
                  x = x * wf(nz)*fi(nz_half)/DSQRT(2.D0**DBLE(nz));
                  if( a_pairing .gt. b0*bz ) then
                      x = x * ( + a2 - b0*b0*bz*bz )**DBLE(nz_half);
                  else
                      x = x * ( - a2 + b0*b0*bz*bz )**DBLE(nz_half);
                      x = x * iv(nz_half);
                  endif
                  x = x / DSQRT( ( a2 + b0*b0*bz*bz )**DBLE(nz+1) );

                  Vz( nz1 , nz2 , NZZ ) = x;

              enddo
          enddo
      enddo






      Vr = 0.D0;
      do nr1 = 0 , NRMAX
          do ml1 = -MLMAX , MLMAX
              if( 2*nr1+abs(ml1) .gt. N0FX ) CYCLE;
              do nr2 = 0 , NRMAX
                  do ml2 = -MLMAX , MLMAX
                      if( 2*nr2+abs(ml2) .gt. N0FX ) CYCLE;

                      if( abs(ml1-ml2) .ne. K_multipole ) CYCLE;

                      lmt = nr1+nr2+(abs(ml1)+abs(ml2)-abs(ml1-ml2))/2;
                      do NRR = 0 , lmt

                          nr = nr1 + nr2 - NRR +
     &                        (abs(ml1)+abs(ml2)-abs(ml1-ml2))/2;

                          x = TalmiMoshinsky_2d( nr1 , +ml1     ,
     &                                           nr2 , -ml2     ,
     &                                           NRR ,  ml1-ml2 ,
     &                                           nr  ,  0        );
                          if( a_pairing .lt. b0*bp ) then
                              x = x * ( + b0*b0*bp*bp - a2 )**DBLE(nr);
                          else
                              x = x * ( - b0*b0*bp*bp + a2 )**DBLE(nr);
                              x = x * iv(nr);
                          endif
                          x = x / ( ( a2 + b0*b0*bp*bp )**DBLE(nr+1) );

                          Vr( nr1 , ml1 , nr2 , ml2 , NRR ) = x;

                      enddo

                  enddo
              enddo
          enddo
      enddo






      fac = DSQRT(b0*b0*b0*bp*bp*bz);
      fac = fac / ( (2*pi)**0.75D0 );

      il = 0;
      K  = K_multipole;
      do ib1 = 1 , N_blocks
         i0 = 1;
         i1 = id_spx(ib1);
         do i = i0 , i1
            if( fg_spx(i,ib1) .ne. 'f' ) CYCLE;

            fg1 = fg_spx(i,ib1);
            nz1 = nz_spx(i,ib1);
            nr1 = nr_spx(i,ib1);
            ml1 = ml_spx(i,ib1);

            do ib2 = 1 , ib1
               j0 = 1;
               j1 = id_spx(ib2);
               if( ib2 .eq. ib1 ) j1 = i;
               do j = j0 , j1
                  if( fg_spx(j,ib2) .ne. 'f' ) CYCLE;

                  fg2 = fg_spx(j,ib2);
                  nz2 = nz_spx(j,ib2);
                  nr2 = nr_spx(j,ib2);
                  ml2 = ml_spx(j,ib2);

                  if( abs(ml1-ml2) .ne. K_multipole ) CYCLE;

                  do Nr = 0 , nr1+nr2+(abs(ml1)+abs(ml2)-K)/2
                     do Nz = 0 , nz1+nz2
                        if( mod(Nz,2) .ne. mod(nz1+nz2,2) ) CYCLE;

                        il = il + 1;
                        W(il) = fac * Vz( nz1 , nz2 , Nz )
     &                              * Vr( nr1 , ml1 , nr2 , ml2 , Nr );

                     enddo
                  enddo

               enddo
            enddo

         enddo
      enddo
      if( il .ne. NWMAX ) then
          stop 'Error: il =/= WSIZE !';
      endif






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END init_pairing() ***************************';
      write(6,*) '';
      endif

      return;
      end;
