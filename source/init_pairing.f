c======================================================================c

      subroutine init_pairing( lpr )

c======================================================================c

      USE dirqfampar;
      USE gfvmod;
      USE fam;
      USE simplex;
      USE pairparams;
      USE Wpairing;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;

      common /baspar/ hom, hb0, b0;
      common /defbas/ beta0, q, bp, bz;



      DOUBLE PRECISION Vz( 0: N0FX , 0:N0FX , 0:2*N0FX );
      DOUBLE PRECISION Vr( 0:(N0FX/2) , -N0FX:N0FX ,
     &                     0:(N0FX/2) , -N0FX:N0FX ,
     &                     0:(2*(N0FX/2)+N0FX)         );

      CHARACTER fg1, fg2;
      pi = 3.14159265358979324D0;
      NZMAX = N0FX;
      NRMAX = N0FX/2;
      MLMAX = N0FX;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN init_pairing() ************************';
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

                  nz      = nz1+nz2-NZZ;
                  nz_half = nz/2;

                  call assert( nz     .le.IGFV , 'IGFV too small' );
                  call assert( nz_half.le.IGFV , 'IGFV too small' );

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

                      !See Appendix (E.4) and (E.5)
                      if( abs(ml1-ml2) .ne. K_multipole ) CYCLE;


                      lmt = nr1+nr2+(abs(ml1)+abs(ml2)-abs(ml1-ml2))/2;
                      do NRR = 0 , lmt

                          nr = nr1 + nr2 - NRR +
     &                        (abs(ml1)+abs(ml2)-abs(ml1-ml2))/2;

                          call assert( nr.le.IGFV , 'IGFV too small' );

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
                          x = x / ( a2 + b0*b0*bp*bp )**DBLE(nr+1);

                          call assert(NRR.le.2*NRMAX+MLMAX,'OutOfBnd');
                          Vr( nr1 , ml1 , nr2 , ml2 , NRR ) = x;

                      enddo

                  enddo
              enddo
          enddo
      enddo






      fac = DSQRT(b0*b0*b0*bp*bp*bz);
      fac = fac / ( (2*pi)**0.75D0 );

      il = 0;
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

                  do Nr = 0 , nr1+nr2+(abs(ml1)+abs(ml2)-abs(ml1-ml2))/2
                     do Nz = 0 , nz1+nz2
                        if( mod(Nz,2) .ne. mod(nz1+nz2,2) ) CYCLE;

                        call assert( nz1.le.NZMAX        , 'OutOfBnd' );
                        call assert( nz2.le.NZMAX        , 'OutOfBnd' );
                        call assert( nr1.le.NRMAX        , 'OutOfBnd' );
                        call assert( nr2.le.NRMAX        , 'OutOfBnd' );
                        call assert( abs(ml1).le.MLMAX   , 'OutOfBnd' );
                        call assert( abs(ml2).le.MLMAX   , 'OutOfBnd' );
                        call assert( Nz.le.2*NZMAX       , 'OutOfBnd' );
                        call assert( Nr.le.2*NRMAX+MLMAX , 'OutOfBnd' );

                        il = il + 1;
                        W(il) = fac * Vz( nz1 , nz2 , Nz )
     &                              * Vr( nr1 , ml1 , nr2 , ml2 , Nr );

                     enddo
                  enddo

               enddo
            enddo

         enddo
      enddo
      call assert( il.eq.NWMAX , 'il =/= NWMAX' );






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END init_pairing() **************************';
      write(6,*) '';
      endif

      return;
      end;
