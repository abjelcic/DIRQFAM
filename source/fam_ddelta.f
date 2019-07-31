c======================================================================c

      subroutine fam_ddelta( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /fam/ omega_start, omega_end, delta_omega, omega_print,
     &             omega, gamma_smear,
     &             i_calculation_type, i_coulomb, i_pairing,
     &             J_multipole, K_multipole, ISO;

      CHARACTER fg_spx;
      common /simplex/ N_total         , N_blocks        ,
     &                 ia_spx(NBX)     , id_spx(NBX)     ,
     &                 nf_size(NBX)    , ng_size(NBX)    ,
     &                 nz_spx(NBSX,NBX), nr_spx(NBSX,NBX),
     &                 ml_spx(NBSX,NBX), fg_spx(NBSX,NBX);

      common /TMR_param/ G_pairing, a_pairing;

      LOGICAL dh_nnz, dDelta_nnz, dkappa_nnz, f_nnz;
      common /nnz_blocks/ dh_nnz    ( NBX , NBX ),
     &                    dDelta_nnz( NBX , NBX ),
     &                    dkappa_nnz( NBX , NBX ),
     &                    f_nnz     ( NBX , NBX );

      common /W_pairing/ W( NWMAX );

      COMPLEX*16 dDelta_pl, dDelta_mi;
      common /dDelta/ dDelta_pl( NTX , NTX , 2 ),
     &                dDelta_mi( NTX , NTX , 2 );

      COMPLEX*16 dkappa_pl, dkappa_mi;
      common /delta_kappa/ dkappa_pl( NTX , NTX , 2 ),
     &                     dkappa_mi( NTX , NTX , 2 );



      COMPLEX*16 P_pl( 0:2*N0FX , 0:2*N0FX , 2 );
      COMPLEX*16 P_mi( 0:2*N0FX , 0:2*N0FX , 2 );

      CHARACTER fg1, fg2;
      COMPLEX*16 z, z1, z2;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_ddelta() ***************************';
      write(6,*) '';
      endif

      if( i_pairing .eq. 0 ) then
          return;
      endif






c-----Selection rules test
      if( .false. ) then

      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dkappa_nnz(ib1,ib2) ) then

                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);

                  do j = j0 , j0-1+id_spx(ib2)
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);
                      do i = i0 , i0-1+id_spx(ib1)
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          ! We work in only Delta_ff approximation
                          if( fg1.ne.'f' .or. fg2.ne.'f') CYCLE;


                          z = dkappa_pl(i,j,it) + dkappa_pl(j,i,it);
                          if( ABS(z) .gt. 1.D-10 ) then
                              if( abs(ml1-ml2) .ne. K_multipole ) then
                                  write(6,*)'Should be zero: ', z;
                                  stop 'Error: |ml1-ml2| =/= K!';
                              endif
                          endif
                          z = dkappa_mi(i,j,it) + dkappa_mi(j,i,it);
                          if( ABS(z) .gt. 1.D-10 ) then
                              if( abs(ml1-ml2) .ne. K_multipole ) then
                                  write(6,*)'Should be zero: ', z;
                                  stop 'Error: |ml1-ml2| =/= K!';
                              endif
                          endif


                      enddo
                  enddo


                  endif
              enddo
          enddo
      enddo

      endif






c-----Symmetrization of dkappa
      do it = 1 , 2
          do i = 1 , N_total
              do j = 1 , i-1
                  dkappa_pl(i,j,it) = + dkappa_pl(i,j,it)
     &                                + dkappa_pl(j,i,it);
                  dkappa_mi(i,j,it) = + dkappa_mi(i,j,it)
     &                                + dkappa_mi(j,i,it);
              enddo
          enddo
      enddo






c-----Calculation of P_pl and P_mi
      P_pl = COMPLEX( 0.D0 , 0.D0 );
      P_mi = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2

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


                     ii = ia_spx(ib1) + i - 1;
                     jj = ia_spx(ib2) + j - 1;
                     z1 = dkappa_pl( ii , jj , it );
                     z2 = dkappa_mi( ii , jj , it );

                     do Nr = 0 , nr1+nr2+(abs(ml1)+abs(ml2)-K)/2
                        do Nz = 0 , nz1+nz2
                           if( mod(Nz,2) .ne. mod(nz1+nz2,2) ) CYCLE;

                           il = il + 1;

                           P_pl(Nz,Nr,it) = P_pl(Nz,Nr,it) + W(il)*z1;

                           P_mi(Nz,Nr,it) = P_mi(Nz,Nr,it) + W(il)*z2;

                        enddo
                     enddo

                  enddo
               enddo

            enddo
         enddo

      enddo






c-----Calculation of dDelta_pl and dDelta_mi
      dDelta_pl = COMPLEX( 0.D0 , 0.D0 );
      dDelta_mi = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2

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

                     ii = ia_spx(ib1) + i - 1;
                     jj = ia_spx(ib2) + j - 1;

                     do Nr = 0 , nr1+nr2+(abs(ml1)+abs(ml2)-K)/2
                        do Nz = 0 , nz1+nz2
                           if( mod(Nz,2) .ne. mod(nz1+nz2,2) ) CYCLE;

                           il = il + 1;

                           dDelta_pl(ii,jj,it) = dDelta_pl(ii,jj,it)
     &                                           + P_pl(Nz,Nr,it)*W(il);

                           dDelta_mi(ii,jj,it) = dDelta_mi(ii,jj,it)
     &                                           + P_mi(Nz,Nr,it)*W(il);

                        enddo
                     enddo

                  enddo
               enddo

            enddo
         enddo

      enddo






c-----Calculation of full dDelta_pl and dDelta_mi
      fac = - G_pairing;
      if( K_multipole .ne. 0 ) then
          fac = 0.5D0 * fac;
      endif
      do it = 1 , 2
          do i = 1 , N_total
              do j = 1 , i
                  dDelta_pl(i,j,it) = fac * dDelta_pl(i,j,it);
                  dDelta_mi(i,j,it) = fac * dDelta_mi(i,j,it);

                  dDelta_pl(j,i,it) = dDelta_pl(i,j,it);
                  dDelta_mi(j,i,it) = dDelta_mi(i,j,it);
              enddo
          enddo
      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_ddelta() *****************************';
      write(6,*) '';
      endif

      return;
      end;
