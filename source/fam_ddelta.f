c======================================================================c

      subroutine fam_ddelta( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE simplex;
      USE pairparams;
      USE nnz_blocks;
      USE Wpairing;
      USE dDelta;
      USE dkappa;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      DOUBLE COMPLEX P_pl( 0:2*N0FX , 0:2*N0FX , 2 );
      DOUBLE COMPLEX P_mi( 0:2*N0FX , 0:2*N0FX , 2 );

      CHARACTER fg1, fg2;
      DOUBLE COMPLEX z, z1, z2;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_ddelta() **************************';
      write(6,*) '';
      endif

      if( i_pairing .eq. 0 ) then
          dDelta_pl = COMPLEX( 0.D0 , 0.D0 );
          dDelta_mi = COMPLEX( 0.D0 , 0.D0 );
          return;
      endif






#ifdef DEBUG
c-----Selection rules test
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dkappa_nnz(ib1,ib2) .eqv. .false. ) CYCLE;

                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);

                  do j = j0 , j0-1+id_spx(ib2)
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);
                      do i = i0 , i0-1+id_spx(ib1)
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          ! We work in Delta_ff only approximation
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


              enddo
          enddo
      enddo
#endif






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


                     ii = i-1+ia_spx(ib1);
                     jj = j-1+ia_spx(ib2);
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

                     ii = i-1+ia_spx(ib1);
                     jj = j-1+ia_spx(ib2);

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
      if( K_multipole .ne. 0 ) fac = 0.5D0 * fac;

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
      write(6,*) '****** END fam_ddelta() ****************************';
      write(6,*) '';
      endif

      return;
      end;
