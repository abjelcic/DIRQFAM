c======================================================================c

      subroutine fam_broyden( lpr , initialize )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE simplex;
      USE fam_iter;
      USE dh;
      USE dDelta;
      USE fambroydenmod;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;
      LOGICAL initialize;



      CHARACTER fg1, fg2;
      DATA bmix /0.3D0/;
      DATA xmi  /0.1D0/;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_broyden() *************************';
      write(6,*) '';
      endif






c-----Initialization
      if( initialize .eqv. .true. ) then

         vin = 0.D0;
         il  = 0;
         do it = 1 , 2
            do ib2 = 1 , N_blocks
               do ib1 = 1 , ib2

                  j0 = ia_spx(ib2);
                  j1 = ia_spx(ib2)+id_spx(ib2)-1;
                  do j = j0 , j1
                     fg2 = fg_spx(j-ia_spx(ib2)+1,ib2);
                     ml2 = ml_spx(j-ia_spx(ib2)+1,ib2);

                     i0 = ia_spx(ib1);
                     i1 = ia_spx(ib1)+id_spx(ib1)-1;
                     if( ib1 .eq. ib2 ) i1 = j;
                     do i = i0 , i1
                        fg1 = fg_spx(i-ia_spx(ib1)+1,ib1);
                        ml1 = ml_spx(i-ia_spx(ib1)+1,ib1);


                        if( abs(ml1-ml2).eq.K_multipole ) then

                           il = il + 1;
                           vin(il) = DREAL(dh_1(i,j,it));

                           il = il + 1;
                           vin(il) = DIMAG(dh_1(i,j,it));

                        elseif( abs(ml1+ml2+1).eq.K_multipole .and.
     &                          fg1.ne.fg2                ) then

                           il = il + 1;
                           vin(il) = DREAL(dh_1(i,j,it));

                           il = il + 1;
                           vin(il) = DIMAG(dh_1(i,j,it));

                        endif


                        if( fg1.eq.'f' .and. fg2.eq.'f' .and.
     &                      abs(ml1-ml2).eq.K_multipole  ) then

                            il = il + 1;
                            vin(il) = DREAL(dDelta_pl(i,j,it));

                            il = il + 1;
                            vin(il) = DIMAG(dDelta_pl(i,j,it));

                            il = il + 1;
                            vin(il) = DREAL(dDelta_mi(i,j,it));

                            il = il + 1;
                            vin(il) = DIMAG(dDelta_mi(i,j,it));

                        endif


                     enddo

                  enddo

               enddo
            enddo
         enddo

         call assert( il.eq.nn , 'il =/= nn in fam_broyden()' );

         return;

      endif






c-----Inserting the new vector in Broyden's mixing procedure
      if( initialize .eqv. .false. ) then

         il = 0;
         do it = 1 , 2
            do ib2 = 1 , N_blocks
               do ib1 = 1 , ib2

                  j0 = ia_spx(ib2);
                  j1 = ia_spx(ib2)+id_spx(ib2)-1;
                  do j = j0 , j1
                     fg2 = fg_spx(j-ia_spx(ib2)+1,ib2);
                     ml2 = ml_spx(j-ia_spx(ib2)+1,ib2);

                     i0 = ia_spx(ib1);
                     i1 = ia_spx(ib1)+id_spx(ib1)-1;
                     if( ib1 .eq. ib2 ) i1 = j;
                     do i = i0 , i1
                        fg1 = fg_spx(i-ia_spx(ib1)+1,ib1);
                        ml1 = ml_spx(i-ia_spx(ib1)+1,ib1);


                        if( abs(ml1-ml2).eq.K_multipole ) then

                           il = il + 1;
                           vou(il) = DREAL(dh_1(i,j,it)) - vin(il);

                           il = il + 1;
                           vou(il) = DIMAG(dh_1(i,j,it)) - vin(il);

                        elseif( abs(ml1+ml2+1).eq.K_multipole .and.
     &                          fg1.ne.fg2                ) then

                           il = il + 1;
                           vou(il) = DREAL(dh_1(i,j,it)) - vin(il);

                           il = il + 1;
                           vou(il) = DIMAG(dh_1(i,j,it)) - vin(il);

                        endif


                        if( fg1.eq.'f' .and. fg2.eq.'f' .and.
     &                      abs(ml1-ml2).eq.K_multipole  ) then

                           il = il + 1;
                           vou(il) = DREAL(dDelta_pl(i,j,it)) - vin(il);

                           il = il + 1;
                           vou(il) = DIMAG(dDelta_pl(i,j,it)) - vin(il);

                           il = il + 1;
                           vou(il) = DREAL(dDelta_mi(i,j,it)) - vin(il);

                           il = il + 1;
                           vou(il) = DIMAG(dDelta_mi(i,j,it)) - vin(il);

                        endif


                     enddo

                  enddo

               enddo
            enddo
         enddo

         call assert( il.eq.nn , 'il =/= nn in fam_broyden()' );

      endif






c-----Calculation of the relative difference (Frobenious norm)
c-----between two consecutive Broyden's vectors
      s1 = 0.D0;
      s2 = 0.D0;
      do i = 1 , nn
         s1 = s1 + vin(i)**2.D0;
         s2 = s2 + vou(i)**2.D0;
      enddo
      if( s1 .gt. 1.D-10 ) then
          error = DSQRT(s2/s1);
      else
          error = 1.D+99;
      endif






c-----Broyden's mixing procedure starts here
      if ( mm.eq.0 .or. iter.eq.1 ) then
         do i = 1 , nn
            vin(i) = vin(i) + xmi*vou(i);
         enddo
         ilast = 0;
      else
         iuse = min( iter-1 - 1, mm );
         ipos1 = iter - 2 - ( (iter-3)/mm )*mm;
         inex = iter - 1 - ( (iter-2)/mm )*mm;

         if( iter .eq. 2 ) then
            do j = 1 , mm
               do i = 1 , nn
                  df(i,j) = 0.D0;
                  dv(i,j) = 0.D0;
               enddo
            enddo
            bw0 = 0.01D0;
         else
            do i = 1 , nn
               df(i,ipos1) = vou(i) - df(i,ipos1);
               dv(i,ipos1) = vin(i) - dv(i,ipos1);
            enddo

            dnorm = sqrt( dnrm2(nn,df(1,ipos1),1)**2.0D0 );
            call dscal( nn, 1.0D0 / dnorm, df(1,ipos1), 1 );
            call dscal( nn, 1.0D0 / dnorm, dv(1,ipos1), 1 );
         endif
         do i = 1 , iuse
            do j = i+1 , iuse
               bbeta(i,j) = ddot( nn, df(1,j), 1, df(1,i), 1 );
            enddo
            bbeta(i,i) = 1.0D0 + bw0*bw0;
         enddo

         call dsytrf( 'U', iuse, bbeta, mm, ibwork, bwork, mm, info );
         if( info .ne. 0 ) stop '   in broyden: info at DSYTRF  ';

         call dsytri( 'U', iuse, bbeta, mm, ibwork, bwork, info );
         if( info .ne. 0 ) stop '   in broyden: info at DSYTRI  ';

         do i = 1 , iuse
            do j = i+1 , iuse
               bbeta(j,i) = bbeta(i,j);
            enddo
            bwork(i) = ddot( nn, df(1,i), 1, vou, 1 );
         enddo

         do i = 1 , nn
            curv(i) = bmix * vou(i);
         enddo
         do i = 1 , iuse
            gamma = 0.0D0;
            do j = 1 , iuse
               gamma = gamma + bbeta(j,i)*bwork(j);
            enddo
            do k = 1 , nn
               curv(k) = curv(k) - gamma*( dv(k,i) + bmix * df(k,i) );
            enddo
         enddo

         call dcopy( nn, vou, 1, df(1,inex), 1 );
         call dcopy( nn, vin, 1, dv(1,inex), 1 );

         curvature = ddot( nn, vou, 1, curv, 1 );

         ilast = 1
         do i = 1 , nn
            vin(i) = vin(i) + curv(i);
         enddo

      endif
c-----Broyden's mixing procedure ends here






c-----Extracting the new mixed vector from Broyden's mixing procedure
      il = 1;
      do it = 1 , 2
         do ib2 = 1 , N_blocks
            do ib1 = 1 , ib2

               j0 = ia_spx(ib2);
               j1 = ia_spx(ib2)+id_spx(ib2)-1;
               do j = j0 , j1
                  fg2 = fg_spx(j-ia_spx(ib2)+1,ib2);
                  ml2 = ml_spx(j-ia_spx(ib2)+1,ib2);

                  i0 = ia_spx(ib1);
                  i1 = ia_spx(ib1)+id_spx(ib1)-1;
                  if( ib1 .eq. ib2 ) i1 = j;
                  do i = i0 , i1
                     fg1 = fg_spx(i-ia_spx(ib1)+1,ib1);
                     ml1 = ml_spx(i-ia_spx(ib1)+1,ib1);


                     if( abs(ml1-ml2).eq.K_multipole ) then

                        dh_1(i,j,it) = COMPLEX(vin(il),vin(il+1));
                        il = il + 2;

                     elseif( abs(ml1+ml2+1).eq.K_multipole .and.
     &                       fg1.ne.fg2                ) then

                        dh_1(i,j,it) = COMPLEX(vin(il),vin(il+1));
                        il = il + 2;

                     endif


                     if( fg1.eq.'f' .and. fg2.eq.'f' .and.
     &                   abs(ml1-ml2).eq.K_multipole  ) then

                        dDelta_pl(i,j,it) = COMPLEX(vin(il),vin(il+1));
                        il = il + 2;

                        dDelta_mi(i,j,it) = COMPLEX(vin(il),vin(il+1));
                        il = il + 2;

                     endif


                  enddo

               enddo

            enddo
         enddo
      enddo






c-----Full construction of matrix dh_1
      do it = 1 , 2
         do ib2 = 1 , N_blocks
            do ib1 = ib2 , N_blocks

               j0 = ia_spx(ib2);
               j1 = ia_spx(ib2)+id_spx(ib2)-1;
               do j = j0 , j1
                  fg2 = fg_spx(j-ia_spx(ib2)+1,ib2);
                  ml2 = ml_spx(j-ia_spx(ib2)+1,ib2);

                  i0 = ia_spx(ib1);
                  i1 = ia_spx(ib1)+id_spx(ib1)-1;
                  if( ib1 .eq. ib2 ) i0 = j+1;
                  do i = i0 , i1
                     fg1 = fg_spx(i-ia_spx(ib1)+1,ib1);
                     ml1 = ml_spx(i-ia_spx(ib1)+1,ib1);

                     if( abs(ml1-ml2).eq.K_multipole ) then
                         dh_1(i,j,it) = + dh_1(j,i,it);
                     elseif( abs(ml1+ml2+1).eq.K_multipole .and.
     &                       fg1.ne.fg2                ) then
                         dh_1(i,j,it) = - dh_1(j,i,it);
                     endif

                  enddo

               enddo

            enddo
         enddo
      enddo






c-----Construction of matrix dh_2
      do it = 1 , 2
         do ib2 = 1 , N_blocks
           do ib1 = 1 , N_blocks

              j0 = ia_spx(ib2);
              j1 = j0+id_spx(ib2)-1;
              do j = j0 , j1
                 fg2 = fg_spx(j-j0+1,ib2);
                 ml2 = ml_spx(j-j0+1,ib2);

                 i0 = ia_spx(ib1);
                 i1 = i0+id_spx(ib1)-1;
                 do i = i0 , i1
                    fg1 = fg_spx(i-i0+1,ib1);
                    ml1 = ml_spx(i-i0+1,ib1);

                    if( abs(ml1-ml2) .eq. K_multipole ) then
                        if( fg1 .eq. fg2 ) then
                            dh_2(i,j,it) = + dh_1(i,j,it);
                        else
                            dh_2(i,j,it) = - dh_1(j,i,it);
                        endif
                    elseif( abs(ml1+ml2+1) .eq. K_multipole ) then
                        if( fg1 .ne. fg2 ) then
                            dh_2(i,j,it) = - dh_1(j,i,it);
                        endif
                    endif

              enddo

            enddo

          enddo
        enddo
      enddo






c-----Full construction of matrix dDelta_pl
      do it = 1 , 2
         do ib2 = 1 , N_blocks
            do ib1 = ib2 , N_blocks

               j0 = ia_spx(ib2);
               j1 = ia_spx(ib2)+id_spx(ib2)-1;
               do j = j0 , j1
                  fg2 = fg_spx(j-ia_spx(ib2)+1,ib2);
                  ml2 = ml_spx(j-ia_spx(ib2)+1,ib2);

                  i0 = ia_spx(ib1);
                  i1 = ia_spx(ib1)+id_spx(ib1)-1;
                  if( ib1 .eq. ib2 ) i0 = j+1;
                  do i = i0 , i1
                     fg1 = fg_spx(i-ia_spx(ib1)+1,ib1);
                     ml1 = ml_spx(i-ia_spx(ib1)+1,ib1);

                     if( fg1.eq.'f' .and. fg2.eq.'f' .and.
     &                   abs(ml1-ml2).eq.K_multipole  ) then
                        dDelta_pl(i,j,it) = dDelta_pl(j,i,it);
                     endif

                  enddo

               enddo

            enddo
         enddo
      enddo






c-----Full construction of matrix dDelta_mi
      do it = 1 , 2
         do ib2 = 1 , N_blocks
            do ib1 = ib2 , N_blocks

               j0 = ia_spx(ib2);
               j1 = ia_spx(ib2)+id_spx(ib2)-1;
               do j = j0 , j1
                  fg2 = fg_spx(j-ia_spx(ib2)+1,ib2);
                  ml2 = ml_spx(j-ia_spx(ib2)+1,ib2);

                  i0 = ia_spx(ib1);
                  i1 = ia_spx(ib1)+id_spx(ib1)-1;
                  if( ib1 .eq. ib2 ) i0 = j+1;
                  do i = i0 , i1
                     fg1 = fg_spx(i-ia_spx(ib1)+1,ib1);
                     ml1 = ml_spx(i-ia_spx(ib1)+1,ib1);

                     if( fg1.eq.'f' .and. fg2.eq.'f' .and.
     &                   abs(ml1-ml2).eq.K_multipole  ) then
                        dDelta_mi(i,j,it) = dDelta_mi(j,i,it);
                     endif

                  enddo

               enddo

            enddo
         enddo
      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_broyden() ***************************';
      write(6,*) '';
      endif

      return;
      end;
