c======================================================================c

      subroutine fam_ddensdcurr( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE simplex;
      USE PHI;
      USE nnz_blocks;
      USE drho;
      USE ddens;
      USE dcurr;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      DOUBLE PRECISION Ar( NTX , NTX );
      DOUBLE PRECISION Ai( NTX , NTX );

      DOUBLE PRECISION AUX1(    NTX , KTRUNC );
      DOUBLE PRECISION AUX2( KTRUNC , KTRUNC );

      DOUBLE COMPLEX D( NCOORD );
      DOUBLE COMPLEX z;
      CHARACTER fg1, fg2;
      pi = 3.14159265358979324D0;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_ddensdcurr() **********************';
      write(6,*) '';
      endif

c-----Note that drho_1 and drho_2, have block structure given by dh_nnz






#ifdef DEBUG
c-----Selection rules test
      K = K_multipole;
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh_nnz(ib1,ib2) .eqv. .false. ) CYCLE;

                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);

                  do j = j0 , j0-1+id_spx(ib2)
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);
                      do i = i0 , i0-1+id_spx(ib1)
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1 .eq. fg2 ) then
                              z = drho_1(i,j,it) + drho_2(i,j,it);
                              if( ABS(z) .gt. 1.D-10 ) then
                                  if( abs(ml1-ml2) .ne. K ) then
                                      write(6,*)'Should be zero: ', z;
                                      stop 'Error: |ml1-ml2| =/= K!';
                                  endif
                              endif
                          endif

                          if( fg1 .ne. fg2 ) then
                              z = + drho_1(i,j,it) + drho_2(i,j,it);
                              if( ABS(z) .gt. 1.D-10 ) then
                                  if( abs(ml1+ml2+1) .ne. K ) then
                                      write(6,*)'Should be zero: ', z;
                                      stop 'Error: |ml1+ml2+1| =/= K!';
                                  endif
                              endif
                              z = + drho_1(i,j,it) - drho_2(i,j,it);
                              if( ABS(z) .gt. 1.D-10 ) then
                                  if( abs(ml1-ml2) .ne. K ) then
                                      write(6,*)'Should be zero: ', z;
                                      stop 'Error: |ml1-ml2| =/= K!';
                                  endif
                              endif
                          endif


                      enddo
                  enddo


              enddo
          enddo
      enddo
#endif






c-----Calculation of drho_v
      drho_v = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2

          Ar = 0.D0;
          Ai = 0.D0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , ib2
                  if( dh_nnz(ib1,ib2) .eqv. .false. ) CYCLE;

                  j0 = ia_spx(ib2);
                  j1 = ia_spx(ib2)+id_spx(ib2)-1;
                  do j = j0 , j1
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);

                      i0 = ia_spx(ib1);
                      i1 = ia_spx(ib1)+id_spx(ib1)-1;
                      if( ib1 .eq. ib2 ) i1 = j;
                      do i = i0 , i1
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1.ne.fg2 ) CYCLE;
                          if( abs(ml1-ml2).ne.K_multipole ) CYCLE;

                          z = + drho_1(i,j,it) + drho_2(i,j,it)
     &                        + drho_1(j,i,it) + drho_2(j,i,it);
                          if( i .eq. j ) z = 0.5D0 * z;

                          Ar(i,j) = DREAL(z);
                          Ai(i,j) = DIMAG(z);
                      enddo
                  enddo

              enddo
          enddo


          call diagUVtAUV(  N_total , NCOORD ,   k_PHI , N_blocks ,
     &                           Ar ,    NTX ,      Ai ,      NTX ,
     &                       dh_nnz ,    NBX ,
     &                       ia_spx , id_spx ,
     &                        PHI_U ,    NTX , PHI_SVt ,   KTRUNC , D ,
     &                         AUX1 ,    NTX ,    AUX2 ,   KTRUNC     );

          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  drho_v(ih,il,it) = D(ihl) / (2.D0*pi);
              enddo
          enddo

      enddo






c-----Calculation of drho_s (protons + neutrons)
      drho_s = COMPLEX( 0.D0 , 0.D0 );
      do it = 3 , 3

          Ar = 0.D0;
          Ai = 0.D0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , ib2
                  if( dh_nnz(ib1,ib2) .eqv. .false. ) CYCLE;

                  j0 = ia_spx(ib2);
                  j1 = ia_spx(ib2)+id_spx(ib2)-1;
                  do j = j0 , j1
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);

                      i0 = ia_spx(ib1);
                      i1 = ia_spx(ib1)+id_spx(ib1)-1;
                      if( ib1 .eq. ib2 ) i1 = j;
                      do i = i0 , i1
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1.ne.fg2 ) CYCLE;
                          if( abs(ml1-ml2).ne.K_multipole ) CYCLE;

                          z = + drho_1(i,j,1) + drho_2(i,j,1)
     &                        + drho_1(i,j,2) + drho_2(i,j,2)
     &                        + drho_1(j,i,1) + drho_2(j,i,1)
     &                        + drho_1(j,i,2) + drho_2(j,i,2);
                          if( i .eq. j ) z = 0.5D0 * z;
                          if( fg1.eq.'g' .and. fg2.eq.'g' ) z = -z;

                          Ar(i,j) = DREAL(z);
                          Ai(i,j) = DIMAG(z);
                      enddo
                  enddo

              enddo
          enddo

          call diagUVtAUV(  N_total , NCOORD ,   k_PHI , N_blocks ,
     &                           Ar ,    NTX ,      Ai ,      NTX ,
     &                       dh_nnz ,    NBX ,
     &                       ia_spx , id_spx ,
     &                        PHI_U ,    NTX , PHI_SVt ,   KTRUNC , D ,
     &                         AUX1 ,    NTX ,    AUX2 ,   KTRUNC     );

          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  drho_s(ih,il) = D(ihl) / (2.D0*pi);
              enddo
          enddo

      enddo






c-----Calculation of dj_1
      dj_1 = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2

          Ar = 0.D0;
          Ai = 0.D0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , ib2
                  if( dh_nnz(ib1,ib2) .eqv. .false. ) CYCLE;

                  j0 = ia_spx(ib2);
                  j1 = ia_spx(ib2)+id_spx(ib2)-1;
                  do j = j0 , j1
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);

                      i0 = ia_spx(ib1);
                      i1 = ia_spx(ib1)+id_spx(ib1)-1;
                      if( ib1 .eq. ib2 ) i1 = j;
                      do i = i0 , i1
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1.eq.fg2 ) CYCLE;
                          if( ml1+ml2+1 .ne. +K_multipole ) CYCLE;

                          z = + drho_1(i,j,it) + drho_2(i,j,it)
     &                        - drho_1(j,i,it) - drho_2(j,i,it);
                          if( fg1.eq.'g' .and. fg2.eq.'f' ) z = -z;

                          Ar(i,j) = DREAL(z);
                          Ai(i,j) = DIMAG(z);
                      enddo
                  enddo

              enddo
          enddo

          call diagUVtAUV( N_total , NCOORD ,    k_PHI , N_blocks ,
     &                          Ar ,    NTX ,       Ai ,      NTX ,
     &                      dh_nnz ,    NBX ,
     &                      ia_spx , id_spx ,
     &                       PHI_U ,    NTX ,  PHI_SVt ,   KTRUNC , D ,
     &                        AUX1 ,    NTX ,     AUX2 ,   KTRUNC     );

          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  dj_1(ih,il,it) = D(ihl) / COMPLEX( 0.D0 , 2.D0*pi );

                  if( K_multipole.eq.0 ) then
                      dj_1(ih,il,it) = dj_1(ih,il,it) / DSQRT(2.D0);
                  else
                      dj_1(ih,il,it) = dj_1(ih,il,it) * DSQRT(2.D0);
                  endif
              enddo
          enddo

      enddo






c-----Calculation of dj_2
      dj_2 = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2

          Ar = 0.D0;
          Ai = 0.D0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , ib2
                  if( dh_nnz(ib1,ib2) .eqv. .false. ) CYCLE;

                  j0 = ia_spx(ib2);
                  j1 = ia_spx(ib2)+id_spx(ib2)-1;
                  do j = j0 , j1
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);

                      i0 = ia_spx(ib1);
                      i1 = ia_spx(ib1)+id_spx(ib1)-1;
                      if( ib1 .eq. ib2 ) i1 = j;
                      do i = i0 , i1
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1.eq.fg2 ) CYCLE;
                          if( ml1+ml2+1 .ne. -K_multipole ) CYCLE;

                          z = + drho_1(i,j,it) + drho_2(i,j,it)
     &                        - drho_1(j,i,it) - drho_2(j,i,it);
                          if( fg1.eq.'g' .and. fg2.eq.'f' ) z = -z;

                          Ar(i,j) = DREAL(z);
                          Ai(i,j) = DIMAG(z);
                      enddo
                  enddo

              enddo
          enddo

          call diagUVtAUV( N_total , NCOORD ,    k_PHI , N_blocks ,
     &                          Ar ,    NTX ,       Ai ,      NTX ,
     &                      dh_nnz ,    NBX ,
     &                      ia_spx , id_spx ,
     &                       PHI_U ,    NTX ,  PHI_SVt ,   KTRUNC , D ,
     &                        AUX1 ,    NTX ,     AUX2 ,   KTRUNC     );

          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  dj_2(ih,il,it) = D(ihl) / COMPLEX( 0.D0 , 2.D0*pi );

                  if( K_multipole.eq.0 ) then
                      dj_2(ih,il,it) = dj_2(ih,il,it) / DSQRT(2.D0);
                  else
                      dj_2(ih,il,it) = dj_2(ih,il,it) * DSQRT(2.D0);
                  endif
              enddo
          enddo

      enddo






c-----Calculation of dj_3
      dj_3 = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2

          Ar = 0.D0;
          Ai = 0.D0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , ib2
                  if( dh_nnz(ib1,ib2) .eqv. .false. ) CYCLE;

                  j0 = ia_spx(ib2);
                  j1 = ia_spx(ib2)+id_spx(ib2)-1;
                  do j = j0 , j1
                      fg2 = fg_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);

                      i0 = ia_spx(ib1);
                      i1 = ia_spx(ib1)+id_spx(ib1)-1;
                      if( ib1 .eq. ib2 ) i1 = j;
                      do i = i0 , i1
                          fg1 = fg_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);

                          if( fg1.eq.fg2 ) CYCLE;
                          if( abs(ml1-ml2).ne.K_multipole ) CYCLE;

                          z = + drho_1(i,j,it) - drho_2(i,j,it)
     &                        + drho_1(j,i,it) - drho_2(j,i,it);

                          Ar(i,j) = DREAL(z);
                          Ai(i,j) = DIMAG(z);
                      enddo
                  enddo

              enddo
          enddo

          call diagUVtAUV( N_total , NCOORD ,    k_PHI , N_blocks ,
     &                          Ar ,    NTX ,       Ai ,      NTX ,
     &                      dh_nnz ,    NBX ,
     &                      ia_spx , id_spx ,
     &                       PHI_U ,    NTX ,  PHI_SVt ,   KTRUNC , D ,
     &                        AUX1 ,    NTX ,     AUX2 ,   KTRUNC     );

          ihl = 0;
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;
                  ihl = ihl + 1;
                  dj_3(ih,il,it) = D(ihl) / (2.D0*pi);
              enddo
          enddo

      enddo






c-----Calculation of {dj_r,dj_p,dj_z}
      do it = 1 , 2
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  fac = 1.D0 / DSQRT(2.D0);

                  dj_r(ih,il,it) = fac*(+dj_1(ih,il,it)+dj_2(ih,il,it));

                  dj_p(ih,il,it) = fac*(-dj_1(ih,il,it)+dj_2(ih,il,it));

                  dj_z(ih,il,it) = dj_3(ih,il,it);

              enddo
          enddo
      enddo






c-----Calculation of the Laplacians for induced currents and densities
      call fam_laplacian( .false. );






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_ddensdcurr() ************************';
      write(6,*) '';
      endif

      return;
      end;






c======================================================================c

      subroutine fam_laplacian( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE quadrature;
      USE ddens;
      USE dcurr;
      USE dlaplace;
      USE famlaplacianmod;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;

      common /baspar/ hom, hb0, b0;
      common /defbas/ beta0, q, bp, bz;



      LOGICAL first_call /.true./;
      SAVE first_call;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_laplacian() ***********************';
      write(6,*) '';
      endif


      if( first_call .eqv. .true. ) then

          first_call = .false.;

          il = 0;
          do Nr = 0 , NSH/2
              do Nz = 0 , NSH - 2*Nr
                  il = il + 1;
                  NRR(il) = Nr;
                  NZZ(il) = Nz;
              enddo
          enddo
          call assert( il.eq.NSIZE , 'NSIZE too small' );

          do Nz = 0 , NSH
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;
                      z = DBLE(isign(1,ih)) * zb_fam(abs(ih));
                      b = b0*bz / DSQRT(2.D0);
                      phiz(ih,Nz) = phi_nz(Nz,b,z);
              enddo
          enddo

          do K = 0 , J_MAX+1
              do Nr = 0 , NSH/2
                  do il = 1 , NGL
                      r = rb_fam(il);
                      b = b0*bp / DSQRT(2.D0);
                      phirK(il,Nr,K) = phi_nr_ml(Nr,K,b,r);
                  enddo
              enddo
          enddo

      endif






c-----Calculation of Delta_{K} drho_vp(z,r)
      ldrho_vp = COMPLEX( 0.D0 , 0.D0 );
      do i = 1 , NSIZE
          Nz = NZZ(i);
          Nr = NRR(i);

          c(i) = COMPLEX( 0.D0 , 0.D0 );
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;
                  K = iabs( K_multipole );
                  c(i) = c(i) + wzwr(iabs(ih),il)
     &                        * phiz(ih,Nz)*phirK(il,Nr,K)
     &                        * drho_v(ih,il,2);
              enddo
          enddo

      enddo
      do i = 1 , NSIZE
          Nz = NZZ(i);
          Nr = NRR(i);

          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  z = DBLE(isign(1,ih)) * zb_fam(abs(ih));
                  r = rb_fam(il);
                  b = b0 / DSQRT(2.D0);
                  K = iabs( K_multipole );

                  fac = 0.D0;
                  fac = fac +            z**2.D0 / (b*bz)**4.D0;
                  fac = fac +            r**2.D0 / (b*bp)**4.D0;
                  fac = fac -       DBLE(2*Nz+1) / (b*bz)**2.D0;
                  fac = fac - DBLE(2*(2*Nr+K+1)) / (b*bp)**2.D0;
                  fac = fac * phiz(ih,Nz) * phirK(il,Nr,K);

                  ldrho_vp(ih,il) = ldrho_vp(ih,il) + c(i)*fac;
              enddo
          enddo
      enddo






c-----Calculation of Delta_{K} drho_s(z,r)
      ldrho_s = COMPLEX( 0.D0 , 0.D0 );
      do i = 1 , NSIZE
          Nz = NZZ(i);
          Nr = NRR(i);

          c(i) = COMPLEX( 0.D0 , 0.D0 );
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;
                  K = iabs( K_multipole );
                  c(i) = c(i) + wzwr(iabs(ih),il)
     &                        * phiz(ih,Nz)*phirK(il,Nr,K)
     &                        * drho_s(ih,il);
              enddo
          enddo

      enddo
      do i = 1 , NSIZE
          Nz = NZZ(i);
          Nr = NRR(i);

          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  z = DBLE(isign(1,ih)) * zb_fam(abs(ih));
                  r = rb_fam(il);
                  b = b0 / DSQRT(2.D0);
                  K = iabs( K_multipole );

                  fac = 0.D0;
                  fac = fac +            z**2.D0 / (b*bz)**4.D0;
                  fac = fac +            r**2.D0 / (b*bp)**4.D0;
                  fac = fac -       DBLE(2*Nz+1) / (b*bz)**2.D0;
                  fac = fac - DBLE(2*(2*Nr+K+1)) / (b*bp)**2.D0;
                  fac = fac * phiz(ih,Nz) * phirK(il,Nr,K);

                  ldrho_s(ih,il) = ldrho_s(ih,il) + c(i)*fac;
              enddo
          enddo
      enddo






c-----Calculation of Delta_{|K-1|} dj_1p(z,r)
      ldj_1p = COMPLEX( 0.D0 , 0.D0 );
      do i = 1 , NSIZE
          Nz = NZZ(i);
          Nr = NRR(i);

          c(i) = COMPLEX( 0.D0 , 0.D0 );
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;
                  K = iabs( K_multipole - 1 );
                  c(i) = c(i) + wzwr(iabs(ih),il)
     &                        * phiz(ih,Nz)*phirK(il,Nr,K)
     &                        * dj_1(ih,il,2);
              enddo
          enddo

      enddo
      do i = 1 , NSIZE
          Nz = NZZ(i);
          Nr = NRR(i);

          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  z = DBLE(isign(1,ih)) * zb_fam(abs(ih));
                  r = rb_fam(il);
                  b = b0 / DSQRT(2.D0);
                  K = iabs( K_multipole - 1 );

                  fac = 0.D0;
                  fac = fac +            z**2.D0 / (b*bz)**4.D0;
                  fac = fac +            r**2.D0 / (b*bp)**4.D0;
                  fac = fac -       DBLE(2*Nz+1) / (b*bz)**2.D0;
                  fac = fac - DBLE(2*(2*Nr+K+1)) / (b*bp)**2.D0;
                  fac = fac * phiz(ih,Nz) * phirK(il,Nr,K);

                  ldj_1p(ih,il) = ldj_1p(ih,il) + c(i)*fac;
              enddo
          enddo
      enddo






c-----Calculation of Delta_{|K+1|} dj_2p(z,r)
      ldj_2p = COMPLEX( 0.D0 , 0.D0 );
      do i = 1 , NSIZE
          Nz = NZZ(i);
          Nr = NRR(i);

          c(i) = COMPLEX( 0.D0 , 0.D0 );
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;
                  K = iabs( K_multipole + 1 );
                  c(i) = c(i) + wzwr(iabs(ih),il)
     &                        * phiz(ih,Nz)*phirK(il,Nr,K)
     &                        * dj_2(ih,il,2);
              enddo
          enddo

      enddo
      do i = 1 , NSIZE
          Nz = NZZ(i);
          Nr = NRR(i);

          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  z = DBLE(isign(1,ih)) * zb_fam(abs(ih));
                  r = rb_fam(il);
                  b = b0 / DSQRT(2.D0);
                  K = iabs( K_multipole + 1 );

                  fac = 0.D0;
                  fac = fac +            z**2.D0 / (b*bz)**4.D0;
                  fac = fac +            r**2.D0 / (b*bp)**4.D0;
                  fac = fac -       DBLE(2*Nz+1) / (b*bz)**2.D0;
                  fac = fac - DBLE(2*(2*Nr+K+1)) / (b*bp)**2.D0;
                  fac = fac * phiz(ih,Nz) * phirK(il,Nr,K);

                  ldj_2p(ih,il) = ldj_2p(ih,il) + c(i)*fac;
              enddo
          enddo
      enddo






c-----Calculation of Delta_{K} dj_3p(z,r)
      ldj_3p = COMPLEX( 0.D0 , 0.D0 );
      do i = 1 , NSIZE
          Nz = NZZ(i);
          Nr = NRR(i);

          c(i) = COMPLEX( 0.D0 , 0.D0 );
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;
                  K = iabs( K_multipole );
                  c(i) = c(i) + wzwr(iabs(ih),il)
     &                        * phiz(ih,Nz)*phirK(il,Nr,K)
     &                        * dj_3(ih,il,2);
              enddo
          enddo

      enddo
      do i = 1 , NSIZE
          Nz = NZZ(i);
          Nr = NRR(i);

          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih .eq. 0 ) CYCLE;

                  z = DBLE(isign(1,ih)) * zb_fam(abs(ih));
                  r = rb_fam(il);
                  b = b0 / DSQRT(2.D0);
                  K = iabs( K_multipole );

                  fac = 0.D0;
                  fac = fac +            z**2.D0 / (b*bz)**4.D0;
                  fac = fac +            r**2.D0 / (b*bp)**4.D0;
                  fac = fac -       DBLE(2*Nz+1) / (b*bz)**2.D0;
                  fac = fac - DBLE(2*(2*Nr+K+1)) / (b*bp)**2.D0;
                  fac = fac * phiz(ih,Nz) * phirK(il,Nr,K);

                  ldj_3p(ih,il) = ldj_3p(ih,il) + c(i)*fac;
              enddo
          enddo
      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_laplacian() *************************';
      write(6,*) '';
      endif

      return;
      end;






c======================================================================c

      subroutine diagUVtAUV(    N ,      M ,      K ,      NB ,
     &                         AR ,   LDAR ,
     &                         AI ,   LDAI ,
     &                       Annz , LDAnnz ,
     &                     ia_blk , id_blk ,
     &                          U ,    LDU ,
     &                          V ,    LDV ,      D ,
     &                       AUX1 , LDAUX1 ,
     &                       AUX2 , LDAUX2                      )

c======================================================================c
      ! Calculates the diagonal D = diag( (U*V)^T * A * (U*V) )
      ! elements of a complex upper triangular matrix A whose block
      ! structure is given by Annz
      !
      !  [in]: A         is complex  NxN  matrix, A = AR + j*AI
      !        Annz      is logical NBxNB matrix
      !        U         is real     NxK  matrix
      !        V         is real     KxM  matrix
      !        AUX1      is real     NxK  matrix (auxiliary)
      !        AUX2      is real     KxK  matrix (auxiliary)
      !        ia_blk(i) is the first index of i-th block
      !        id_blk(i) is the size of i-th block
      !
      !
      ! [out]: D is complex vector od length M

      IMPLICIT NONE;
      EXTERNAL         ddot;
      DOUBLE PRECISION ddot;

      INTEGER N, M, K, NB;
      INTEGER LDAR, LDAI, LDAnnz, LDU, LDV, LDAUX1, LDAUX2;
      INTEGER ia_blk( NB );
      INTEGER id_blk( NB );
      DOUBLE PRECISION AR( LDAR , N );
      DOUBLE PRECISION AI( LDAI , N );
      DOUBLE PRECISION  U(  LDU , K );
      DOUBLE PRECISION  V(  LDV , M );
      DOUBLE PRECISION AUX1( LDAUX1 , K );
      DOUBLE PRECISION AUX2( LDAUX2 , K );
      LOGICAL Annz( LDAnnz , NB );
      DOUBLE COMPLEX D( M );



      DOUBLE PRECISION Y( LDV );
      INTEGER i, ib, jb, i0, j0;

      D = COMPLEX( 0.D0 , 0.D0 );



c-----Real part
      ! Calculation of AUX1 = AR*U
      AUX1 = 0.D0;
      do ib = 1 , NB
          do jb = ib , NB
              i0 = ia_blk(ib);
              j0 = ia_blk(jb);
              if( Annz(ib,jb) ) then
                  call dgemm(        'N' ,   'N' ,
     &                        id_blk(ib) ,     K , id_blk(jb) ,
     &                              1.D0 ,
     &                         AR(i0,j0) ,  LDAR ,
     &                           U(j0,1) ,   LDU ,
     &                              1.D0 ,
     &                        AUX1(i0,1) , LDAUX1  );

              endif
          enddo
      enddo

      ! Calculation of AUX2 = U^T * AUX1 = U^T * AR * U
      call dgemm(    'T' ,    'N' ,
     &                 K ,      K ,      N ,
     &              1.D0 ,
     &                 U ,    LDU ,
     &              AUX1 , LDAUX1 ,
     &              0.D0 ,
     &              AUX2 , LDAUX2   );

      ! Calculation of < V_.i | AUX2 | V_.i >
      do i = 1 , M
          ! Calculation of Y = AUX2 * V_.i
          call dgemv(  'N' ,
     &                   K ,       K ,
     &                1.D0 ,
     &                AUX2 ,  LDAUX2 ,
     &              V(1,i) ,       1 ,
     &                0.D0 ,
     &                   Y ,       1  );

          ! Calculation of < V_.i | Y >
          D(i) = D(i) + COMPLEX( ddot(K,V(1,i),1,Y,1) , 0.D0 );
      enddo






c-----Imaginary part
      ! Calculation of AUX1 = AI*U
      AUX1 = 0.D0;
      do ib = 1 , NB
          do jb = ib , NB
              i0 = ia_blk(ib);
              j0 = ia_blk(jb);
              if( Annz(ib,jb) ) then
                  call dgemm(        'N' ,   'N' ,
     &                        id_blk(ib) ,     K , id_blk(jb) ,
     &                              1.D0 ,
     &                         AI(i0,j0) ,  LDAI ,
     &                           U(j0,1) ,   LDU ,
     &                              1.D0 ,
     &                        AUX1(i0,1) , LDAUX1  );

              endif
          enddo
      enddo

      ! Calculation of AUX2 = U^T * AUX1 = U^T * AI * U
      call dgemm(    'T' ,    'N' ,
     &                 K ,      K ,      N ,
     &              1.D0 ,
     &                 U ,    LDU ,
     &              AUX1 , LDAUX1 ,
     &              0.D0 ,
     &              AUX2 , LDAUX2   );

      ! Calculation of < V_.i | AUX2 | V_.i >
      do i = 1 , M
          ! Calculation of Y = AUX2 * V_.i
          call dgemv(  'N' ,
     &                   K ,       K ,
     &                1.D0 ,
     &                AUX2 ,  LDAUX2 ,
     &              V(1,i) ,       1 ,
     &                0.D0 ,
     &                   Y ,       1  );

          ! Calculation of < V_.i | Y >
          D(i) = D(i) + COMPLEX( 0.D0 , ddot(K,V(1,i),1,Y,1) );
      enddo






      return;
      end;
