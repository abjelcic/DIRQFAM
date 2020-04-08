c======================================================================c

      subroutine init_spurious( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE simplex;
      USE v_matrix;
      USE u_matrix;
      USE quadrature;
      USE basis;
      USE spurious;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;

      CHARACTER*2 nucnam;
      common /nucnuc/ amas, nneu, npro, nmas, nucnam;
      common /basnnn/ n0f, n0b;



      CHARACTER fg1, fg2;
      DOUBLE COMPLEX rcm( NTX , NTX ),
     &               pcm( NTX , NTX ),
     &               Tmp( NTX , NTX );
      LOGICAL     cm_nnz( NTX , NTX );
      DOUBLE COMPLEX trace;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN init_spurious() ***********************';
      write(6,*) '';
      endif


c--------------------------------------------------------------c
c     Rcm operator has the following form:                     c
c     Rcm = z/A for K=0, and Rcm = x/A for K=1.                c
c                                                              c
c     Pcm operator has the following form:                     c
c     Pcm = 1/i * d/dz for K=0, and Pcm = 1/i * d/dx for K=1.  c
c--------------------------------------------------------------c


      if( J_multipole.ne.1 .and. J_multipole.ne.3 ) then
          return;
      endif
      if( K_multipole.ne.0 .and. K_multipole.ne.1 ) then
          return;
      endif


      rcm    = COMPLEX( 0.D0 , 0.D0 );
      pcm    = COMPLEX( 0.D0 , 0.D0 );
      Tmp    = COMPLEX( 0.D0 , 0.D0 );
      cm_nnz = .false.;
      r20    = COMPLEX( 0.D0 , 0.D0 );
      p20    = COMPLEX( 0.D0 , 0.D0 );






c-----Calculation of cm_nnz
      do ib2 = 1 , N_blocks
          do ib1 = 1 , N_blocks

              j0 = ia_spx(ib2);
              j1 = ia_spx(ib2)+id_spx(ib2)-1;
              do j = j0 , j1
                  fg2 = fg_spx(j-j0+1,ib2);
                  nz2 = nz_spx(j-j0+1,ib2);
                  nr2 = nr_spx(j-j0+1,ib2);
                  ml2 = ml_spx(j-j0+1,ib2);

                  i0 = ia_spx(ib1);
                  i1 = ia_spx(ib1)+id_spx(ib1)-1;
                  do i = i0 , i1
                      fg1 = fg_spx(i-i0+1,ib1);
                      nz1 = nz_spx(i-i0+1,ib1);
                      nr1 = nr_spx(i-i0+1,ib1);
                      ml1 = ml_spx(i-i0+1,ib1);

                      if( K_multipole .eq. 0 ) then
                          if( fg1 .ne. fg2 ) CYCLE;
                          if( nr1 .ne. nr2 ) CYCLE;
                          if( abs(ml1-ml2) .ne. 0 ) CYCLE;
                          cm_nnz( ib1 , ib2 ) = .true.;
                      endif

                      if( K_multipole .eq. 1 ) then
                          if( fg1 .ne. fg2 ) CYCLE;
                          if( nz1 .ne. nz2 ) CYCLE;
                          if( abs(ml1-ml2) .ne. 1 ) CYCLE;
                          cm_nnz( ib1 , ib2 ) = .true.;
                      endif

                  enddo
               enddo

           enddo
       enddo






c-----Calculation of rcm
      if( K_multipole .eq. 0 ) then

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  j0 = ia_spx(ib2);
                  j1 = ia_spx(ib2)+id_spx(ib2)-1;
                  do j = j0 , j1
                      fg2 = fg_spx(j-j0+1,ib2);
                      nz2 = nz_spx(j-j0+1,ib2);
                      nr2 = nr_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);


                      i0 = ia_spx(ib1);
                      i1 = ia_spx(ib1)+id_spx(ib1)-1;
                      do i = i0 , i1
                          fg1 = fg_spx(i-i0+1,ib1);
                          nz1 = nz_spx(i-i0+1,ib1);
                          nr1 = nr_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);


                          if( fg1 .ne. fg2 ) CYCLE;
                          if( nr1 .ne. nr2 ) CYCLE;
                          if( ml1 .ne. ml2 ) CYCLE;


                          acc = 0.D0;
                          do ih = -NGH , +NGH
                              if( ih .eq. 0 ) CYCLE;
                              z = DBLE(isign(1,ih)) * zb_fam(abs(ih));

                              acc = acc + z
     &                                  * wz(abs(ih))
     &                                  * phi_z(ih,i)*phi_z(ih,j);
                          enddo
                          acc = acc / DBLE(nneu+npro);



                          rcm(i,j) = COMPLEX( acc , 0.D0 );

                      enddo
                   enddo

               enddo
           enddo

      else if( K_multipole .eq. 1 ) then

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  j0 = ia_spx(ib2);
                  j1 = ia_spx(ib2)+id_spx(ib2)-1;
                  do j = j0 , j1
                      fg2 = fg_spx(j-j0+1,ib2);
                      nz2 = nz_spx(j-j0+1,ib2);
                      nr2 = nr_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);


                      i0 = ia_spx(ib1);
                      i1 = ia_spx(ib1)+id_spx(ib1)-1;
                      do i = i0 , i1
                          fg1 = fg_spx(i-i0+1,ib1);
                          nz1 = nz_spx(i-i0+1,ib1);
                          nr1 = nr_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);


                          if( fg1 .ne. fg2 ) CYCLE;
                          if( nz1 .ne. nz2 ) CYCLE;
                          if( abs(ml1-ml2) .ne. 1 ) CYCLE;


                          acc = 0.D0;
                          do il = 1 , NGL
                              r = rb_fam(il);

                              acc = acc + r
     &                                  * wr(il)
     &                                  * phi_r(il,i)*phi_r(il,j);
                          enddo
                          acc =  0.5D0 * acc / DBLE(nneu+npro);



                          rcm(i,j) = COMPLEX( acc , 0.D0 );

                      enddo
                   enddo

               enddo
           enddo

      else
          stop 'Error: K_multipole wrong in init_spurious()!';
      endif






c-----Calculation of pcm
      if( K_multipole .eq. 0 ) then

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  j0 = ia_spx(ib2);
                  j1 = ia_spx(ib2)+id_spx(ib2)-1;
                  do j = j0 , j1
                      fg2 = fg_spx(j-j0+1,ib2);
                      nz2 = nz_spx(j-j0+1,ib2);
                      nr2 = nr_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);


                      i0 = ia_spx(ib1);
                      i1 = ia_spx(ib1)+id_spx(ib1)-1;
                      do i = i0 , i1
                          fg1 = fg_spx(i-i0+1,ib1);
                          nz1 = nz_spx(i-i0+1,ib1);
                          nr1 = nr_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);


                          if( fg1 .ne. fg2 ) CYCLE;
                          if( nr1 .ne. nr2 ) CYCLE;
                          if( ml1 .ne. ml2 ) CYCLE;


                          acc = 0.D0;
                          do ih = -NGH , +NGH
                              if( ih .eq. 0 ) CYCLE;

                              acc = acc + ( + phi_z(ih,i)*dphi_z(ih,j)
     &                                      - phi_z(ih,j)*dphi_z(ih,i) )
     &                                    * wz(abs(ih));
                          enddo
                          acc = -0.5D0 * acc;



                          pcm(i,j) = COMPLEX( 0.D0 , acc );

                      enddo
                   enddo

               enddo
           enddo

      else if( K_multipole .eq. 1 ) then

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  j0 = ia_spx(ib2);
                  j1 = ia_spx(ib2)+id_spx(ib2)-1;
                  do j = j0 , j1
                      fg2 = fg_spx(j-j0+1,ib2);
                      nz2 = nz_spx(j-j0+1,ib2);
                      nr2 = nr_spx(j-j0+1,ib2);
                      ml2 = ml_spx(j-j0+1,ib2);


                      i0 = ia_spx(ib1);
                      i1 = ia_spx(ib1)+id_spx(ib1)-1;
                      do i = i0 , i1
                          fg1 = fg_spx(i-i0+1,ib1);
                          nz1 = nz_spx(i-i0+1,ib1);
                          nr1 = nr_spx(i-i0+1,ib1);
                          ml1 = ml_spx(i-i0+1,ib1);


                          if( fg1 .ne. fg2 ) CYCLE;
                          if( nz1 .ne. nz2 ) CYCLE;
                          if( abs(ml1-ml2) .ne. 1 ) CYCLE;


                          acc = 0.D0;
                          do il = 1 , NGL
                              r = rb_fam(il);

                              acc = acc + ( + phi_r(il,i)*dphi_r(il,j)
     &                                      - phi_r(il,j)*dphi_r(il,i) )
     &                                    * wr(il);

                              acc = acc + wr(il)
     &                                  * DBLE(ml2*ml2-ml1*ml1)
     &                                  * phi_r(il,i)*phi_r(il,j) / r;
                          enddo
                          acc =  - 0.25D0 * acc;



                          pcm(i,j) = COMPLEX( 0.D0 , acc );

                      enddo
                   enddo

               enddo
           enddo

      else
          stop 'Error: K_multipole wrong in init_spurious()!';
      endif






#ifdef DEBUG
c-----Verify [rcm,pcm] = (i/A)*I
      if( .true. ) then

      ! Tmp = rcm*pcm
      do ib = 1 , N_blocks
          do jb = 1 , N_blocks

              do kb = 1 , N_blocks
                  if( cm_nnz(ib,kb) .and. cm_nnz(kb,jb) ) then

                  i0 = ia_spx(ib);
                  j0 = ia_spx(jb);
                  k0 = ia_spx(kb);

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib)  , id_spx(jb)  , id_spx(kb) ,
     &                        COMPLEX( +1.D0 , 0.D0 )   ,
     &                        rcm(i0,k0)                , NTX        ,
     &                        pcm(k0,j0)                , NTX        ,
     &                        COMPLEX( +1.D0 , 0.D0 )   ,
     &                        Tmp(i0,j0)                , NTX        );

                  endif
              enddo

          enddo
      enddo
      ! Tmp = rcm*pcm - pcm*rcm = [rcm,pcm]
      do ib = 1 , N_blocks
          do jb = 1 , N_blocks

              do kb = 1 , N_blocks
                  if( cm_nnz(ib,kb) .and. cm_nnz(kb,jb) ) then

                  i0 = ia_spx(ib);
                  j0 = ia_spx(jb);
                  k0 = ia_spx(kb);

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib)  , id_spx(jb)  , id_spx(kb) ,
     &                        COMPLEX( -1.D0 , 0.D0 )   ,
     &                        pcm(i0,k0)                , NTX        ,
     &                        rcm(k0,j0)                , NTX        ,
     &                        COMPLEX( +1.D0 , 0.D0 )   ,
     &                        Tmp(i0,j0)                , NTX        );

                  endif
              enddo

          enddo
      enddo
      ! Tmp = [rcm,pcm] - (i/A) * I
      do i = 1 , N_total
          Tmp(i,i) = Tmp(i,i) - COMPLEX( 0.D0 , 1.D0/DBLE(nneu+npro) );
      enddo

      do ib2 = 1 , N_blocks
          do ib1 = 1 , N_blocks

              j0 = ia_spx(ib2);
              j1 = ia_spx(ib2)+id_spx(ib2)-1;
              do j = j0 , j1
                  fg2  = fg_spx(j-j0+1,ib2);
                  nz2  = nz_spx(j-j0+1,ib2);
                  nr2  = nr_spx(j-j0+1,ib2);
                  ml2  = ml_spx(j-j0+1,ib2);
                  Nsh2 = 2*nr2 + abs(ml2) + nz2;

                  i0 = ia_spx(ib1);
                  i1 = ia_spx(ib1)+id_spx(ib1)-1;
                  do i = i0 , i1
                      fg1  = fg_spx(i-i0+1,ib1);
                      nz1  = nz_spx(i-i0+1,ib1);
                      nr1  = nr_spx(i-i0+1,ib1);
                      ml1  = ml_spx(i-i0+1,ib1);
                      Nsh1 = 2*nr1 + abs(ml1) + nz1;


                      if( fg1.eq.'f' .and. Nsh1.eq.n0f   ) CYCLE;
                      if( fg1.eq.'g' .and. Nsh1.eq.n0f+1 ) CYCLE;

                      if( fg2.eq.'f' .and. Nsh2.eq.n0f   ) CYCLE;
                      if( fg2.eq.'g' .and. Nsh2.eq.n0f+1 ) CYCLE;


                      call assert( ABS(Tmp(i,j)) .lt. 1.D-12 ,
     &                             '[rcm,pcm] =/= (i/A)*I'     );


                  enddo
               enddo

           enddo
       enddo

      endif
#endif






c-----Calculation of r20 = hermconjg(u)*rcm*v + h.c.
      do it = 1 , 2

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  if( cm_nnz(ib1,ib2) .eqv. .false. ) CYCLE;

                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);

                  call zgemm( 'c'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        u(1,1,ib1,it)             , NBSX        ,
     &                        rcm(i0,j0)                , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Tmp(i0,j0)                , NTX         );

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Tmp(i0,j0)                , NTX         ,
     &                        v(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        r20(i0,j0,it)             , NTX         );


              enddo
          enddo

          do i = 1 , N_total
              do j = i , N_total
                  r20(i,j,it) = r20(i,j,it) + DCONJG(r20(j,i,it));
              enddo
          enddo
          do i = 1 , N_total
              do j = 1 , i-1
                  r20(i,j,it) = DCONJG( r20(j,i,it) );
              enddo
          enddo

      enddo






c-----Calculation of p20 = hermconjg(u)*pcm*v - h.c.
      do it = 1 , 2

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  if( cm_nnz(ib1,ib2) .eqv. .false. ) CYCLE;

                  i0 = ia_spx(ib1);
                  j0 = ia_spx(ib2);

                  call zgemm( 'c'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib1) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        u(1,1,ib1,it)             , NBSX        ,
     &                        pcm(i0,j0)                , NTX         ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        Tmp(i0,j0)                , NTX         );

                  call zgemm( 'n'         , 'n'         ,
     &                        id_spx(ib1) , id_spx(ib2) , id_spx(ib2) ,
     &                        COMPLEX( 1.D0 , 0.D0 )    ,
     &                        Tmp(i0,j0)                , NTX         ,
     &                        v(1,1,ib2,it)             , NBSX        ,
     &                        COMPLEX( 0.D0 , 0.D0 )    ,
     &                        p20(i0,j0,it)             , NTX         );


              enddo
          enddo

          do i = 1 , N_total
              do j = i , N_total
                  p20(i,j,it) = p20(i,j,it) - DCONJG(p20(j,i,it));
              enddo
          enddo
          do i = 1 , N_total
              do j = 1 , i-1
                  p20(i,j,it) = - DCONJG( p20(j,i,it) );
              enddo
          enddo

      enddo






c-----Calculation of <Phi| [Rcm,Pcm] |Phi> = 2*Tr[ hermconj(r20)*p20 ]
      do it = 1 , 2

          trace = COMPLEX( 0.D0 , 0.D0 );
          do j = 1 , N_total
              do i = 1 , N_total
                  trace = trace + DCONJG(r20(i,j,it)) * p20(i,j,it);
              enddo
          enddo

          RcmPcm_commutator(it) = 2.D0 * trace;
      enddo
#ifdef DEBUG
      write(6,*) '[Rcm,Pcm] = ',
     &           RcmPcm_commutator(1) + RcmPcm_commutator(2);
#endif





      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END init_spurious() *************************';
      write(6,*) '';
      endif

      return;
      end;
