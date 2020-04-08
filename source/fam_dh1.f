c======================================================================c

      subroutine fam_dh1( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE simplex;
      USE wbasis;
      USE dpot;
      USE dh;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      CHARACTER fg1, fg2;
      DOUBLE COMPLEX z;
      DOUBLE COMPLEX dVpS_odd  ( 1:NGH , 1:NGL , 2 ),
     &               dVpS_evn  ( 1:NGH , 1:NGL , 2 ),
     &               dVmS_odd  ( 1:NGH , 1:NGL , 2 ),
     &               dVmS_evn  ( 1:NGH , 1:NGL , 2 ),
     &               dSig_z_odd( 1:NGH , 1:NGL , 2 ),
     &               dSig_z_evn( 1:NGH , 1:NGL , 2 ),
     &               dSig_r_odd( 1:NGH , 1:NGL , 2 ),
     &               dSig_r_evn( 1:NGH , 1:NGL , 2 ),
     &               dSig_p_odd( 1:NGH , 1:NGL , 2 ),
     &               dSig_p_evn( 1:NGH , 1:NGL , 2 );



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_dh1() *****************************';
      write(6,*) '';
      endif






c-----Initializing odd/even parts of the induced potentials
      do it = 1 , 2
        do il = 1 , NGL
          do ih = 1 , NGH

              dVpS_odd(ih,il,it) =   dVpS(ih,il,it) -   dVps(-ih,il,it);
              dVpS_evn(ih,il,it) =   dVpS(ih,il,it) +   dVps(-ih,il,it);

              dVmS_odd(ih,il,it) =   dVmS(ih,il,it) -   dVms(-ih,il,it);
              dVmS_evn(ih,il,it) =   dVmS(ih,il,it) +   dVms(-ih,il,it);

            dSig_z_odd(ih,il,it) = dSig_z(ih,il,it) - dSig_z(-ih,il,it);
            dSig_z_evn(ih,il,it) = dSig_z(ih,il,it) + dSig_z(-ih,il,it);

            dSig_r_odd(ih,il,it) = dSig_r(ih,il,it) - dSig_r(-ih,il,it);
            dSig_r_evn(ih,il,it) = dSig_r(ih,il,it) + dSig_r(-ih,il,it);

            dSig_p_odd(ih,il,it) = dSig_p(ih,il,it) - dSig_p(-ih,il,it);
            dSig_p_evn(ih,il,it) = dSig_p(ih,il,it) + dSig_p(-ih,il,it);

          enddo
        enddo
      enddo






c-----Calculation of the upper triangle of dh_1 matrix
      dh_1 = COMPLEX( 0.D0 , 0.D0 );
      K    = K_multipole;
      do it = 1 , 2
        do ib2 = 1 , N_blocks
          do ib1 = 1 , ib2

            j0 = ia_spx(ib2);
            j1 = ia_spx(ib2)+id_spx(ib2)-1;
            do j = j0 , j1
              fg2 = fg_spx(j-j0+1,ib2);
              nz2 = nz_spx(j-j0+1,ib2);
              ml2 = ml_spx(j-j0+1,ib2);


              i0 = ia_spx(ib1);
              i1 = ia_spx(ib1)+id_spx(ib1)-1;
              if( ib1 .eq. ib2 ) i1 = j;
              do i = i0 , i1
                fg1 = fg_spx(i-i0+1,ib1);
                nz1 = nz_spx(i-i0+1,ib1);
                ml1 = ml_spx(i-i0+1,ib1);




                if( fg1.eq.'f' .and. fg2.eq.'f' ) then

                    if( abs(ml1-ml2) .ne. K ) CYCLE;

                    z = COMPLEX( 0.D0 , 0.D0 );
                    if( MOD(nz1+nz2,2) .eq. 0 ) then
                        do il = 1 , NGL
                            do ih = 1 , NGH
                                z = z + dVpS_evn(ih,il,it)
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                            enddo
                        enddo
                    else
                        do il = 1 , NGL
                            do ih = 1 , NGH
                                z = z + dVpS_odd(ih,il,it)
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                            enddo
                        enddo
                    endif

                    if( K .ne. 0 ) z = 0.5D0 * z;

                    dh_1(i,j,it) = z;

                endif


                if( fg1.eq.'g' .and. fg2.eq.'g' ) then

                    if( abs(ml1-ml2) .ne. K ) CYCLE;

                    z = COMPLEX( 0.D0 , 0.D0 );
                    if( MOD(nz1+nz2,2) .eq. 0 ) then
                        do il = 1 , NGL
                            do ih = 1 , NGH
                                z = z + dVmS_evn(ih,il,it)
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                            enddo
                        enddo
                    else
                        do il = 1 , NGL
                            do ih = 1 , NGH
                                z = z + dVmS_odd(ih,il,it)
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                            enddo
                        enddo
                    endif

                    if( K .ne. 0 ) z = 0.5D0 * z;

                    dh_1(i,j,it) = z;

                endif


                if( fg1.eq.'g' .and. fg2.eq.'f' ) then

                    if( abs(ml1+ml2+1) .eq. K ) then

                        z = COMPLEX( 0.D0 , 0.D0 );

                        if( K .eq. 0 ) then
                          if( MOD(nz1+nz2,2) .eq. 0 ) then
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + dSig_r_evn(ih,il,it)
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                              enddo
                            enddo
                          else
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + dSig_r_odd(ih,il,it)
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                              enddo
                            enddo
                          endif
                        endif

                        if( K.ne.0 .and. (ml1+ml2+1).eq.(+K) ) then
                          if( MOD(nz1+nz2,2) .eq. 0 ) then
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + ( + dSig_r_evn(ih,il,it)
     &                                    - dSig_p_evn(ih,il,it) )
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                              enddo
                            enddo
                          else
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + ( + dSig_r_odd(ih,il,it)
     &                                    - dSig_p_odd(ih,il,it) )
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                              enddo
                            enddo
                          endif
                        endif

                        if( K.ne.0 .and. (ml1+ml2+1).eq.(-K) ) then
                          if( MOD(nz1+nz2,2) .eq. 0 ) then
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + ( + dSig_r_evn(ih,il,it)
     &                                    + dSig_p_evn(ih,il,it) )
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                              enddo
                            enddo
                          else
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + ( + dSig_r_odd(ih,il,it)
     &                                    + dSig_p_odd(ih,il,it) )
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                              enddo
                            enddo
                          endif
                        endif


                        z = z * COMPLEX( 0.D0 , 1.D0 );
                        if( K .ne. 0 ) z = 0.5D0 * z;

                        dh_1(i,j,it) = + z;

                    endif

                    if( abs(ml1-ml2) .eq. K ) then

                        z = COMPLEX( 0.D0 , 0.D0 );
                        if( MOD(nz1+nz2,2) .eq. 0 ) then
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + dSig_z_evn(ih,il,it)
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                              enddo
                            enddo
                        else
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + dSig_z_odd(ih,il,it)
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                              enddo
                            enddo
                        endif

                        if( K .ne. 0 ) z = 0.5D0 * z;

                        dh_1(i,j,it) = - z;

                    endif

                endif


                if( fg1.eq.'f' .and. fg2.eq.'g' ) then

                    if( abs(ml1+ml2+1) .eq. K ) then

                        z = COMPLEX( 0.D0 , 0.D0 );

                        if( K .eq. 0 ) then
                          if( MOD(nz1+nz2,2) .eq. 0 ) then
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + dSig_r_evn(ih,il,it)
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                              enddo
                            enddo
                          else
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + dSig_r_odd(ih,il,it)
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                              enddo
                            enddo
                          endif
                        endif

                        if( K.ne.0 .and. (ml1+ml2+1).eq.(+K) ) then
                          if( MOD(nz1+nz2,2) .eq. 0 ) then
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + ( + dSig_r_evn(ih,il,it)
     &                                    - dSig_p_evn(ih,il,it) )
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                              enddo
                            enddo
                          else
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + ( + dSig_r_odd(ih,il,it)
     &                                    - dSig_p_odd(ih,il,it) )
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                              enddo
                            enddo
                          endif
                        endif

                        if( K.ne.0 .and. (ml1+ml2+1).eq.(-K) ) then
                          if( MOD(nz1+nz2,2) .eq. 0 ) then
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + ( + dSig_r_evn(ih,il,it)
     &                                    + dSig_p_evn(ih,il,it) )
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                              enddo
                            enddo
                          else
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + ( + dSig_r_odd(ih,il,it)
     &                                    + dSig_p_odd(ih,il,it) )
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                              enddo
                            enddo
                          endif
                        endif


                        z = z * COMPLEX( 0.D0 , 1.D0 );
                        if( K .ne. 0 ) z = 0.5D0 * z;

                        dh_1(i,j,it) = - z;

                    endif

                    if( abs(ml1-ml2) .eq. K ) then

                        z = COMPLEX( 0.D0 , 0.D0 );
                        if( MOD(nz1+nz2,2) .eq. 0 ) then
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + dSig_z_evn(ih,il,it)
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                              enddo
                            enddo
                        else
                            do il = 1 , NGL
                              do ih = 1 , NGH
                                z = z + dSig_z_odd(ih,il,it)
     &                                * wPhi(ih,il,i)*wPhi(ih,il,j);
                              enddo
                            enddo
                        endif

                        if( K .ne. 0 ) z = 0.5D0 * z;

                        dh_1(i,j,it) = - z;

                    endif

                endif


              enddo
            enddo

          enddo
        enddo
      enddo






c-----Construction of full dh_1 matrix
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






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_dh1() *******************************';
      write(6,*) '';
      endif

      return;
      end;
