!======================================================================!

      subroutine fam_ddelta()

!======================================================================!
      use fam_input;
      use simplex;
      use pairparams;
      use Wpairing;
      use dDelta;
      use dkappa;
      implicit none;

      integer                     :: it, ib1, ib2, i1, i2;
      character                   :: fg1, fg2;
      integer                     :: nz1, nz2;
      integer                     :: nr1, nr2;
      integer                     :: ml1, ml2;
      integer                     :: N_z, N_r;
      integer                     :: iW;
      double complex              :: dkappa_pl, dkappa_mi;
      double precision            :: normkappa;
      double precision , external :: frobeniousNormOfBlockMatrix;
      double complex              :: P_pl( 0:2*n0f , 0:2*n0f , 1:2 );
      double complex              :: P_mi( 0:2*n0f , 0:2*n0f , 1:2 );




      if( include_pairing == 0 ) then
          do it = 1 , 2
              call setToZeroBlockMatrix( dDelta1_pl(it) );
              call setToZeroBlockMatrix( dDelta1_mi(it) );
          end do
          return;
      end if




#ifdef DEBUG
      ! Selection rules test.
      ! (dkappa1_pl+dkappa1_pl^T)_{i1,i2} for fg1==fg2=='f' are selected by |ml1-ml2| = K.
      ! (dkappa1_mi+dkappa1_mi^T)_{i1,i2} for fg1==fg2=='f' are selected by |ml1-ml2| = K.
      do it = 1 , 2
          normkappa = max( frobeniousNormOfBlockMatrix(dkappa1_pl(it)) , &
                           frobeniousNormOfBlockMatrix(dkappa1_mi(it))   );

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  call assert( dkappa1_pl(it)%nnzblocks(ib1,ib2) .eqv. dkappa1_mi(it)%nnzblocks(ib1,ib2) , 'dkappa1_pl and dkappa1_mi nnzblocks error!' );

                  if( dkappa1_pl(it)%nnzblocks(ib1,ib2) .and. dkappa1_mi(it)%nnzblocks(ib1,ib2) ) then
                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , id_spx(ib1)

                              fg1 = fg_spx(ib1)%index(i1);
                              ml1 = ml_spx(ib1)%index(i1);

                              fg2 = fg_spx(ib2)%index(i2);
                              ml2 = ml_spx(ib2)%index(i2);

                              if( all([fg1,fg2]==['f','f']) .and. abs(ml1-ml2)/=K_multipole ) &
                                  call assert( abs( dkappa1_pl(it)%blocks(ib1,ib2)%mat(i1,i2) + dkappa1_pl(it)%blocks(ib2,ib1)%mat(i2,i1) ) / normkappa < 1.d-12 , 'dkappa1_pl+dkappa1_pl^T selection rules failed!' );

                              if( all([fg1,fg2]==['f','f']) .and. abs(ml1-ml2)/=K_multipole ) &
                                  call assert( abs( dkappa1_mi(it)%blocks(ib1,ib2)%mat(i1,i2) + dkappa1_mi(it)%blocks(ib2,ib1)%mat(i2,i1) ) / normkappa < 1.d-12 , 'dkappa1_mi+dkappa1_mi^T selection rules failed!' );

                          end do
                      end do
                  end if

              end do
          end do

      end do
#endif




      ! Calculation of P_pl(Nz,Nr) and P_mi(Nz,Nr), see Eq. (E.3).
      do it = 1 , 2

          P_pl(:,:,it) = cmplx( 0.d0 , 0.d0 , kind=8 );
          P_mi(:,:,it) = cmplx( 0.d0 , 0.d0 , kind=8 );

          iW = 0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , ib2
                  if( dDelta1_pl(it)%nnzblocks(ib1,ib2) .and. dDelta1_mi(it)%nnzblocks(ib1,ib2) ) then

                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , merge( i2 , id_spx(ib1) , ib1==ib2 ) ! Loop over the upper triangle.

                              fg1 = fg_spx(ib1)%index(i1);
                              nz1 = nz_spx(ib1)%index(i1);
                              nr1 = nr_spx(ib1)%index(i1);
                              ml1 = ml_spx(ib1)%index(i1);

                              fg2 = fg_spx(ib2)%index(i2);
                              nz2 = nz_spx(ib2)%index(i2);
                              nr2 = nr_spx(ib2)%index(i2);
                              ml2 = ml_spx(ib2)%index(i2);

                              if( all([fg1,fg2]==['f','f']) .and. abs(ml1-ml2)==K_multipole ) then

                                  dkappa_pl = merge( 0.5d0 , 1.d0 , all([ib1,i1]==[ib2,i2]) ) * ( dkappa1_pl(it)%blocks(ib1,ib2)%mat(i1,i2) + dkappa1_pl(it)%blocks(ib2,ib1)%mat(i2,i1) );
                                  dkappa_mi = merge( 0.5d0 , 1.d0 , all([ib1,i1]==[ib2,i2]) ) * ( dkappa1_mi(it)%blocks(ib1,ib2)%mat(i1,i2) + dkappa1_mi(it)%blocks(ib2,ib1)%mat(i2,i1) );

                                  do N_r = 0 , ( 2*nr1+abs(ml1) + 2*nr2+abs(ml2) - K_multipole )/2
                                      do N_z = 0 , nz1+nz2
                                          if( mod(N_z,2) == mod(nz1+nz2,2) ) then
                                              iW = iW + 1;

                                              P_pl(N_z,N_r,it) = P_pl(N_z,N_r,it) + W(iW) * dkappa_pl;
                                              P_mi(N_z,N_r,it) = P_mi(N_z,N_r,it) + W(iW) * dkappa_mi;

                                          end if
                                      end do
                                  end do

                              end if

                          end do
                      end do

                  end if
              end do
          end do
          call assert( iW == size(W) , 'W size error.' );

      end do




      ! Calculation of dDelta1_pl and dDelta1_mi, see Eq. (E.5).
      do it = 1 , 2

          iW = 0;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , ib2
                  if( dDelta1_pl(it)%nnzblocks(ib1,ib2) .and. dDelta1_mi(it)%nnzblocks(ib1,ib2) ) then

                      dDelta1_pl(it)%blocks(ib1,ib2)%mat(:,:) = cmplx( 0.d0 , 0.d0 , kind=8 );
                      dDelta1_mi(it)%blocks(ib1,ib2)%mat(:,:) = cmplx( 0.d0 , 0.d0 , kind=8 );

                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , merge( i2 , id_spx(ib1) , ib1==ib2 ) ! Loop over the upper triangle.

                              fg1 = fg_spx(ib1)%index(i1);
                              nz1 = nz_spx(ib1)%index(i1);
                              nr1 = nr_spx(ib1)%index(i1);
                              ml1 = ml_spx(ib1)%index(i1);

                              fg2 = fg_spx(ib2)%index(i2);
                              nz2 = nz_spx(ib2)%index(i2);
                              nr2 = nr_spx(ib2)%index(i2);
                              ml2 = ml_spx(ib2)%index(i2);

                              if( all([fg1,fg2]==['f','f']) .and. abs(ml1-ml2)==K_multipole ) then

                                  do N_r = 0 , ( 2*nr1+abs(ml1) + 2*nr2+abs(ml2) - K_multipole )/2
                                      do N_z = 0 , nz1+nz2
                                          if( mod(N_z,2) == mod(nz1+nz2,2) ) then
                                              iW = iW + 1;

                                              dDelta1_pl(it)%blocks(ib1,ib2)%mat(i1,i2) = dDelta1_pl(it)%blocks(ib1,ib2)%mat(i1,i2) + W(iW) * P_pl(N_z,N_r,it);
                                              dDelta1_mi(it)%blocks(ib1,ib2)%mat(i1,i2) = dDelta1_mi(it)%blocks(ib1,ib2)%mat(i1,i2) + W(iW) * P_mi(N_z,N_r,it);

                                          end if
                                      end do
                                  end do

                              end if

                          end do
                      end do

                      dDelta1_pl(it)%blocks(ib1,ib2)%mat(:,:) = - G_pairing * merge( 1.d0 , 0.5d0 , K_multipole==0 ) * dDelta1_pl(it)%blocks(ib1,ib2)%mat(:,:);
                      dDelta1_mi(it)%blocks(ib1,ib2)%mat(:,:) = - G_pairing * merge( 1.d0 , 0.5d0 , K_multipole==0 ) * dDelta1_mi(it)%blocks(ib1,ib2)%mat(:,:);

                  end if
              end do
          end do
          call assert( iW == size(W) , 'W size error.' );

          call constructFulldDelta1FromItsUpperTriangle( dDelta1_pl(it) );
          call constructFulldDelta1FromItsUpperTriangle( dDelta1_mi(it) );

      end do




      return;
      end subroutine fam_ddelta






!======================================================================!

      subroutine constructFulldDelta1FromItsUpperTriangle( dDelta1 )

!======================================================================!
      use dataTypes;
      implicit none;
      type(complexBlockMatrix) , intent(inout) :: dDelta1;
      integer                                  :: ib1, ib2, i1, i2;

      do ib2 = 1 , size(dDelta1%nnzblocks,2)
          do ib1 = ib2 , size(dDelta1%nnzblocks,1)
              if( dDelta1%nnzblocks(ib1,ib2) ) then
                  do i2 = 1 , size(dDelta1%blocks(ib1,ib2)%mat,2)
                      do i1 = merge( i2 , 1 , ib1==ib2 ) , size(dDelta1%blocks(ib1,ib2)%mat,1) ! Loop over the lower triangle.

                          dDelta1%blocks(ib1,ib2)%mat(i1,i2) = dDelta1%blocks(ib2,ib1)%mat(i2,i1);

                      end do
                  end do
              end if
          end do
      end do

      return;
      end subroutine constructFulldDelta1FromItsUpperTriangle
