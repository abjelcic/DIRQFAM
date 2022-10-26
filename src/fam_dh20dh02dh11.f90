!======================================================================!

      subroutine fam_dh20dh02()

!======================================================================!
      use simplex;
      use u_matrix;
      use v_matrix;
      use dh02dh20matrix;
      use dh;
      use dDelta;
      use tempBlockMatrix;
      implicit none;

      double complex , parameter :: cone  = cmplx( 1.d0 , 0.d0 , kind=8 );
      double complex , parameter :: czero = cmplx( 0.d0 , 0.d0 , kind=8 );
      integer                    :: it, ib1, ib2;




      ! Calculation of dh20 = - u' * dh1         * v
      !                       - v' * dh2^T       * u
      !                       + u' * dDelta1_pl  * u
      !                       - v' * dDelta1_mi' * v , see Eq. (B.14).
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh20(it)%nnzblocks(ib1,ib2) ) then

                      dh20(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = czero;

                      if( dh1(it)%nnzblocks(ib1,ib2) ) then

                          call zgemm( 'c' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                      , &
                                      -cone                                                                    , &
                                      u(it)%blocks(ib1,ib1)%mat(1,1)    , size(u(it)%blocks(ib1,ib1)%mat,1)    , &
                                      dh1(it)%blocks(ib1,ib2)%mat(1,1)  , size(dh1(it)%blocks(ib1,ib2)%mat,1)  , &
                                      czero                                                                    , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)   , size(tmpMat%blocks(ib1,ib2)%mat,1)     );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                      , &
                                      +cone                                                                    , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)   , size(tmpMat%blocks(ib1,ib2)%mat,1)   , &
                                      v(it)%blocks(ib2,ib2)%mat(1,1)    , size(v(it)%blocks(ib2,ib2)%mat,1)    , &
                                      cone                                                                     , &
                                      dh20(it)%blocks(ib1,ib2)%mat(1,1) , size(dh20(it)%blocks(ib1,ib2)%mat,1)   );

                      end if

                      if( dh2(it)%nnzblocks(ib2,ib1) ) then ! Notice here (ib2,ib1).

                          call zgemm( 'c' , 't' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                       , &
                                      -cone                                                                     , &
                                      v(it)%blocks(ib1,ib1)%mat(1,1)    , size(v(it)%blocks(ib1,ib1)%mat,1)     , &
                                      dh2(it)%blocks(ib2,ib1)%mat(1,1)  , size(dh2(it)%blocks(ib2,ib1)%mat,1)   , & ! Notice here (ib2,ib1).
                                      czero                                                                     , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)   , size(tmpMat%blocks(ib1,ib2)%mat,1)      );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                       , &
                                      +cone                                                                     , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)   , size(tmpMat%blocks(ib1,ib2)%mat,1)    , &
                                      u(it)%blocks(ib2,ib2)%mat(1,1)    , size(u(it)%blocks(ib2,ib2)%mat,1)     , &
                                      cone                                                                      , &
                                      dh20(it)%blocks(ib1,ib2)%mat(1,1) , size(dh20(it)%blocks(ib1,ib2)%mat,1)    );

                      end if

                      if( dDelta1_pl(it)%nnzblocks(ib1,ib2) ) then

                          call zgemm( 'c' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                                  , &
                                      +cone                                                                                , &
                                      u(it)%blocks(ib1,ib1)%mat(1,1)          , size(u(it)%blocks(ib1,ib1)%mat,1)          , &
                                      dDelta1_pl(it)%blocks(ib1,ib2)%mat(1,1) , size(dDelta1_pl(it)%blocks(ib1,ib2)%mat,1) , &
                                      czero                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)           );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                                  , &
                                      +cone                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)         , &
                                      u(it)%blocks(ib2,ib2)%mat(1,1)          , size(u(it)%blocks(ib2,ib2)%mat,1)          , &
                                      cone                                                                                 , &
                                      dh20(it)%blocks(ib1,ib2)%mat(1,1)       , size(dh20(it)%blocks(ib1,ib2)%mat,1)         );

                      end if

                      if( dDelta1_mi(it)%nnzblocks(ib2,ib1) ) then ! Notice here (ib2,ib1).

                          call zgemm( 'c' , 'c' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                                  , &
                                      -cone                                                                                , &
                                      v(it)%blocks(ib1,ib1)%mat(1,1)          , size(v(it)%blocks(ib1,ib1)%mat,1)          , &
                                      dDelta1_mi(it)%blocks(ib2,ib1)%mat(1,1) , size(dDelta1_mi(it)%blocks(ib2,ib1)%mat,1) , & ! Notice here (ib2,ib1).
                                      czero                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)           );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                                  , &
                                      +cone                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)         , &
                                      v(it)%blocks(ib2,ib2)%mat(1,1)          , size(v(it)%blocks(ib2,ib2)%mat,1)          , &
                                      cone                                                                                 , &
                                      dh20(it)%blocks(ib1,ib2)%mat(1,1)       , size(dh20(it)%blocks(ib1,ib2)%mat,1)         );

                      end if

                  end if
              end do
          end do
      end do




      ! Calculation of dh02 = ( - v' * dh1         * u
      !                         - u' * dh2^T       * v
      !                         - v' * dDelta1_pl  * v
      !                         + u' * dDelta1_mi' * u )^T , see Eq. (B.15).
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh02(it)%nnzblocks(ib1,ib2) ) then

                      dh02(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = czero;

                      if( dh1(it)%nnzblocks(ib1,ib2) ) then

                          call zgemm( 'c' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                      , &
                                      -cone                                                                    , &
                                      v(it)%blocks(ib1,ib1)%mat(1,1)    , size(v(it)%blocks(ib1,ib1)%mat,1)    , &
                                      dh1(it)%blocks(ib1,ib2)%mat(1,1)  , size(dh1(it)%blocks(ib1,ib2)%mat,1)  , &
                                      czero                                                                    , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)   , size(tmpMat%blocks(ib1,ib2)%mat,1)     );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                      , &
                                      +cone                                                                    , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)   , size(tmpMat%blocks(ib1,ib2)%mat,1)   , &
                                      u(it)%blocks(ib2,ib2)%mat(1,1)    , size(u(it)%blocks(ib2,ib2)%mat,1)    , &
                                      cone                                                                     , &
                                      dh02(it)%blocks(ib1,ib2)%mat(1,1) , size(dh02(it)%blocks(ib1,ib2)%mat,1)   );

                      end if

                      if( dh2(it)%nnzblocks(ib2,ib1) ) then ! Notice here (ib2,ib1).

                          call zgemm( 'c' , 't' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                       , &
                                      -cone                                                                     , &
                                      u(it)%blocks(ib1,ib1)%mat(1,1)    , size(u(it)%blocks(ib1,ib1)%mat,1)     , &
                                      dh2(it)%blocks(ib2,ib1)%mat(1,1)  , size(dh2(it)%blocks(ib2,ib1)%mat,1)   , & ! Notice here (ib2,ib1).
                                      czero                                                                     , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)   , size(tmpMat%blocks(ib1,ib2)%mat,1)      );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                       , &
                                      +cone                                                                     , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)   , size(tmpMat%blocks(ib1,ib2)%mat,1)    , &
                                      v(it)%blocks(ib2,ib2)%mat(1,1)    , size(v(it)%blocks(ib2,ib2)%mat,1)     , &
                                      cone                                                                      , &
                                      dh02(it)%blocks(ib1,ib2)%mat(1,1) , size(dh02(it)%blocks(ib1,ib2)%mat,1)    );

                      end if

                      if( dDelta1_pl(it)%nnzblocks(ib1,ib2) ) then

                          call zgemm( 'c' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                                  , &
                                      -cone                                                                                , &
                                      v(it)%blocks(ib1,ib1)%mat(1,1)          , size(v(it)%blocks(ib1,ib1)%mat,1)          , &
                                      dDelta1_pl(it)%blocks(ib1,ib2)%mat(1,1) , size(dDelta1_pl(it)%blocks(ib1,ib2)%mat,1) , &
                                      czero                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)           );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                                  , &
                                      +cone                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)         , &
                                      v(it)%blocks(ib2,ib2)%mat(1,1)          , size(v(it)%blocks(ib2,ib2)%mat,1)          , &
                                      cone                                                                                 , &
                                      dh02(it)%blocks(ib1,ib2)%mat(1,1)       , size(dh02(it)%blocks(ib1,ib2)%mat,1)         );

                      end if

                      if( dDelta1_mi(it)%nnzblocks(ib2,ib1) ) then ! Notice here (ib2,ib1).

                          call zgemm( 'c' , 'c' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                                  , &
                                      +cone                                                                                , &
                                      u(it)%blocks(ib1,ib1)%mat(1,1)          , size(u(it)%blocks(ib1,ib1)%mat,1)          , &
                                      dDelta1_mi(it)%blocks(ib2,ib1)%mat(1,1) , size(dDelta1_mi(it)%blocks(ib2,ib1)%mat,1) , & ! Notice here (ib2,ib1).
                                      czero                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)           );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                                  , &
                                      +cone                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)         , &
                                      u(it)%blocks(ib2,ib2)%mat(1,1)          , size(u(it)%blocks(ib2,ib2)%mat,1)          , &
                                      cone                                                                                 , &
                                      dh02(it)%blocks(ib1,ib2)%mat(1,1)       , size(dh02(it)%blocks(ib1,ib2)%mat,1)         );

                      end if

                  end if
              end do
          end do

          call transposeBlockMatrix( dh02(it) );

      end do




      return;
      end subroutine fam_dh20dh02






!======================================================================!

      subroutine calculate_dh11()

!======================================================================!
      use simplex;
      use u_matrix;
      use v_matrix;
      use dh11matrix;
      use dh;
      use dDelta;
      use tempBlockMatrix;
      implicit none;

      double complex , parameter :: cone  = cmplx( 1.d0 , 0.d0 , kind=8 );
      double complex , parameter :: czero = cmplx( 0.d0 , 0.d0 , kind=8 );
      integer                    :: it, ib1, ib2;




      ! Calculation of dh11_1 = + u' * dh1         * u
      !                         - v' * dh2^T       * v
      !                         + u' * dDelta1_pl  * v
      !                         + v' * dDelta1_mi' * u .
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh11_1(it)%nnzblocks(ib1,ib2) ) then

                      dh11_1(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = czero;

                      if( dh1(it)%nnzblocks(ib1,ib2) ) then

                          call zgemm( 'c' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                          , &
                                      +cone                                                                        , &
                                      u(it)%blocks(ib1,ib1)%mat(1,1)      , size(u(it)%blocks(ib1,ib1)%mat,1)      , &
                                      dh1(it)%blocks(ib1,ib2)%mat(1,1)    , size(dh1(it)%blocks(ib1,ib2)%mat,1)    , &
                                      czero                                                                        , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)     , size(tmpMat%blocks(ib1,ib2)%mat,1)       );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                          , &
                                      +cone                                                                        , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)     , size(tmpMat%blocks(ib1,ib2)%mat,1)     , &
                                      u(it)%blocks(ib2,ib2)%mat(1,1)      , size(u(it)%blocks(ib2,ib2)%mat,1)      , &
                                      cone                                                                         , &
                                      dh11_1(it)%blocks(ib1,ib2)%mat(1,1) , size(dh11_1(it)%blocks(ib1,ib2)%mat,1)   );

                      end if

                      if( dh2(it)%nnzblocks(ib2,ib1) ) then ! Notice here (ib2,ib1).

                          call zgemm( 'c' , 't' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                           , &
                                      -cone                                                                         , &
                                      v(it)%blocks(ib1,ib1)%mat(1,1)      , size(v(it)%blocks(ib1,ib1)%mat,1)       , &
                                      dh2(it)%blocks(ib2,ib1)%mat(1,1)    , size(dh2(it)%blocks(ib2,ib1)%mat,1)     , & ! Notice here (ib2,ib1).
                                      czero                                                                         , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)     , size(tmpMat%blocks(ib1,ib2)%mat,1)        );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                           , &
                                      +cone                                                                         , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)     , size(tmpMat%blocks(ib1,ib2)%mat,1)      , &
                                      v(it)%blocks(ib2,ib2)%mat(1,1)      , size(v(it)%blocks(ib2,ib2)%mat,1)       , &
                                      cone                                                                          , &
                                      dh11_1(it)%blocks(ib1,ib2)%mat(1,1) , size(dh11_1(it)%blocks(ib1,ib2)%mat,1)    );

                      end if

                      if( dDelta1_pl(it)%nnzblocks(ib1,ib2) ) then

                          call zgemm( 'c' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                                  , &
                                      +cone                                                                                , &
                                      u(it)%blocks(ib1,ib1)%mat(1,1)          , size(u(it)%blocks(ib1,ib1)%mat,1)          , &
                                      dDelta1_pl(it)%blocks(ib1,ib2)%mat(1,1) , size(dDelta1_pl(it)%blocks(ib1,ib2)%mat,1) , &
                                      czero                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)           );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                                  , &
                                      +cone                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)         , &
                                      v(it)%blocks(ib2,ib2)%mat(1,1)          , size(v(it)%blocks(ib2,ib2)%mat,1)          , &
                                      cone                                                                                 , &
                                      dh11_1(it)%blocks(ib1,ib2)%mat(1,1)     , size(dh11_1(it)%blocks(ib1,ib2)%mat,1)       );

                      end if

                      if( dDelta1_mi(it)%nnzblocks(ib2,ib1) ) then ! Notice here (ib2,ib1).

                          call zgemm( 'c' , 'c' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                                  , &
                                      +cone                                                                                , &
                                      v(it)%blocks(ib1,ib1)%mat(1,1)          , size(v(it)%blocks(ib1,ib1)%mat,1)          , &
                                      dDelta1_mi(it)%blocks(ib2,ib1)%mat(1,1) , size(dDelta1_mi(it)%blocks(ib2,ib1)%mat,1) , & ! Notice here (ib2,ib1).
                                      czero                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)           );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                                  , &
                                      +cone                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)         , &
                                      u(it)%blocks(ib2,ib2)%mat(1,1)          , size(u(it)%blocks(ib2,ib2)%mat,1)          , &
                                      cone                                                                                 , &
                                      dh11_1(it)%blocks(ib1,ib2)%mat(1,1)     , size(dh11_1(it)%blocks(ib1,ib2)%mat,1)       );

                      end if

                  end if
              end do
          end do
      end do




      ! Calculation of dh11_2 = ( - v' * dh1         * v
      !                           + u' * dh2^T       * u
      !                           + v' * dDelta1_pl  * u
      !                           + u' * dDelta1_mi' * v )^T .
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh11_2(it)%nnzblocks(ib1,ib2) ) then

                      dh11_2(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = czero;

                      if( dh1(it)%nnzblocks(ib1,ib2) ) then

                          call zgemm( 'c' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                          , &
                                      -cone                                                                        , &
                                      v(it)%blocks(ib1,ib1)%mat(1,1)      , size(v(it)%blocks(ib1,ib1)%mat,1)      , &
                                      dh1(it)%blocks(ib1,ib2)%mat(1,1)    , size(dh1(it)%blocks(ib1,ib2)%mat,1)    , &
                                      czero                                                                        , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)     , size(tmpMat%blocks(ib1,ib2)%mat,1)       );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                          , &
                                      +cone                                                                        , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)     , size(tmpMat%blocks(ib1,ib2)%mat,1)     , &
                                      v(it)%blocks(ib2,ib2)%mat(1,1)      , size(v(it)%blocks(ib2,ib2)%mat,1)      , &
                                      cone                                                                         , &
                                      dh11_2(it)%blocks(ib1,ib2)%mat(1,1) , size(dh11_2(it)%blocks(ib1,ib2)%mat,1)   );

                      end if

                      if( dh2(it)%nnzblocks(ib2,ib1) ) then ! Notice here (ib2,ib1).

                          call zgemm( 'c' , 't' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                           , &
                                      +cone                                                                         , &
                                      u(it)%blocks(ib1,ib1)%mat(1,1)      , size(u(it)%blocks(ib1,ib1)%mat,1)       , &
                                      dh2(it)%blocks(ib2,ib1)%mat(1,1)    , size(dh2(it)%blocks(ib2,ib1)%mat,1)     , & ! Notice here (ib2,ib1).
                                      czero                                                                         , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)     , size(tmpMat%blocks(ib1,ib2)%mat,1)        );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                           , &
                                      +cone                                                                         , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)     , size(tmpMat%blocks(ib1,ib2)%mat,1)      , &
                                      u(it)%blocks(ib2,ib2)%mat(1,1)      , size(u(it)%blocks(ib2,ib2)%mat,1)       , &
                                      cone                                                                          , &
                                      dh11_2(it)%blocks(ib1,ib2)%mat(1,1) , size(dh11_2(it)%blocks(ib1,ib2)%mat,1)    );

                      end if

                      if( dDelta1_pl(it)%nnzblocks(ib1,ib2) ) then

                          call zgemm( 'c' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                                  , &
                                      +cone                                                                                , &
                                      v(it)%blocks(ib1,ib1)%mat(1,1)          , size(v(it)%blocks(ib1,ib1)%mat,1)          , &
                                      dDelta1_pl(it)%blocks(ib1,ib2)%mat(1,1) , size(dDelta1_pl(it)%blocks(ib1,ib2)%mat,1) , &
                                      czero                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)           );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                                  , &
                                      +cone                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)         , &
                                      u(it)%blocks(ib2,ib2)%mat(1,1)          , size(u(it)%blocks(ib2,ib2)%mat,1)          , &
                                      cone                                                                                 , &
                                      dh11_2(it)%blocks(ib1,ib2)%mat(1,1)     , size(dh11_2(it)%blocks(ib1,ib2)%mat,1)       );

                      end if

                      if( dDelta1_mi(it)%nnzblocks(ib2,ib1) ) then ! Notice here (ib2,ib1).

                          call zgemm( 'c' , 'c' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                                  , &
                                      +cone                                                                                , &
                                      u(it)%blocks(ib1,ib1)%mat(1,1)          , size(u(it)%blocks(ib1,ib1)%mat,1)          , &
                                      dDelta1_mi(it)%blocks(ib2,ib1)%mat(1,1) , size(dDelta1_mi(it)%blocks(ib2,ib1)%mat,1) , & ! Notice here (ib2,ib1).
                                      czero                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)           );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                                  , &
                                      +cone                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)         , &
                                      v(it)%blocks(ib2,ib2)%mat(1,1)          , size(v(it)%blocks(ib2,ib2)%mat,1)          , &
                                      cone                                                                                 , &
                                      dh11_2(it)%blocks(ib1,ib2)%mat(1,1)     , size(dh11_2(it)%blocks(ib1,ib2)%mat,1)       );

                      end if

                  end if
              end do
          end do

          call transposeBlockMatrix( dh11_2(it) );

      end do




      return;
      end subroutine calculate_dh11
