!======================================================================!

      subroutine fam_drhodkappa()

!======================================================================!
      use simplex;
      use u_matrix;
      use v_matrix;
      use drho;
      use dkappa;
      use xyfam;
      use tempBlockMatrix;
      implicit none;

      double complex , parameter :: cone  = cmplx( 1.d0 , 0.d0 , kind=8 );
      double complex , parameter :: czero = cmplx( 0.d0 , 0.d0 , kind=8 );
      integer                    :: it, ib1, ib2;

      ! Induced density matrix and pairing tensor are given in Eqs. (B.8), (B.10) and (B.11):
      !
      !     drho1      = ( -  u * x * v'  -  v * y^T * u' )   ,
      !     dkappa1_pl = ( +  u * x * u'  -  v * y^T * v' )   ,
      !     drho2      = ( -  v * x * u'  -  u * y^T * v' )^T ,
      !     dkappa1_mi = ( -  v * x * v'  +  u * y^T * u' )'  .
      !
      ! Straightforward calculation gives drho1,drho2,dkappa1_pl,dkappa1_mi with 16 matrix
      ! multiplications. However, they can be calculated with 12 matrix multiplications.
      !
      ! First we calculate (u*x) and (v*y^T) and then:
      ! drho1      =   - (u*x)*v' - (v*y^T)*u'     ,
      ! dkappa1_pl =   + (u*x)*u' - (v*y^T)*v'     .
      !
      ! First we calculate (v*x) and (u*y^T) and then:
      ! drho2      = ( - (v*x)*u' - (u*y^T)*v' )^T ,
      ! dkappa1_mi = ( - (v*x)*v' + (u*y^T)*u' )'  .




      ! Calculation of drho1 and dkappa1_pl.
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  call assert( drho1(it)%nnzblocks(ib1,ib2) .eqv. dkappa1_pl(it)%nnzblocks(ib1,ib2) , 'drho1 and dkappa1_pl nnzblocks error!' );

                  if( drho1(it)%nnzblocks(ib1,ib2) .and. dkappa1_pl(it)%nnzblocks(ib1,ib2) ) then

                      drho1     (it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = czero;
                      dkappa1_pl(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = czero;

                      if( x_fam(it)%nnzblocks(ib1,ib2) ) then

                          ! (u*x).
                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                                  , &
                                      +cone                                                                                , &
                                      u(it)%blocks(ib1,ib1)%mat(1,1)          , size(u(it)%blocks(ib1,ib1)%mat,1)          , &
                                      x_fam(it)%blocks(ib1,ib2)%mat(1,1)      , size(x_fam(it)%blocks(ib1,ib2)%mat,1)      , &
                                      czero                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)           );

                          ! drho1 = drho1 - (u*x)*v'.
                          call zgemm( 'n' , 'c' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                                  , &
                                      -cone                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)         , &
                                      v(it)%blocks(ib2,ib2)%mat(1,1)          , size(v(it)%blocks(ib2,ib2)%mat,1)          , &
                                      cone                                                                                 , &
                                      drho1(it)%blocks(ib1,ib2)%mat(1,1)      , size(drho1(it)%blocks(ib1,ib2)%mat,1)        );

                          ! dkappa1_pl = dkappa1_pl + (u*x)*u'.
                          call zgemm( 'n' , 'c' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                                  , &
                                      +cone                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)         , &
                                      u(it)%blocks(ib2,ib2)%mat(1,1)          , size(u(it)%blocks(ib2,ib2)%mat,1)          , &
                                      cone                                                                                 , &
                                      dkappa1_pl(it)%blocks(ib1,ib2)%mat(1,1) , size(dkappa1_pl(it)%blocks(ib1,ib2)%mat,1)   );

                      end if

                      if( y_fam(it)%nnzblocks(ib2,ib1) ) then ! Notice here (ib2,ib1).

                          ! (v*y^T).
                          call zgemm( 'n' , 't' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                                  , &
                                      +cone                                                                                , &
                                      v(it)%blocks(ib1,ib1)%mat(1,1)          , size(v(it)%blocks(ib1,ib1)%mat,1)          , &
                                      y_fam(it)%blocks(ib2,ib1)%mat(1,1)      , size(y_fam(it)%blocks(ib2,ib1)%mat,1)      , & ! Notice here (ib2,ib1).
                                      czero                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)           );

                          ! drho1 = drho1 - (v*y^T)*u'.
                          call zgemm( 'n' , 'c' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                                  , &
                                      -cone                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)         , &
                                      u(it)%blocks(ib2,ib2)%mat(1,1)          , size(u(it)%blocks(ib2,ib2)%mat,1)          , &
                                      cone                                                                                 , &
                                      drho1(it)%blocks(ib1,ib2)%mat(1,1)      , size(drho1(it)%blocks(ib1,ib2)%mat,1)        );

                          ! dkappa1_pl = dkappa1_pl - (v*y^T)*v'.
                          call zgemm( 'n' , 'c' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                                  , &
                                      -cone                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)         , &
                                      v(it)%blocks(ib2,ib2)%mat(1,1)          , size(v(it)%blocks(ib2,ib2)%mat,1)          , &
                                      cone                                                                                 , &
                                      dkappa1_pl(it)%blocks(ib1,ib2)%mat(1,1) , size(dkappa1_pl(it)%blocks(ib1,ib2)%mat,1)   );

                      end if

                  end if

              end do
          end do
      end do




      ! Calculation of drho2 and dkappa1_mi.
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  call assert( drho2(it)%nnzblocks(ib1,ib2) .eqv. dkappa1_mi(it)%nnzblocks(ib1,ib2) , 'drho2 and dkappa1_mi nnzblocks error!' );

                  if( drho2(it)%nnzblocks(ib1,ib2) .and. dkappa1_mi(it)%nnzblocks(ib1,ib2) ) then

                      drho2     (it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = czero;
                      dkappa1_mi(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = czero;

                      if( x_fam(it)%nnzblocks(ib1,ib2) ) then

                          ! (v*x).
                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                                  , &
                                      +cone                                                                                , &
                                      v(it)%blocks(ib1,ib1)%mat(1,1)          , size(v(it)%blocks(ib1,ib1)%mat,1)          , &
                                      x_fam(it)%blocks(ib1,ib2)%mat(1,1)      , size(x_fam(it)%blocks(ib1,ib2)%mat,1)      , &
                                      czero                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)           );

                          ! drho2 = drho2 - (v*x)*u'.
                          call zgemm( 'n' , 'c' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                                  , &
                                      -cone                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)         , &
                                      u(it)%blocks(ib2,ib2)%mat(1,1)          , size(u(it)%blocks(ib2,ib2)%mat,1)          , &
                                      cone                                                                                 , &
                                      drho2(it)%blocks(ib1,ib2)%mat(1,1)      , size(drho2(it)%blocks(ib1,ib2)%mat,1)        );

                          ! dkappa1_mi = dkappa1_mi - (v*x)*v'.
                          call zgemm( 'n' , 'c' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                                  , &
                                      -cone                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)         , &
                                      v(it)%blocks(ib2,ib2)%mat(1,1)          , size(v(it)%blocks(ib2,ib2)%mat,1)          , &
                                      cone                                                                                 , &
                                      dkappa1_mi(it)%blocks(ib1,ib2)%mat(1,1) , size(dkappa1_mi(it)%blocks(ib1,ib2)%mat,1)   );

                      end if

                      if( y_fam(it)%nnzblocks(ib2,ib1) ) then ! Notice here (ib2,ib1).

                          ! (u*y^T).
                          call zgemm( 'n' , 't' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                                  , &
                                      +cone                                                                                , &
                                      u(it)%blocks(ib1,ib1)%mat(1,1)          , size(u(it)%blocks(ib1,ib1)%mat,1)          , &
                                      y_fam(it)%blocks(ib2,ib1)%mat(1,1)      , size(y_fam(it)%blocks(ib2,ib1)%mat,1)      , & ! Notice here (ib2,ib1).
                                      czero                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)           );

                          ! drho2 = drho2 - (u*y^T)*v'.
                          call zgemm( 'n' , 'c' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                                  , &
                                      -cone                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)         , &
                                      v(it)%blocks(ib2,ib2)%mat(1,1)          , size(v(it)%blocks(ib2,ib2)%mat,1)          , &
                                      cone                                                                                 , &
                                      drho2(it)%blocks(ib1,ib2)%mat(1,1)      , size(drho2(it)%blocks(ib1,ib2)%mat,1)        );

                          ! dkappa1_mi = dkappa1_mi + (u*y^T)*u'.
                          call zgemm( 'n' , 'c' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                                  , &
                                      +cone                                                                                , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)         , size(tmpMat%blocks(ib1,ib2)%mat,1)         , &
                                      u(it)%blocks(ib2,ib2)%mat(1,1)          , size(u(it)%blocks(ib2,ib2)%mat,1)          , &
                                      cone                                                                                 , &
                                      dkappa1_mi(it)%blocks(ib1,ib2)%mat(1,1) , size(dkappa1_mi(it)%blocks(ib1,ib2)%mat,1)   );

                      end if

                  end if

              end do
          end do

          call transposeBlockMatrix( drho2(it) );
          call hermitianTransposeBlockMatrix( dkappa1_mi(it) );

      end do




      return;
      end subroutine fam_drhodkappa
