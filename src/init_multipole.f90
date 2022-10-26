!======================================================================!

      subroutine init_multipole()

!======================================================================!
      use dataTypes;
      use fam_input;
      use simplex;
      use u_matrix;
      use v_matrix;
      use quadrature;
      use wbasis;
      use fmatrix;
      use f02f20matrix;
      use tempBlockMatrix;
      implicit none;

      double precision , external  :: tau3;
      double precision             :: f( -NGH:+NGH , 1:NGL , 1:2 , 0:J_MAX , 0:J_MAX );
      double precision             :: fac_iso(1:2);
      double precision             :: r, z;
      integer                      :: ih, il;
      double precision , parameter :: pi    = 4.d0*atan(1.d0);
      double complex   , parameter :: cone  = cmplx( 1.d0 , 0.d0 , kind=8 );
      double complex   , parameter :: czero = cmplx( 0.d0 , 0.d0 , kind=8 );
      integer                      :: it, ib1, ib2, i1, i2;
      integer                      :: J, K;
      character                    :: fg1, fg2;
      integer                      :: nz1, nz2;
      integer                      :: ml1, ml2;
      double precision             :: integral;

!-------------------------------------------------------------------------------!
!     Excitation operator has the following form:                               !
!                                                                               !
!     f = fac_iso(it) * 1/sqrt(2+2*delta_{K,0})                                 !
!                     * |R|^J                                                   !
!                     * ( Y_{J,K}(theta,phi) + (-1)^K * Y_{J,-K}(theta,phi) ),  !
!                                                                               !
!     where (R,theta,phi) are spherical coordinates.                            !
!                                                                               !
!     Only (J=0,K=0) is an exception: f = fac_iso(it) * |R|^2.                  !
!                                                                               !
!     We coded it=1 for neutrons and it=2 for protons.                          !
!                                                                               !
!     For isoscalar excitation:                                                 !
!         fac_iso(it) = +1.                                                     !
!                                                                               !
!     For isovector excitation:                                                 !
!         fac_iso(it) = tau3(it).                                               !
!                                                                               !
!     with only exception for isovector J=1 excitation, where we use:           !
!         fac_iso(it=1) = tau3(it=1) * Z/(Z+N) (neutrons),                      !
!         fac_iso(it=2) = tau3(it=2) * N/(Z+N) (protons).                       !
!                                                                               !
!     tau3(it) is equal to +1 or -1 depending on the isospin convention.        !
!-------------------------------------------------------------------------------!

      call assert( J_multipole <= J_MAX       , 'J >  J_MAX            ' );
      call assert( J_multipole >= 0           , 'J <  0                ' );
      call assert( K_multipole <= J_multipole , 'K >  J                ' );
      call assert( K_multipole >= 0           , 'K <  0                ' );
      call assert( any(Isospin==[0,1])        , 'Isospin must be 0 or 1' );




      ! Initializing multipole operators in coordinate space.
      if( Isospin == 0 ) then
          do it = 1 , 2
              fac_iso(it) = +1.d0;
          end do
      end if
      if( Isospin == 1 ) then
          do it = 1 , 2
              fac_iso(it) = tau3(it);
          end do
          if( J_multipole == 1 ) then
              fac_iso(1) = ( tau3(1) * nucleusZ ) / ( nucleusZ + nucleusN );
              fac_iso(2) = ( tau3(2) * nucleusN ) / ( nucleusZ + nucleusN );
          end if
      end if

      do it = 1 , 2
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;

                  ! Cylindrical coordinates.
                  r = rb_fam(il);
                  z = sign(1,ih) * zb_fam(abs(ih));


                  ! J=0, K=0.
                  f(ih,il,it,0,0) = fac_iso(it) * ( + 1               ) * ( z**2 + r**2 );                         ! cos(0*phi).

                  ! J=1, K=0,1.
                  f(ih,il,it,1,0) = fac_iso(it) * ( +sqrt(   3/pi)/ 2 ) * z;                                       ! cos(0*phi).
                  f(ih,il,it,1,1) = fac_iso(it) * ( -sqrt(   3/pi)/ 2 ) * r;                                       ! cos(1*phi).

                  ! J=2, K=0,1,2
                  f(ih,il,it,2,0) = fac_iso(it) * ( +sqrt(   5/pi)/ 4 ) * ( 2*z**2 - r**2 );                       ! cos(0*phi).
                  f(ih,il,it,2,1) = fac_iso(it) * ( -sqrt(  15/pi)/ 2 ) * r**1 * z;                                ! cos(1*phi).
                  f(ih,il,it,2,2) = fac_iso(it) * ( +sqrt(  15/pi)/ 4 ) * r**2;                                    ! cos(2*phi).

                  ! J=3, K=0,1,2,3.
                  f(ih,il,it,3,0) = fac_iso(it) * ( +sqrt(   7/pi)/ 4 ) * z * ( 2*z**2 - 3*r**2 );                 ! cos(0*phi).
                  f(ih,il,it,3,1) = fac_iso(it) * ( -sqrt(  42/pi)/ 8 ) * r**1 * ( 4*z**2 - r**2 );                ! cos(1*phi).
                  f(ih,il,it,3,2) = fac_iso(it) * ( +sqrt( 105/pi)/ 4 ) * r**2 * z;                                ! cos(2*phi).
                  f(ih,il,it,3,3) = fac_iso(it) * ( -sqrt(  70/pi)/ 8 ) * r**3;                                    ! cos(3*phi).

                  ! J=4, K=0,1,2,3,4.
                  f(ih,il,it,4,0) = fac_iso(it) * ( +sqrt(   9/pi)/16 ) * ( 8*z**4 - 24*z**2*r**2 + 3*r**4 );      ! cos(0*phi).
                  f(ih,il,it,4,1) = fac_iso(it) * ( -sqrt(  90/pi)/ 8 ) * r**1 * z * ( 4*z**2 - 3*r**2 );          ! cos(1*phi).
                  f(ih,il,it,4,2) = fac_iso(it) * ( +sqrt(  45/pi)/ 8 ) * r**2 * ( 6*z**2 - r**2 );                ! cos(2*phi).
                  f(ih,il,it,4,3) = fac_iso(it) * ( -sqrt( 630/pi)/ 8 ) * r**3 * z;                                ! cos(3*phi).
                  f(ih,il,it,4,4) = fac_iso(it) * ( +sqrt( 315/pi)/16 ) * r**4;                                    ! cos(4*phi).

                  ! J=5, K=0,1,2,3,4,5.
                  f(ih,il,it,5,0) = fac_iso(it) * ( +sqrt(  11/pi)/16 ) * z * ( 8*z**4 - 40*z**2*r**2 + 15*r**4 ); ! cos(0*phi).
                  f(ih,il,it,5,1) = fac_iso(it) * ( -sqrt( 165/pi)/16 ) * r**1 * ( 8*z**4 - 12*z**2*r**2 + r**4 ); ! cos(1*phi).
                  f(ih,il,it,5,2) = fac_iso(it) * ( +sqrt(1155/pi)/ 8 ) * r**2 * z * ( 2*z**2 - r**2 );            ! cos(2*phi).
                  f(ih,il,it,5,3) = fac_iso(it) * ( -sqrt( 770/pi)/32 ) * r**3 * ( 8*z**2 - r**2 );                ! cos(3*phi).
                  f(ih,il,it,5,4) = fac_iso(it) * ( +sqrt(3465/pi)/16 ) * r**4 * z;                                ! cos(4*phi).
                  f(ih,il,it,5,5) = fac_iso(it) * ( -sqrt(1386/pi)/32 ) * r**5;                                    ! cos(5*phi).

              end do
          end do
      end do




!-----Calculation of f1_JK matrix.
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  if( f1_JK(it)%nnzblocks(ib1,ib2) ) f1_JK(it)%blocks(ib1,ib2)%mat(:,:) = czero;

                  do i2 = 1 , id_spx(ib2)
                      do i1 = 1 , id_spx(ib1)

                          fg1 = fg_spx(ib1)%index(i1);
                          nz1 = nz_spx(ib1)%index(i1);
                          ml1 = ml_spx(ib1)%index(i1);

                          fg2 = fg_spx(ib2)%index(i2);
                          nz2 = nz_spx(ib2)%index(i2);
                          ml2 = ml_spx(ib2)%index(i2);

                          if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then

                              J = J_multipole;
                              K = K_multipole;

                              integral = 0.d0;
                              do il = 1 , NGL
                                  do ih = 1 , NGH
                                      integral = integral + ( f(+ih,il,it,J,K) + merge(+1.d0,-1.d0,mod(nz1+nz2,2)==0)*f(-ih,il,it,J,K) ) * wPhi(ib1)%arr(ih,il,i1) * wPhi(ib2)%arr(ih,il,i2);
                                  end do
                              end do
                              integral = merge( 1.d0 , 0.5d0 , K_multipole==0 ) * integral;

                              f1_JK(it)%blocks(ib1,ib2)%mat(i1,i2) = cmplx( integral , 0.d0 , kind=8 );

                          end if

                      end do
                  end do
              end do
          end do
      end do

!-----Calculation of f2_JK matrix.
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks

                  if( f2_JK(it)%nnzblocks(ib1,ib2) ) f2_JK(it)%blocks(ib1,ib2)%mat(:,:) = czero;

                  do i2 = 1 , id_spx(ib2)
                      do i1 = 1 , id_spx(ib1)

                          fg1 = fg_spx(ib1)%index(i1);
                          nz1 = nz_spx(ib1)%index(i1);
                          ml1 = ml_spx(ib1)%index(i1);

                          fg2 = fg_spx(ib2)%index(i2);
                          nz2 = nz_spx(ib2)%index(i2);
                          ml2 = ml_spx(ib2)%index(i2);

                          if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then

                              J = J_multipole;
                              K = K_multipole;

                              integral = 0.d0;
                              do il = 1 , NGL
                                  do ih = 1 , NGH
                                      integral = integral + ( f(+ih,il,it,J,K) + merge(+1.d0,-1.d0,mod(nz1+nz2,2)==0)*f(-ih,il,it,J,K) ) * wPhi(ib1)%arr(ih,il,i1) * wPhi(ib2)%arr(ih,il,i2);
                                  end do
                              end do
                              integral = merge( 1.d0 , 0.5d0 , K_multipole==0 ) * integral;

                              f2_JK(it)%blocks(ib1,ib2)%mat(i1,i2) = cmplx( integral , 0.d0 , kind=8 );

                          end if

                      end do
                  end do
              end do
          end do
      end do




      ! Calculation of f20 = - u' * f1_JK   * v
      !                      - v' * f2_JK^T * u , see Eq. (B.6).
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( f20(it)%nnzblocks(ib1,ib2) ) then

                      f20(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = czero;

                      if( f1_JK(it)%nnzblocks(ib1,ib2) ) then

                          call zgemm( 'c' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                        , &
                                      -cone                                                                      , &
                                      u(it)%blocks(ib1,ib1)%mat(1,1)     , size(u(it)%blocks(ib1,ib1)%mat,1)     , &
                                      f1_JK(it)%blocks(ib1,ib2)%mat(1,1) , size(f1_JK(it)%blocks(ib1,ib2)%mat,1) , &
                                      czero                                                                      , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)    , size(tmpMat%blocks(ib1,ib2)%mat,1)      );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                        , &
                                      +cone                                                                      , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)    , size(tmpMat%blocks(ib1,ib2)%mat,1)    , &
                                      v(it)%blocks(ib2,ib2)%mat(1,1)     , size(v(it)%blocks(ib2,ib2)%mat,1)     , &
                                      cone                                                                       , &
                                      f20(it)%blocks(ib1,ib2)%mat(1,1)   , size(f20(it)%blocks(ib1,ib2)%mat,1)     );

                      end if

                      if( f2_JK(it)%nnzblocks(ib2,ib1) ) then ! Notice here (ib2,ib1).

                          call zgemm( 'c' , 't' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                        , &
                                      -cone                                                                      , &
                                      v(it)%blocks(ib1,ib1)%mat(1,1)     , size(v(it)%blocks(ib1,ib1)%mat,1)     , &
                                      f2_JK(it)%blocks(ib2,ib1)%mat(1,1) , size(f2_JK(it)%blocks(ib2,ib1)%mat,1) , & ! Notice here (ib2,ib1).
                                      czero                                                                      , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)    , size(tmpMat%blocks(ib1,ib2)%mat,1)      );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                        , &
                                      +cone                                                                      , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)    , size(tmpMat%blocks(ib1,ib2)%mat,1)    , &
                                      u(it)%blocks(ib2,ib2)%mat(1,1)     , size(u(it)%blocks(ib2,ib2)%mat,1)     , &
                                      cone                                                                       , &
                                      f20(it)%blocks(ib1,ib2)%mat(1,1)   , size(f20(it)%blocks(ib1,ib2)%mat,1)     );

                      end if

                  end if
              end do
          end do
      end do

      ! Calculation of f02 = ( - v' * f1_JK   * u
      !                        - u' * f2_JK^T * v  )^T , see Eq. (B.6).
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( f02(it)%nnzblocks(ib1,ib2) ) then

                      f02(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = czero;

                      if( f1_JK(it)%nnzblocks(ib1,ib2) ) then

                          call zgemm( 'c' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                        , &
                                      -cone                                                                      , &
                                      v(it)%blocks(ib1,ib1)%mat(1,1)     , size(v(it)%blocks(ib1,ib1)%mat,1)     , &
                                      f1_JK(it)%blocks(ib1,ib2)%mat(1,1) , size(f1_JK(it)%blocks(ib1,ib2)%mat,1) , &
                                      czero                                                                      , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)    , size(tmpMat%blocks(ib1,ib2)%mat,1)      );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                        , &
                                      +cone                                                                      , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)    , size(tmpMat%blocks(ib1,ib2)%mat,1)    , &
                                      u(it)%blocks(ib2,ib2)%mat(1,1)     , size(u(it)%blocks(ib2,ib2)%mat,1)     , &
                                      cone                                                                       , &
                                      f02(it)%blocks(ib1,ib2)%mat(1,1)   , size(f02(it)%blocks(ib1,ib2)%mat,1)     );

                      end if

                      if( f2_JK(it)%nnzblocks(ib2,ib1) ) then ! Notice here (ib2,ib1).

                          call zgemm( 'c' , 't' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                        , &
                                      -cone                                                                      , &
                                      u(it)%blocks(ib1,ib1)%mat(1,1)     , size(u(it)%blocks(ib1,ib1)%mat,1)     , &
                                      f2_JK(it)%blocks(ib2,ib1)%mat(1,1) , size(f2_JK(it)%blocks(ib2,ib1)%mat,1) , & ! Notice here (ib2,ib1).
                                      czero                                                                      , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)    , size(tmpMat%blocks(ib1,ib2)%mat,1)      );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                        , &
                                      +cone                                                                      , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)    , size(tmpMat%blocks(ib1,ib2)%mat,1)    , &
                                      v(it)%blocks(ib2,ib2)%mat(1,1)     , size(v(it)%blocks(ib2,ib2)%mat,1)     , &
                                      cone                                                                       , &
                                      f02(it)%blocks(ib1,ib2)%mat(1,1)   , size(f02(it)%blocks(ib1,ib2)%mat,1)     );

                      end if

                  end if
              end do
          end do

          call transposeBlockMatrix( f02(it) );

      end do




      return;
      end subroutine init_multipole
