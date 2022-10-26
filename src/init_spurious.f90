!======================================================================!

      subroutine init_spurious()

!======================================================================!
      use dataTypes;
      use fam_input;
      use simplex;
      use u_matrix;
      use v_matrix;
      use quadrature;
      use basis;
      use spurious;
      use tempBlockMatrix;
      implicit none;

      double precision             :: r, z;
      integer                      :: ih, il;
      double complex   , parameter :: cone  = cmplx( 1.d0 , 0.d0 , kind=8 );
      double complex   , parameter :: czero = cmplx( 0.d0 , 0.d0 , kind=8 );
      double complex   , parameter :: II    = cmplx( 0.d0 , 1.d0 , kind=8 );
      integer                      :: it, ib1, ib2, i1, i2;
      character                    :: fg1, fg2;
      integer                      :: nz1, nz2;
      integer                      :: nr1, nr2;
      integer                      :: ml1, ml2;
      double precision             :: integral;
      double complex               :: trace;
      type(complexBlockMatrix)     :: rcm;
      type(complexBlockMatrix)     :: pcm;

!--------------------------------------------------------------!
!     Rcm operator has the following form:                     !
!     Rcm = z/A for K=0, and Rcm = x/A for K=1.                !
!                                                              !
!     Pcm operator has the following form:                     !
!     Pcm = -i * d/dz for K=0, and Pcm = -i * d/dx for K=1.    !
!--------------------------------------------------------------!

      ! Translational spurious mode removal only if J is odd and K=0 or K=1.
      if( .not.( mod(J_multipole,2)==1 .and. any(K_multipole==[0,1]) ) ) then
          return;
      end if




      ! Calculation of rcm.
      allocate( rcm%nnzblocks(1:N_blocks,1:N_blocks) );
      allocate( rcm%blocks   (1:N_blocks,1:N_blocks) );
      rcm%nnzblocks(:,:) = .false.;
      do ib1 = 1 , N_blocks
          do ib2 = 1 , N_blocks
              if( r20(1)%nnzblocks(ib1,ib2) .or. r20(2)%nnzblocks(ib1,ib2) ) then
                  rcm%nnzblocks(ib1,ib2) = .true.;
                  allocate( rcm%blocks(ib1,ib2)%mat(1:id_spx(ib1),1:id_spx(ib2)) );
              end if
          end do
      end do

      if( K_multipole == 0 ) then
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( rcm%nnzblocks(ib1,ib2) ) then

                      rcm%blocks(ib1,ib2)%mat(:,:) = czero;

                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , id_spx(ib1)

                              fg1 = fg_spx(ib1)%index(i1);
                              nz1 = nz_spx(ib1)%index(i1);
                              nr1 = nr_spx(ib1)%index(i1);
                              ml1 = ml_spx(ib1)%index(i1);

                              fg2 = fg_spx(ib2)%index(i2);
                              nz2 = nz_spx(ib2)%index(i2);
                              nr2 = nr_spx(ib2)%index(i2);
                              ml2 = ml_spx(ib2)%index(i2);

                              if( fg1==fg2 .and. nr1==nr2 .and. abs(ml1-ml2)==K_multipole ) then

                                  integral = 0.d0;
                                  do ih = -NGH , +NGH
                                      if( ih == 0 ) cycle;
                                      z = sign(1,ih) * zb_fam(abs(ih));

                                      integral = integral + wz(abs(ih)) * z * phi_z(ib1)%mat(ih,i1) * phi_z(ib2)%mat(ih,i2);
                                  end do

                                  rcm%blocks(ib1,ib2)%mat(i1,i2) = 1.d0/(nucleusZ+nucleusN) * integral;

                              end if

                          end do
                      end do
                  end if
              end do
          end do
      end if

      if( K_multipole == 1 ) then
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( rcm%nnzblocks(ib1,ib2) ) then

                      rcm%blocks(ib1,ib2)%mat(:,:) = czero;

                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , id_spx(ib1)

                              fg1 = fg_spx(ib1)%index(i1);
                              nz1 = nz_spx(ib1)%index(i1);
                              nr1 = nr_spx(ib1)%index(i1);
                              ml1 = ml_spx(ib1)%index(i1);

                              fg2 = fg_spx(ib2)%index(i2);
                              nz2 = nz_spx(ib2)%index(i2);
                              nr2 = nr_spx(ib2)%index(i2);
                              ml2 = ml_spx(ib2)%index(i2);

                              if( fg1==fg2 .and. nz1==nz2 .and. abs(ml1-ml2)==K_multipole ) then

                                  integral = 0.d0;
                                  do il = 1 , NGL
                                      r = rb_fam(il);

                                      integral = integral + wr(il) * r * phi_r(ib1)%mat(il,i1) * phi_r(ib2)%mat(il,i2);
                                  end do

                                  rcm%blocks(ib1,ib2)%mat(i1,i2) = 0.5d0/(nucleusZ+nucleusN) * integral;

                              end if

                          end do
                      end do
                  end if
              end do
          end do
      end if




      ! Calculation of pcm.
      allocate( pcm%nnzblocks(1:N_blocks,1:N_blocks) );
      allocate( pcm%blocks   (1:N_blocks,1:N_blocks) );
      pcm%nnzblocks(:,:) = .false.;
      do ib1 = 1 , N_blocks
          do ib2 = 1 , N_blocks
              if( p20(1)%nnzblocks(ib1,ib2) .or. p20(2)%nnzblocks(ib1,ib2) ) then
                  pcm%nnzblocks(ib1,ib2) = .true.;
                  allocate( pcm%blocks(ib1,ib2)%mat(1:id_spx(ib1),1:id_spx(ib2)) );
              end if
          end do
      end do

      if( K_multipole == 0 ) then
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( pcm%nnzblocks(ib1,ib2) ) then

                      pcm%blocks(ib1,ib2)%mat(:,:) = czero;

                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , id_spx(ib1)

                              fg1 = fg_spx(ib1)%index(i1);
                              nz1 = nz_spx(ib1)%index(i1);
                              nr1 = nr_spx(ib1)%index(i1);
                              ml1 = ml_spx(ib1)%index(i1);

                              fg2 = fg_spx(ib2)%index(i2);
                              nz2 = nz_spx(ib2)%index(i2);
                              nr2 = nr_spx(ib2)%index(i2);
                              ml2 = ml_spx(ib2)%index(i2);

                              if( fg1==fg2 .and. nr1==nr2 .and. abs(ml1-ml2)==K_multipole ) then

                                  integral = 0.d0;
                                  do ih = -NGH , +NGH
                                      if( ih == 0 ) cycle;

                                      integral = integral + wz(abs(ih)) * ( + phi_z(ib1)%mat(ih,i1)*dphi_z(ib2)%mat(ih,i2) &
                                                                            - phi_z(ib2)%mat(ih,i2)*dphi_z(ib1)%mat(ih,i1) );
                                  end do

                                  pcm%blocks(ib1,ib2)%mat(i1,i2) = -II/2 * integral;

                              end if

                          end do
                      end do
                  end if
              end do
          end do
      end if

      if( K_multipole == 1 ) then
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( pcm%nnzblocks(ib1,ib2) ) then

                      pcm%blocks(ib1,ib2)%mat(:,:) = czero;

                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , id_spx(ib1)

                              fg1 = fg_spx(ib1)%index(i1);
                              nz1 = nz_spx(ib1)%index(i1);
                              nr1 = nr_spx(ib1)%index(i1);
                              ml1 = ml_spx(ib1)%index(i1);

                              fg2 = fg_spx(ib2)%index(i2);
                              nz2 = nz_spx(ib2)%index(i2);
                              nr2 = nr_spx(ib2)%index(i2);
                              ml2 = ml_spx(ib2)%index(i2);

                              if( fg1==fg2 .and. nz1==nz2 .and. abs(ml1-ml2)==K_multipole ) then

                                  integral = 0.d0;
                                  do il = 1 , NGL
                                      r = rb_fam(il);

                                      integral = integral + wr(il) * ( + phi_r(ib1)%mat(il,i1)*dphi_r(ib2)%mat(il,i2) &
                                                                       - phi_r(ib2)%mat(il,i2)*dphi_r(ib1)%mat(il,i1) );

                                      integral = integral + wr(il) * (ml2**2-ml1**2)/r * phi_r(ib1)%mat(il,i1) * phi_r(ib2)%mat(il,i2);
                                  end do

                                  pcm%blocks(ib1,ib2)%mat(i1,i2) = -II/4 * integral;

                              end if

                          end do
                      end do
                  end if
              end do
          end do
      end if




      ! Calculation of r20 = - u' * rcm   * v
      !                      - v' * rcm^T * u.
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( r20(it)%nnzblocks(ib1,ib2) ) then

                      r20(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = czero;

                      if( rcm%nnzblocks(ib1,ib2) ) then

                          call zgemm( 'c' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                        , &
                                      -cone                                                                      , &
                                      u(it)%blocks(ib1,ib1)%mat(1,1)     , size(u(it)%blocks(ib1,ib1)%mat,1)     , &
                                      rcm%blocks(ib1,ib2)%mat(1,1)       , size(rcm%blocks(ib1,ib2)%mat,1)       , &
                                      czero                                                                      , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)    , size(tmpMat%blocks(ib1,ib2)%mat,1)      );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                        , &
                                      +cone                                                                      , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)    , size(tmpMat%blocks(ib1,ib2)%mat,1)    , &
                                      v(it)%blocks(ib2,ib2)%mat(1,1)     , size(v(it)%blocks(ib2,ib2)%mat,1)     , &
                                      cone                                                                       , &
                                      r20(it)%blocks(ib1,ib2)%mat(1,1)   , size(r20(it)%blocks(ib1,ib2)%mat,1)     );

                      end if

                      if( rcm%nnzblocks(ib2,ib1) ) then ! Notice here (ib2,ib1).

                          call zgemm( 'c' , 't' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                        , &
                                      -cone                                                                      , &
                                      v(it)%blocks(ib1,ib1)%mat(1,1)     , size(v(it)%blocks(ib1,ib1)%mat,1)     , &
                                      rcm%blocks(ib2,ib1)%mat(1,1)       , size(rcm%blocks(ib2,ib1)%mat,1)       , & ! Notice here (ib2,ib1).
                                      czero                                                                      , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)    , size(tmpMat%blocks(ib1,ib2)%mat,1)      );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                        , &
                                      +cone                                                                      , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)    , size(tmpMat%blocks(ib1,ib2)%mat,1)    , &
                                      u(it)%blocks(ib2,ib2)%mat(1,1)     , size(u(it)%blocks(ib2,ib2)%mat,1)     , &
                                      cone                                                                       , &
                                      r20(it)%blocks(ib1,ib2)%mat(1,1)   , size(r20(it)%blocks(ib1,ib2)%mat,1)     );

                      end if

                  end if
              end do
          end do
      end do


      ! Calculation of p20 = - u' * pcm   * v
      !                      - v' * pcm^T * u.
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( p20(it)%nnzblocks(ib1,ib2) ) then

                      p20(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = czero;

                      if( pcm%nnzblocks(ib1,ib2) ) then

                          call zgemm( 'c' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                        , &
                                      -cone                                                                      , &
                                      u(it)%blocks(ib1,ib1)%mat(1,1)     , size(u(it)%blocks(ib1,ib1)%mat,1)     , &
                                      pcm%blocks(ib1,ib2)%mat(1,1)       , size(pcm%blocks(ib1,ib2)%mat,1)       , &
                                      czero                                                                      , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)    , size(tmpMat%blocks(ib1,ib2)%mat,1)      );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                        , &
                                      +cone                                                                      , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)    , size(tmpMat%blocks(ib1,ib2)%mat,1)    , &
                                      v(it)%blocks(ib2,ib2)%mat(1,1)     , size(v(it)%blocks(ib2,ib2)%mat,1)     , &
                                      cone                                                                       , &
                                      p20(it)%blocks(ib1,ib2)%mat(1,1)   , size(p20(it)%blocks(ib1,ib2)%mat,1)     );

                      end if

                      if( pcm%nnzblocks(ib2,ib1) ) then ! Notice here (ib2,ib1).

                          call zgemm( 'c' , 't' , id_spx(ib1) , id_spx(ib2) , id_spx(ib1)                        , &
                                      -cone                                                                      , &
                                      v(it)%blocks(ib1,ib1)%mat(1,1)     , size(v(it)%blocks(ib1,ib1)%mat,1)     , &
                                      pcm%blocks(ib2,ib1)%mat(1,1)       , size(pcm%blocks(ib2,ib1)%mat,1)       , & ! Notice here (ib2,ib1).
                                      czero                                                                      , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)    , size(tmpMat%blocks(ib1,ib2)%mat,1)      );

                          call zgemm( 'n' , 'n' , id_spx(ib1) , id_spx(ib2) , id_spx(ib2)                        , &
                                      +cone                                                                      , &
                                      tmpMat%blocks(ib1,ib2)%mat(1,1)    , size(tmpMat%blocks(ib1,ib2)%mat,1)    , &
                                      u(it)%blocks(ib2,ib2)%mat(1,1)     , size(u(it)%blocks(ib2,ib2)%mat,1)     , &
                                      cone                                                                       , &
                                      p20(it)%blocks(ib1,ib2)%mat(1,1)   , size(p20(it)%blocks(ib1,ib2)%mat,1)     );

                      end if

                  end if
              end do
          end do
      end do




      ! Calculation of <Phi| [Rcm,Pcm] |Phi> = 2*Tr[ r20' * p20 ].
      do it = 1 , 2

          trace = czero;
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( r20(it)%nnzblocks(ib1,ib2) .and. p20(it)%nnzblocks(ib1,ib2) ) then
                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , id_spx(ib1)
                              trace = trace + conjg(r20(it)%blocks(ib1,ib2)%mat(i1,i2)) * p20(it)%blocks(ib1,ib2)%mat(i1,i2);
                          end do
                      end do
                  end if
              end do
          end do
          RcmPcm_commutator(it) = 2*trace;

      end do

#ifdef DEBUG
      write(*,'(a,f20.12,a,f20.12,a)') '[Rcm,Pcm] = ' , real(sum(RcmPcm_commutator(1:2))) , ' + ' , imag(sum(RcmPcm_commutator(1:2))) , ' i.';
#endif




      ! Deallocate rcm and pcm.
      do ib2 = 1 , size(rcm%nnzblocks,2)
          do ib1 = 1 , size(rcm%nnzblocks,1)
              if( rcm%nnzblocks(ib1,ib2) ) deallocate( rcm%blocks(ib1,ib2)%mat );
          end do
      end do
      deallocate( rcm%blocks    );
      deallocate( rcm%nnzblocks );

      do ib2 = 1 , size(pcm%nnzblocks,2)
          do ib1 = 1 , size(pcm%nnzblocks,1)
              if( pcm%nnzblocks(ib1,ib2) ) deallocate( pcm%blocks(ib1,ib2)%mat );
          end do
      end do
      deallocate( pcm%blocks    );
      deallocate( pcm%nnzblocks );




      return;
      end subroutine init_spurious
