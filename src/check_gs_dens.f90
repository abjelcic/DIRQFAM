!======================================================================!

      subroutine check_gs_dens()

!======================================================================!
      use dataTypes;
      use fam_input;
      use simplex;
      use v_matrix;
      use herpol;
      use lagpol;
      use gaucor;
      use gaussh;
      use gaussl;
      use dens;
      implicit none;

      integer                                 :: it, ib, il, ih, ihl, i, j;
      character                               :: fg1, fg2;
      integer                                 :: nz1, nz2;
      integer                                 :: nr1, nr2;
      integer                                 :: ml1, ml2;
      double precision                        :: qhql1, qhql2;
      double complex           , parameter    :: cone  = cmplx( 1.d0 , 0.d0 , kind=8 );
      double complex           , parameter    :: czero = cmplx( 0.d0 , 0.d0 , kind=8 );
      type(complexBlockMatrix) , dimension(2) :: rho1;
      integer                                 :: MG;
      double precision                        :: rho_v( (NGH+1)*(NGL+1) , 2 );
      double precision                        :: rho_s( (NGH+1)*(NGL+1) , 2 );
      double precision                        :: err1, err2;
      double precision                        :: rho_new, rho_old, abs_error, rel_error;


      ! Density matrix rho = [ rho1 , 0 ; 0 , rho1^T ] = [ v*v' , 0 ; 0 , (v*v')^T ].
      do it = 1 , 2
          allocate( rho1(it)%blocks(1:N_blocks,1:N_blocks) );
          do ib = 1 , N_blocks
              allocate( rho1(it)%blocks(ib,ib)%mat(1:id_spx(ib),1:id_spx(ib)) );

              call zgemm( 'n' , 'c' , id_spx(ib) , id_spx(ib) , id_spx(ib)                    , &
                          cone                                                                , &
                          v(it)%blocks(ib,ib)%mat(1,1)    , size(v(it)%blocks(ib,ib)%mat,1)   , &
                          v(it)%blocks(ib,ib)%mat(1,1)    , size(v(it)%blocks(ib,ib)%mat,1)   , &
                          czero                                                               , &
                          rho1(it)%blocks(ib,ib)%mat(1,1) , size(rho1(it)%blocks(ib,ib)%mat,1)  );
          end do
      end do


      ! Calculating ground-state densities.
      rho_v(:,:) = 0.d0;
      rho_s(:,:) = 0.d0;
      do it = 1 , 2
          do il = 0 , NGL
              do ih = 0 , NGH

                  ihl = 1+ih + il*(NGH+1);

                  do ib = 1 , N_blocks
                      do i = 1 , id_spx(ib)
                          do j = 1 , id_spx(ib)

                              fg1 = fg_spx(ib)%index(i);
                              nz1 = nz_spx(ib)%index(i);
                              nr1 = nr_spx(ib)%index(i);
                              ml1 = ml_spx(ib)%index(i);

                              fg2 = fg_spx(ib)%index(j);
                              nz2 = nz_spx(ib)%index(j);
                              nr2 = nr_spx(ib)%index(j);
                              ml2 = ml_spx(ib)%index(j);

                              if( fg1==fg2 .and. abs(ml1-ml2)==0 ) then

                                  qhql1 = qh(nz1,ih) * ql(nr1,abs(ml1),il);
                                  qhql2 = qh(nz2,ih) * ql(nr2,abs(ml2),il);

                                  if( all([fg1,fg2]==['f','f']) ) then
                                      rho_v(ihl,it) = rho_v(ihl,it) + 2*real(rho1(it)%blocks(ib,ib)%mat(i,j)) * qhql1 * qhql2;
                                      rho_s(ihl,it) = rho_s(ihl,it) + 2*real(rho1(it)%blocks(ib,ib)%mat(i,j)) * qhql1 * qhql2;
                                  end if
                                  if( all([fg1,fg2]==['g','g']) ) then
                                      rho_v(ihl,it) = rho_v(ihl,it) + 2*real(rho1(it)%blocks(ib,ib)%mat(i,j)) * qhql1 * qhql2;
                                      rho_s(ihl,it) = rho_s(ihl,it) - 2*real(rho1(it)%blocks(ib,ib)%mat(i,j)) * qhql1 * qhql2;
                                  end if

                              end if

                          end do
                      end do
                  end do

                  rho_v(ihl,it) = rho_v(ihl,it) / ww(ihl);
                  rho_s(ihl,it) = rho_s(ihl,it) / ww(ihl);

              end do
          end do
      end do


      ! Comparing densities.
      MG = (NGH+1)*(NGL+1);
      err1 = norm2( sum(rho_v(1:MG,1:2),2) - ro(1:MG,2) ) / norm2(ro(1:MG,2));
      err2 = norm2( sum(rho_s(1:MG,1:2),2) - ro(1:MG,1) ) / norm2(ro(1:MG,1));
      call assert( err1 < 1.d-8 , 'Ground state densities problem.' );
      call assert( err2 < 1.d-8 , 'Ground state densities problem.' );


#ifdef DEBUG
      write(6,'(a)') 'Comparison of the isoscalar-vector densities';
      write(6,'(a)') '';
      write(6,'(a)') '          coordinate [fm]                  rho_new           rho_old           |error|       rel.error';
      write(6,'(a)') '';
      do il = 0 , NGL
          do ih = 0 , NGH
              ihl = 1+ih + il*(NGH+1);

              rho_new   = rho_v(ihl,1) + rho_v(ihl,2);
              rho_old   = ro(ihl,2);
              abs_error = abs( rho_new - rho_old );
              rel_error = abs_error / abs(rho_old);


              write(6,'(a,i2,a,1f6.2,a,a,i2,a,1f6.2,a,2f18.13,2e15.5)') '( r(' , il , ') = ' , rb(il) , ',  '  , &
                                                                        '  z(' , ih , ') = ' , zb(ih) , ' )  ' , &
                                                                        rho_new                                , &
                                                                        rho_old                                , &
                                                                        abs_error                              , &
                                                                        rel_error                              ;
          end do
      end do


      write(6,'(a)') 'Comparison of the isoscalar-scalar densities';
      write(6,'(a)') '';
      write(6,'(a)') '        coordinate [fm]                  rho_new           rho_old           |error|       rel.error';
      write(6,'(a)') '';
      do il = 0 , NGL
          do ih = 0 , NGH
              ihl = 1+ih + il*(NGH+1);

              rho_new   = rho_s(ihl,1) + rho_s(ihl,2);
              rho_old   = ro(ihl,1);
              abs_error = abs( rho_new - rho_old );
              rel_error = abs_error / abs(rho_old);


              write(6,'(a,i2,a,1f6.2,a,a,i2,a,1f6.2,a,2f18.13,2e15.5)') '( r(' , il , ') = ' , rb(il) , ',  '  , &
                                                                        '  z(' , ih , ') = ' , zb(ih) , ' )  ' , &
                                                                        rho_new                                , &
                                                                        rho_old                                , &
                                                                        abs_error                              , &
                                                                        rel_error                              ;
          end do
      end do
#endif


      ! Deallocating density matrix.
      do it = 1 , 2
          do ib = 1 , N_blocks
              deallocate( rho1(it)%blocks(ib,ib)%mat );
          end do
          deallocate( rho1(it)%blocks );
      end do


      return;
      end subroutine check_gs_dens
