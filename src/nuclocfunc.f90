!======================================================================!

      subroutine calculate_nuclocfunc()

!======================================================================!
      use dataTypes;
      use fam_input;
      use simplex;
      use v_matrix;
      use basis;
      use drho;
      use quadrature;
      use nuclearlocfunc;
      implicit none;

      double precision         , parameter    :: pi = 4.d0*atan(1.d0);
      double complex           , parameter    :: cone  = cmplx( 1.d0 , 0.d0 , kind=8 );
      double complex           , parameter    :: czero = cmplx( 0.d0 , 0.d0 , kind=8 );
      type(complexBlockMatrix) , dimension(2) :: rho0;
      integer                                 :: it, ib, ib1, ib2, i1, i2;
      integer                                 :: ih, il;
      character                               :: fg1, fg2;
      integer                                 :: ml1, ml2;
      double precision                        :: r;
      double precision                        :: phiz1, phiz2, dphiz1, dphiz2;
      double precision                        :: phir1, phir2, dphir1, dphir2;
      double precision                        :: rho;
      double complex                          :: drho12;
      double precision                        :: f0z   ( -NGH:+NGH , 1:NGL , 1:2 );
      double precision                        :: f0r   ( -NGH:+NGH , 1:NGL , 1:2 );
      double precision                        :: g0    ( -NGH:+NGH , 1:NGL , 1:2 );
      double precision                        :: h0    ( -NGH:+NGH , 1:NGL , 1:2 );
      double precision                        :: D0    ( -NGH:+NGH , 1:NGL , 1:2 );
      double precision                        :: tau0TF( -NGH:+NGH , 1:NGL , 1:2 );
      double precision                        :: F0    ( -NGH:+NGH , 1:NGL , 1:2 );
      double complex                          :: dF    ( -NGH:+NGH , 1:NGL , 1:2 );
      double complex                          :: dfz   ( -NGH:+NGH , 1:NGL , 1:2 );
      double complex                          :: dfr   ( -NGH:+NGH , 1:NGL , 1:2 );
      double complex                          :: dg    ( -NGH:+NGH , 1:NGL , 1:2 );
      double complex                          :: dh    ( -NGH:+NGH , 1:NGL , 1:2 );
      double complex                          :: dD    ( -NGH:+NGH , 1:NGL , 1:2 );




      ! Calculating ground state density matrix rho0 = v*v'.
      do it = 1 , 2
          allocate( rho0(it)%blocks(1:N_blocks,1:N_blocks) );
          do ib = 1 , N_blocks
              allocate( rho0(it)%blocks(ib,ib)%mat(1:id_spx(ib),1:id_spx(ib)) );

              call zgemm( 'n' , 'c' , id_spx(ib) , id_spx(ib) , id_spx(ib)                    , &
                          cone                                                                , &
                          v(it)%blocks(ib,ib)%mat(1,1)    , size(v(it)%blocks(ib,ib)%mat,1)   , &
                          v(it)%blocks(ib,ib)%mat(1,1)    , size(v(it)%blocks(ib,ib)%mat,1)   , &
                          czero                                                               , &
                          rho0(it)%blocks(ib,ib)%mat(1,1) , size(rho0(it)%blocks(ib,ib)%mat,1)  );
          end do
      end do




      ! Calculating f0z, f0r, g0 and h0.
      f0z(:,:,:) = 0.d0;
      f0r(:,:,:) = 0.d0;
      g0 (:,:,:) = 0.d0;
      h0 (:,:,:) = 0.d0;
      do it = 1 , 2
          do ib = 1 , N_blocks
              do i2 = 1 , id_spx(ib)
                  do i1 = 1 , id_spx(ib)

                      fg1 = fg_spx(ib)%index(i1);
                      ml1 = ml_spx(ib)%index(i1);

                      fg2 = fg_spx(ib)%index(i2);
                      ml2 = ml_spx(ib)%index(i2);

                      if( fg1==fg2 .and. abs(ml1-ml2)==0 ) then

                          do il = 1 , NGL
                              do ih = -NGH , +NGH
                                  if( ih == 0 ) cycle;

                                  r = rb_fam(il);
                                  call assert( r > 1.d-8 , 'Possible division by zero in calculate_nuclocfunc.' );

                                   phiz1 =  phi_z(ib)%mat(ih,i1);
                                  dphiz1 = dphi_z(ib)%mat(ih,i1);
                                   phir1 =  phi_r(ib)%mat(il,i1);
                                  dphir1 = dphi_r(ib)%mat(il,i1);

                                   phiz2 =  phi_z(ib)%mat(ih,i2);
                                  dphiz2 = dphi_z(ib)%mat(ih,i2);
                                   phir2 =  phi_r(ib)%mat(il,i2);
                                  dphir2 = dphi_r(ib)%mat(il,i2);

                                  rho = real(rho0(it)%blocks(ib,ib)%mat(i1,i2));

                                  f0z(ih,il,it) = f0z(ih,il,it) + rho * (               (dphiz1*phiz2+phiz1*dphiz2)/2 *   phir1* phir2                 ) * 1/(2*pi);
                                  f0r(ih,il,it) = f0r(ih,il,it) + rho * (                 phiz1* phiz2                * (dphir1*phir2+phir1*dphir2)/2  ) * 1/(2*pi);
                                  g0 (ih,il,it) = g0 (ih,il,it) + rho * (                 phiz1* phiz2                *   phir1* phir2                 ) * 1/(2*pi);
                                  h0 (ih,il,it) = h0 (ih,il,it) + rho * (                dphiz1*dphiz2                *   phir1* phir2               + &
                                                                                          phiz1* phiz2                *  dphir1*dphir2               + &
                                                                          ml1*ml2/r**2 *  phiz1* phiz2                *   phir1* phir2                 ) * 1/(2*pi);

                              end do
                          end do

                      end if

                  end do
              end do
          end do
      end do




      ! Calculating dfz, dfr, dg and dh.
      dfz(:,:,:) = czero;
      dfr(:,:,:) = czero;
      dg (:,:,:) = czero;
      dh (:,:,:) = czero;
      do it = 1 , 2
          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  call assert( drho1(it)%nnzblocks(ib1,ib2) .eqv. drho2(it)%nnzblocks(ib1,ib2) , 'drho1 and drho2 nnz blocks error.' );
                  if( drho1(it)%nnzblocks(ib1,ib2) .and. drho2(it)%nnzblocks(ib1,ib2) ) then

                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , id_spx(ib1)

                              fg1 = fg_spx(ib1)%index(i1);
                              ml1 = ml_spx(ib1)%index(i1);

                              fg2 = fg_spx(ib2)%index(i2);
                              ml2 = ml_spx(ib2)%index(i2);

                              if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then

                                  do il = 1 , NGL
                                      do ih = -NGH , +NGH
                                          if( ih == 0 ) cycle;

                                          r = rb_fam(il);
                                          call assert( r > 1.d-8 , 'Possible division by zero in calculate_nuclocfunc.' );

                                           phiz1 =  phi_z(ib1)%mat(ih,i1);
                                          dphiz1 = dphi_z(ib1)%mat(ih,i1);
                                           phir1 =  phi_r(ib1)%mat(il,i1);
                                          dphir1 = dphi_r(ib1)%mat(il,i1);

                                           phiz2 =  phi_z(ib2)%mat(ih,i2);
                                          dphiz2 = dphi_z(ib2)%mat(ih,i2);
                                           phir2 =  phi_r(ib2)%mat(il,i2);
                                          dphir2 = dphi_r(ib2)%mat(il,i2);

                                          drho12 = ( drho1(it)%blocks(ib1,ib2)%mat(i1,i2) + drho2(it)%blocks(ib1,ib2)%mat(i1,i2) )/2;

                                          dfz(ih,il,it) = dfz(ih,il,it) + drho12 * (               (dphiz1*phiz2+phiz1*dphiz2)/2 *   phir1* phir2                ) * 1/(2*pi);
                                          dfr(ih,il,it) = dfr(ih,il,it) + drho12 * (                 phiz1* phiz2                * (dphir1*phir2+phir1*dphir2)/2 ) * 1/(2*pi);
                                          dg (ih,il,it) = dg (ih,il,it) + drho12 * (                 phiz1* phiz2                *   phir1* phir2                ) * 1/(2*pi);
                                          dh (ih,il,it) = dh (ih,il,it) + drho12 * (                dphiz1*dphiz2                *   phir1* phir2              + &
                                                                                                     phiz1* phiz2                *  dphir1*dphir2              + &
                                                                                     ml1*ml2/r**2 *  phiz1* phiz2                *   phir1* phir2                ) * 1/(2*pi);

                                      end do
                                  end do

                              end if

                          end do
                      end do

                  end if
              end do
          end do
      end do




      ! Calculating D0 and dD.
      D0 = h0 -  ( f0r**2 + f0z**2 )/g0;
      dD = dh + (( f0r**2 + f0z**2 )/g0**2)*dg - 2*( f0r*dfr + f0z*dfz )/g0;




      ! Calculating F0 and dF.
      tau0TF = (3.d0/5.d0) * (6*pi**2)**(2.d0/3.d0) * g0**(5.d0/3.d0);
      F0     = D0/tau0TF;
      dF     = dD/tau0TF - (5.d0/3.d0)*(D0/tau0TF)*(dg/g0);



      ! Calculating C0 and dC.
      C0 = (   1  / (1+F0**2)    );
      dC = (-2*F0 / (1+F0**2)**2 )*dF;




      ! Deallocating rho0.
      do it = 1 , 2
          do ib = 1 , N_blocks
              deallocate( rho0(it)%blocks(ib,ib)%mat );
          end do
          deallocate( rho0(it)%blocks );
      end do




      return;
      end subroutine calculate_nuclocfunc
