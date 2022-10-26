!======================================================================!

      subroutine fam_dh()

!======================================================================!
      use fam_input;
      use simplex;
      use wbasis;
      use dpotentials;
      use dh;
      implicit none;

      double complex   , parameter :: II = cmplx( 0.d0 , 1.d0 , kind=8 );
      integer                      :: z_parity;
      integer                      :: it, ib1, ib2, i1, i2;
      integer                      :: ih, il;
      character                    :: fg1, fg2;
      integer                      :: nz1, nz2;
      integer                      :: nr1, nr2;
      integer                      :: ml1, ml2;
      double precision             :: fac;
      double complex               :: integral;
      !                                    ( ih    , il    , z-parity , it  )
      double complex               :: dVpS ( 1:NGH , 1:NGL , 0:1      , 1:2 );
      double complex               :: dVmS ( 1:NGH , 1:NGL , 0:1      , 1:2 );
      double complex               :: dSigr( 1:NGH , 1:NGL , 0:1      , 1:2 );
      double complex               :: dSigp( 1:NGH , 1:NGL , 0:1      , 1:2 );
      double complex               :: dSigz( 1:NGH , 1:NGL , 0:1      , 1:2 );
      double complex               :: dSig1( 1:NGH , 1:NGL , 0:1      , 1:2 );
      double complex               :: dSig2( 1:NGH , 1:NGL , 0:1      , 1:2 );
      double complex               :: dSig3( 1:NGH , 1:NGL , 0:1      , 1:2 );




      do it = 1 , 2
          do z_parity = 0 , 1
              do il = 1 , NGL
                  do ih = 1 , NGH

                       dVpS(ih,il,z_parity,it) =  dVplusS(+ih,il,it) + merge(+1.d0,-1.d0,mod(z_parity,2)==0) *  dVplusS(-ih,il,it);
                       dVmS(ih,il,z_parity,it) = dVminusS(+ih,il,it) + merge(+1.d0,-1.d0,mod(z_parity,2)==0) * dVminusS(-ih,il,it);
                      dSigr(ih,il,z_parity,it) = dSigma_r(+ih,il,it) + merge(+1.d0,-1.d0,mod(z_parity,2)==0) * dSigma_r(-ih,il,it);
                      dSigp(ih,il,z_parity,it) = dSigma_p(+ih,il,it) + merge(+1.d0,-1.d0,mod(z_parity,2)==0) * dSigma_p(-ih,il,it);
                      dSigz(ih,il,z_parity,it) = dSigma_z(+ih,il,it) + merge(+1.d0,-1.d0,mod(z_parity,2)==0) * dSigma_z(-ih,il,it);

                  end do
              end do
          end do
      end do
      dSig1(:,:,:,:) = ( + dSigr(:,:,:,:) - dSigp(:,:,:,:) ) / sqrt(2.d0);
      dSig2(:,:,:,:) = ( + dSigr(:,:,:,:) + dSigp(:,:,:,:) ) / sqrt(2.d0);
      dSig3(:,:,:,:) = ( + dSigz(:,:,:,:)                  );

#ifdef DEBUG
      if( K_multipole == 0 ) call assert( maxval(abs(dSig1(:,:,:,:)-dSig2(:,:,:,:))) < 1.d-12 , 'For K=0 dSigma_1 is not equal to dSigma_2!' );
#endif




      ! Calculation of dh1 and dh2.
      do it = 1 , 2

          ! First calculate the upper triangle of dh1.
          do ib2 = 1 , N_blocks
              do ib1 = 1 , ib2
                  if( dh1(it)%nnzblocks(ib1,ib2) ) then

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

                              fac = merge( 1.d0 , 0.5d0 , K_multipole==0 );

                              integral = cmplx( 0.d0 , 0.d0 , kind=8 );

                              if( all([fg1,fg2]==['f','f']) .and. abs(ml1-ml2)==+K_multipole ) integral = fac*(+1.d0         )*sum(  dVpS(:,:,mod(nz1+nz2,2),it) * wPhi(ib1)%arr(:,:,i1)*wPhi(ib2)%arr(:,:,i2) );
                              if( all([fg1,fg2]==['g','g']) .and. abs(ml1-ml2)==+K_multipole ) integral = fac*(+1.d0         )*sum(  dVmS(:,:,mod(nz1+nz2,2),it) * wPhi(ib1)%arr(:,:,i1)*wPhi(ib2)%arr(:,:,i2) );
                              if( all([fg1,fg2]==['f','g']) .and.    ml1+ml2+1==+K_multipole ) integral = fac*(-sqrt(2.d0)*II)*sum( dSig1(:,:,mod(nz1+nz2,2),it) * wPhi(ib1)%arr(:,:,i1)*wPhi(ib2)%arr(:,:,i2) );
                              if( all([fg1,fg2]==['g','f']) .and.    ml1+ml2+1==+K_multipole ) integral = fac*(+sqrt(2.d0)*II)*sum( dSig1(:,:,mod(nz1+nz2,2),it) * wPhi(ib1)%arr(:,:,i1)*wPhi(ib2)%arr(:,:,i2) );
                              if( all([fg1,fg2]==['f','g']) .and.    ml1+ml2+1==-K_multipole ) integral = fac*(-sqrt(2.d0)*II)*sum( dSig2(:,:,mod(nz1+nz2,2),it) * wPhi(ib1)%arr(:,:,i1)*wPhi(ib2)%arr(:,:,i2) );
                              if( all([fg1,fg2]==['g','f']) .and.    ml1+ml2+1==-K_multipole ) integral = fac*(+sqrt(2.d0)*II)*sum( dSig2(:,:,mod(nz1+nz2,2),it) * wPhi(ib1)%arr(:,:,i1)*wPhi(ib2)%arr(:,:,i2) );
                              if( all([fg1,fg2]==['f','g']) .and. abs(ml1-ml2)==+K_multipole ) integral = fac*(-1.d0         )*sum( dSig3(:,:,mod(nz1+nz2,2),it) * wPhi(ib1)%arr(:,:,i1)*wPhi(ib2)%arr(:,:,i2) );
                              if( all([fg1,fg2]==['g','f']) .and. abs(ml1-ml2)==+K_multipole ) integral = fac*(-1.d0         )*sum( dSig3(:,:,mod(nz1+nz2,2),it) * wPhi(ib1)%arr(:,:,i1)*wPhi(ib2)%arr(:,:,i2) );

                              dh1(it)%blocks(ib1,ib2)%mat(i1,i2) = integral;

                          end do
                      end do

                  end if
              end do
          end do

          ! Construct the full dh1 using its upper triangle.
          call constructFulldh1FromItsUpperTriangle( dh1(it) );

          ! Construct dh2 using dh1.
          call constructFulldh2FromFulldh1( dh1(it) , dh2(it) );

      end do




      return;
      end subroutine fam_dh






!======================================================================!

      subroutine constructFulldh1FromItsUpperTriangle( dh1 )

!======================================================================!
      use dataTypes;
      use fam_input;
      use simplex;
      implicit none;
      type(complexBlockMatrix) , intent(inout) :: dh1;
      integer                                  :: ib1, ib2, i1, i2;
      character                                :: fg1, fg2;
      integer                                  :: ml1, ml2;

      do ib2 = 1 , N_blocks
          do ib1 = ib2 , N_blocks
              if( dh1%nnzblocks(ib1,ib2) ) then
                  do i2 = 1 , id_spx(ib2)
                      do i1 = merge( i2 , 1 , ib1==ib2 ) , id_spx(ib1) ! Loop over the lower triangle.

                          fg1 = fg_spx(ib1)%index(i1);
                          ml1 = ml_spx(ib1)%index(i1);

                          fg2 = fg_spx(ib2)%index(i2);
                          ml2 = ml_spx(ib2)%index(i2);

                          dh1%blocks(ib1,ib2)%mat(i1,i2) = merge( -1.d0 , +1.d0 , fg1/=fg2 .and. abs(ml1+ml2+1)==K_multipole ) * dh1%blocks(ib2,ib1)%mat(i2,i1);

                      end do
                  end do
              end if
          end do
      end do

      return;
      end subroutine constructFulldh1FromItsUpperTriangle






!======================================================================!

      subroutine constructFulldh2FromFulldh1( dh1 , dh2 )

!======================================================================!
      use dataTypes;
      use fam_input;
      use simplex;
      implicit none;
      type(complexBlockMatrix) , intent(in)    :: dh1;
      type(complexBlockMatrix) , intent(inout) :: dh2;
      integer                                  :: ib1, ib2, i1, i2;
      character                                :: fg1, fg2;
      integer                                  :: ml1, ml2;

      do ib2 = 1 , N_blocks
          do ib1 = 1 , N_blocks
              if( dh1%nnzblocks(ib1,ib2) ) then
                  do i2 = 1 , id_spx(ib2)
                      do i1 = 1 , id_spx(ib1)

                          fg1 = fg_spx(ib1)%index(i1);
                          ml1 = ml_spx(ib1)%index(i1);

                          fg2 = fg_spx(ib2)%index(i2);
                          ml2 = ml_spx(ib2)%index(i2);

                          dh2%blocks(ib1,ib2)%mat(i1,i2) = merge( -1.d0 , +1.d0 , fg1/=fg2 .and. abs(ml1-ml2)==K_multipole ) * dh1%blocks(ib1,ib2)%mat(i1,i2);

                      end do
                  end do
              end if
          end do
      end do

      return;
      end subroutine constructFulldh2FromFulldh1
