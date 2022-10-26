!======================================================================!

      subroutine init_pairing()

!======================================================================!
      use fam_input;
      use simplex;
      use basis;
      use pairparams;
      use Wpairing;
      use dDelta;
      implicit none;

      integer                                               :: nz1    , nz2    , N_z    , nz;
      integer                                               :: nz1_max, nz2_max, N_z_max;
      integer                                               :: nr1    , nr2    , ml1    , ml2    , N_r    , nr;
      integer                                               :: nr1_max, nr2_max, ml1_max, ml2_max, N_r_max;
      integer                                               :: iW;
      character                                             :: fg1, fg2;
      integer                                               :: it, ib1, ib2, i1, i2;
      double precision ,                        parameter   :: pi = 4.d0*atan(1.d0);
      double precision , dimension(:,:,:)     , allocatable :: Vz;
      double precision , dimension(:,:,:,:,:) , allocatable :: Vr;
      double precision ,                        external    :: factorial;
      double precision ,                        external    :: TalmiMoshinsky_z;
      double precision ,                        external    :: TalmiMoshinsky_r;




      if( include_pairing == 0 ) then
          return;
      end if




      nz1_max = n0f;
      nz2_max = n0f;
      N_z_max = nz1_max + nz2_max;
      allocate( Vz( 0:N_z_max , 0:nz2_max , 0:nz1_max ) );

      Vz(:,:,:) = 0.d0;
      do nz1 = 0 , nz1_max
          do nz2 = 0 , nz2_max
              do N_z = 0 , nz1+nz2
                  if( mod(N_z,2) == mod(nz1+nz2,2) ) then

                      nz = nz1+nz2 - N_z;

                      Vz( N_z , nz2 , nz1 ) =   sqrt(factorial(nz)) / ( 2**(nz/2) * factorial(nz/2) ) &
                                              *       ( a_pairing**2 - bz**2 )**(nz/2)                &
                                              / sqrt( ( a_pairing**2 + bz**2 )**(nz+1) )              &
                                              * TalmiMoshinsky_z( nz1 , nz2 , N_z , nz );
                  end if
              end do
          end do
      end do




      nr1_max = n0f/2;
      nr2_max = n0f/2;
      ml1_max = n0f;
      ml2_max = n0f;
      N_r_max = n0f; ! 2*Nr+|Ml| + 2*nr+|ml| == 2*nr1+|ml1| + 2*nr2+|ml2| <= n0f + n0f = 2*n0f, thus Nr <= n0f.
      allocate( Vr( 0:N_r_max , -ml2_max:+ml2_max , 0:nr2_max , -ml1_max:+ml1_max , 0:nr1_max ) );

      Vr(:,:,:,:,:) = 0.d0;
      do nr1 = 0 , nr1_max
          do ml1 = -ml1_max , +ml1_max
              if( 2*nr1 + abs(ml1) <= n0f ) then
                  do nr2 = 0 , nr2_max
                      do ml2 = -ml2_max , +ml2_max
                          if( 2*nr2 + abs(ml2) <= n0f ) then
                              if( abs(ml1-ml2) == K_multipole ) then ! See Eqs. (E.4) and (E.5).
                                  do N_r = 0 , ( 2*nr1+abs(ml1) + 2*nr2+abs(ml2) - K_multipole )/2

                                      nr = ( 2*nr1+abs(ml1) + 2*nr2+abs(ml2) - (2*N_r+abs(ml1-ml2)) )/2;

                                      Vr( N_r , ml2 , nr2 , ml1 , nr1 ) =   ( bp**2 - a_pairing**2 )**nr       &
                                                                          / ( bp**2 + a_pairing**2 )**(nr+1)   &
                                                                          * TalmiMoshinsky_r( nr1 , +ml1     , &
                                                                                              nr2 , -ml2     , &
                                                                                              N_r , +ml1-ml2 , &
                                                                                              nr  ,  0         );
                                  end do
                              end if
                          end if
                      end do
                  end do
              end if
          end do
      end do




      ! Calculating W coefficients, see Eq. (E.2).
      ! Keep in mind that this exact same loop strategy must be
      ! consistently used in fam_ddelta() subroutine, otherwise
      ! the order of index iW can be different.
      iW = 0;
      do ib2 = 1 , N_blocks
          do ib1 = 1 , ib2
              if( all([(dDelta1_pl(it)%nnzblocks(ib1,ib2),it=1,2)]) .and. all([(dDelta1_mi(it)%nnzblocks(ib1,ib2),it=1,2)]) ) then

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

                                          iW    = iW + 1;
                                          W(iW) =   bp*sqrt(bz)/((2*pi)**(0.75d0))  &
                                                  * Vz( N_z , nz2 , nz1 )           &
                                                  * Vr( N_r , ml2 , nr2 , ml1 , nr1 );
                                      end if
                                  end do
                              end do
                          end if

                      end do
                  end do

              end if
          end do
      end do
      call assert( iW == size(W) , 'iW /= size(W) in init_pairing.' );

      deallocate( Vz );
      deallocate( Vr );




      return;
      end subroutine init_pairing
