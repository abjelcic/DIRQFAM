!======================================================================!

      subroutine base_simplex()

!======================================================================!
      use fam_input;
      use dirhbpar;
      use blokap;
      use bloosc;
      use quaosc;
      use simplex;
      implicit none;

      integer            :: ib, kap, nf, ng, i0f, i0g, il, i, j, mx;
      integer            :: nz1, nz2;
      integer            :: nz_tmp, nr_tmp, ml_tmp;
      integer            :: ib1, ib2, i1, i2;
      character          :: fg;
      integer            :: nzz, nrr, mll;
      integer , external :: index_of_vector;




      ! Total number of basis vectors in simplex s=+-i block.
      N_total = 0;
      do ib = 1 , nb
          nf = id(ib,1);
          ng = id(ib,2);
          N_total = N_total + nf+ng;
      end do

      ! Total number of subblocks.
      N_blocks = 2*n0f + 3;

      call assert( N_blocks == NBX , 'N_blocks =/= NBX' );
      call assert( N_blocks == nb  , 'N_blocks =/= nb'  );
      call assert( N_total  == NTX , 'N_total  =/= NTX' );




      ! Construction of simplex-y quantum numbers.
      allocate(  ia_spx(1:N_blocks) );
      allocate(  id_spx(1:N_blocks) );
      allocate( nf_size(1:N_blocks) );
      allocate( ng_size(1:N_blocks) );
      allocate(  nz_spx(1:N_blocks) );
      allocate(  nr_spx(1:N_blocks) );
      allocate(  ml_spx(1:N_blocks) );
      allocate(  fg_spx(1:N_blocks) );

      do ib = 1 , N_blocks

          kap = kb(ib);
          nf  = id(ib,1);
          ng  = id(ib,2);
          i0f = ia(ib,1);
          i0g = ia(ib,2);


          nf_size(ib) = nf;
          ng_size(ib) = ng;
          id_spx (ib) = nf + ng;

          allocate( nz_spx(ib)%index( 1:id_spx(ib) ) );
          allocate( nr_spx(ib)%index( 1:id_spx(ib) ) );
          allocate( ml_spx(ib)%index( 1:id_spx(ib) ) );
          allocate( fg_spx(ib)%index( 1:id_spx(ib) ) );

          il = 0;
          do i = i0f + 1 , i0f + nf
              il = il + 1;

              mx = 2*( abs(kap) - ml(i) ) - 1; ! mx=+1 if ms=+1/2, and mx=-1 if ms=-1/2.
              nz_spx(ib)%index(il) = nz(i);
              nr_spx(ib)%index(il) = nr(i);
              ml_spx(ib)%index(il) = ml(i) * mx;
              fg_spx(ib)%index(il) = 'f';
          end do
          do i = i0g + 1 , i0g + ng
              il = il + 1;

              mx = 2*( abs(kap) - ml(i) ) - 1; ! mx=+1 if ms=+1/2, and mx=-1 if ms=-1/2.
              nz_spx(ib)%index(il) = nz(i);
              nr_spx(ib)%index(il) = nr(i);
              ml_spx(ib)%index(il) = ml(i) * mx;
              fg_spx(ib)%index(il) = 'g';
          end do

      end do

      ia_spx(1) = 1;
      do ib = 2 , N_blocks
          ia_spx(ib) = ia_spx(ib-1) + id_spx(ib-1);
      end do

      allocate(            getBlock(1:N_total) );
      allocate( getIndexWithinBlock(1:N_total) );
      il = 0;
      do ib = 1 , N_blocks
          do i = 1 , id_spx(ib)
              il = il + 1;

              getBlock(il)            = ib;
              getIndexWithinBlock(il) = i;
          end do
      end do
      call assert( il==N_total , 'il =/= N_total' );




      ! Sort basis blocks with respect to nz parity.
      ! This sort makes U,V matrices to appear with purely real or imaginary blocks.
      ! Current code does not use this info and it is not necessary.
      do ib = 1 , N_blocks

          do i = 1 , nf_size(ib)
              do j = i+1 , nf_size(ib)

                  nz1 = nz_spx(ib)%index(i);
                  nz2 = nz_spx(ib)%index(j);

                  if( mod(nz1,2)==1 .and. mod(nz2,2)==0 ) then
                      nz_tmp              = nz_spx(ib)%index(i);
                      nz_spx(ib)%index(i) = nz_spx(ib)%index(j);
                      nz_spx(ib)%index(j) = nz_tmp;

                      nr_tmp              = nr_spx(ib)%index(i);
                      nr_spx(ib)%index(i) = nr_spx(ib)%index(j);
                      nr_spx(ib)%index(j) = nr_tmp;

                      ml_tmp              = ml_spx(ib)%index(i);
                      ml_spx(ib)%index(i) = ml_spx(ib)%index(j);
                      ml_spx(ib)%index(j) = ml_tmp;
                  end if

              end do
          end do

          do i = nf_size(ib)+1 , nf_size(ib)+ng_size(ib)
              do j = i+1 , nf_size(ib)+ng_size(ib)

                  nz1 = nz_spx(ib)%index(i);
                  nz2 = nz_spx(ib)%index(j);

                  if( mod(nz1,2)==1 .and. mod(nz2,2)==0 ) then
                      nz_tmp              = nz_spx(ib)%index(i);
                      nz_spx(ib)%index(i) = nz_spx(ib)%index(j);
                      nz_spx(ib)%index(j) = nz_tmp;

                      nr_tmp              = nr_spx(ib)%index(i);
                      nr_spx(ib)%index(i) = nr_spx(ib)%index(j);
                      nr_spx(ib)%index(j) = nr_tmp;

                      ml_tmp              = ml_spx(ib)%index(i);
                      ml_spx(ib)%index(i) = ml_spx(ib)%index(j);
                      ml_spx(ib)%index(j) = ml_tmp;
                  end if

              end do
          end do

      end do




      ! Test whether the basis has been correctly constructed.
      if( .true. ) then

#ifdef DEBUG
          do ib = 1 , N_blocks
              write(6,'(a)') '------------------------------------------';
              do i = 1 , id_spx(ib)
                  write(6,'(i4,a,a,a,i3,a,i3,a,i3)') ia_spx(ib)+i-1 , '. fg = ', fg_spx(ib)%index(i), &
                                                                      '  nz = ', nz_spx(ib)%index(i), &
                                                                      '  nr = ', nr_spx(ib)%index(i), &
                                                                      '  ml = ', ml_spx(ib)%index(i);
              end do
          end do
#endif

          do ib1 = 1 , N_blocks
              do ib2 = 1 , N_blocks
                  do i1 = 1 , id_spx(ib1)
                      do i2 = 1 , id_spx(ib2)

                          if( ib1==ib2 .and. i1==i2 ) cycle;

                          if( fg_spx(ib1)%index(i1) == fg_spx(ib2)%index(i2) .and. &
                              nz_spx(ib1)%index(i1) == nz_spx(ib2)%index(i2) .and. &
                              nr_spx(ib1)%index(i1) == nr_spx(ib2)%index(i2) .and. &
                              ml_spx(ib1)%index(i1) == ml_spx(ib2)%index(i2)       ) then
                              stop 'Error: Redundancy in basis!';
                          end if

                      end do
                  end do
              end do
          end do

          do ib = 1 , N_blocks
              do i = 1 , id_spx(ib)

                  fg  = fg_spx(ib)%index(i);
                  nzz = nz_spx(ib)%index(i);
                  nrr = nr_spx(ib)%index(i);
                  mll = ml_spx(ib)%index(i);

                  if( fg == 'f' ) then
                      if( .not.( nzz + 2*nrr + abs(mll) <= n0f   ) ) then
                          stop 'Error: Basis vector out of shell!';
                      end if
                  end if
                  if( fg == 'g' ) then
                      if( .not.( nzz + 2*nrr + abs(mll) <= n0f+1 ) ) then
                          stop 'Error: Basis vector out of shell!';
                      end if
                  end if

              end do
          end do

          do ib = 1 , N_blocks
              do i = 1 , id_spx(ib)
                  fg  = fg_spx(ib)%index(i);
                  nzz = nz_spx(ib)%index(i);
                  nrr = nr_spx(ib)%index(i);
                  mll = ml_spx(ib)%index(i);

                  call assert( i==index_of_vector(fg,nzz,nrr,mll,ib) , 'index_of_vector() problem!' );

              end do
          end do

      end if




      return;
      end subroutine base_simplex
