!======================================================================!

      subroutine construct_u()

!======================================================================!
      use simplex;
      use u_matrix;
      use u_energy;
      use blokap;
      use bloosc;
      use quaosc;
      use blodir;
      use waveuv;
      implicit none;
      integer , external :: index_of_vector;
      integer            :: it, ib;
      integer            :: i, j;
      integer            :: nf, ng, i0f, kap, nh, klp, kla, k, n;
      character          :: fg;
      integer            :: nzz, nrr, mll, mx;
      double complex     :: z;


      do it = 1 , 2
          do ib = 1 , N_blocks

              u(it)%blocks(ib,ib)%mat(:,:) = cmplx( 0.d0 , 0.d0 , kind=8 );

              j = 1; ! Index of column of u(it)%blocks(ib,ib)%mat(:,j).

              nf  = id(ib,1);
              ng  = id(ib,2);
              i0f = ia(ib,1);
              kap = kb(ib);
              nh  = nf + ng;

              ! Particles (positive energies in Fermi sea).
              klp = ka(ib,it);
              do k = 1 , nf
                  klp = klp + 1;

                  E_fam_u(it)%blocks(ib)%vec(j) = equ(klp,it);

                  do n = 1 , nh
                      fg = merge( 'f' , 'g' , n<=nf ); ! (n<=nf) ? 'f' : 'g';

                      nzz = nz( i0f + n );
                      nrr = nr( i0f + n );
                      mll = ml( i0f + n );
                      mx  = 2*( abs(kap) - mll ) - 1; ! abs(kap) = Omega+1/2 = ml+ms+1/2, thus mx=2*ms.

                      ! Eq. (A.17).
                      z = cmplx( +fguv(n,klp,it) , 0.d0 , kind=8 );
                      if( fg=='f' .and. mx==+1 ) z = z * cmplx( +1.d0 , +0.d0 , kind=8 );
                      if( fg=='f' .and. mx==-1 ) z = z * cmplx( +0.d0 , +1.d0 , kind=8 );
                      if( fg=='g' .and. mx==+1 ) z = z * cmplx( +0.d0 , +1.d0 , kind=8 );
                      if( fg=='g' .and. mx==-1 ) z = z * cmplx( +1.d0 , +0.d0 , kind=8 );


                      i = index_of_vector( fg , nzz , nrr , mx*mll , ib );

                      u(it)%blocks(ib,ib)%mat(i,j) = z;

                  end do

                  j = j + 1;
              end do

              ! Antiparticles (negative energies in Dirac sea).
              kla = ka(ib,it+2);
              do k = 1 , ng
                  kla = kla + 1;

                  E_fam_u(it)%blocks(ib)%vec(j) = equ(kla,it+2);

                  do n = 1 , nh
                      fg = merge( 'f' , 'g' , n<=nf ); ! (n<=nf) ? 'f' : 'g';

                      nzz = nz( i0f + n );
                      nrr = nr( i0f + n );
                      mll = ml( i0f + n );
                      mx  = 2*( abs(kap) - mll ) - 1; ! abs(kap) = Omega+1/2 = ml+ms+1/2, thus mx=2*ms.

                      ! Eq. (A.17).
                      z = cmplx( +fguv(n,kla,it+2) , 0.d0 , kind=8 );
                      if( fg=='f' .and. mx==+1 ) z = z * cmplx( +1.d0 , +0.d0 , kind=8 );
                      if( fg=='f' .and. mx==-1 ) z = z * cmplx( +0.d0 , +1.d0 , kind=8 );
                      if( fg=='g' .and. mx==+1 ) z = z * cmplx( +0.d0 , +1.d0 , kind=8 );
                      if( fg=='g' .and. mx==-1 ) z = z * cmplx( +1.d0 , +0.d0 , kind=8 );


                      i = index_of_vector( fg , nzz , nrr , mx*mll , ib );

                      u(it)%blocks(ib,ib)%mat(i,j) = z;

                  end do

                  j = j + 1;
              end do

              call assert( j-1 == id_spx(ib) , 'Error in construct_u.' );

          end do
      end do


      return;
      end subroutine construct_u






!======================================================================!

      subroutine construct_v()

!======================================================================!
      use simplex;
      use v_matrix;
      use v_energy;
      use blokap;
      use bloosc;
      use quaosc;
      use blodir;
      use waveuv;
      implicit none;
      integer , external :: index_of_vector;
      integer            :: it, ib;
      integer            :: i, j;
      integer            :: nf, ng, i0f, kap, nh, klp, kla, k, n;
      character          :: fg;
      integer            :: nzz, nrr, mll, mx;
      double complex     :: z;


      do it = 1 , 2
          do ib = 1 , N_blocks

              v(it)%blocks(ib,ib)%mat(:,:) = cmplx( 0.d0 , 0.d0 , kind=8 );

              j = 1; ! Index of column of v(it)%blocks(ib,ib)%mat(:,j).

              nf  = id(ib,1);
              ng  = id(ib,2);
              i0f = ia(ib,1);
              kap = kb(ib);
              nh  = nf + ng;

              ! Particles (positive energies in Fermi sea).
              klp = ka(ib,it);
              do k = 1 , nf
                  klp = klp + 1;

                  E_fam_v(it)%blocks(ib)%vec(j) = equ(klp,it);

                  do n = nh+1 , nh+nh
                      fg = merge( 'f' , 'g' , n<=nh+nf ); ! (n<=nf) ? 'f' : 'g';

                      nzz = nz( i0f + n-nh );
                      nrr = nr( i0f + n-nh );
                      mll = ml( i0f + n-nh );
                      mx  = 2*( abs(kap) - mll ) - 1; ! abs(kap) = Omega+1/2 = ml+ms+1/2, thus mx=2*ms.

                      ! Eq. (A.18).
                      z = cmplx( +fguv(n,klp,it) , 0.d0 , kind=8 );
                      if( fg=='f' .and. mx==+1 ) z = z * cmplx( -1.d0 , +0.d0 , kind=8 );
                      if( fg=='f' .and. mx==-1 ) z = z * cmplx( +0.d0 , -1.d0 , kind=8 );
                      if( fg=='g' .and. mx==+1 ) z = z * cmplx( +0.d0 , -1.d0 , kind=8 );
                      if( fg=='g' .and. mx==-1 ) z = z * cmplx( -1.d0 , +0.d0 , kind=8 );


                      i = index_of_vector( fg , nzz , nrr , mx*mll , ib );

                      v(it)%blocks(ib,ib)%mat(i,j) = z;

                  end do

                  j = j + 1;
              end do

              ! Antiparticles (negative energies in Dirac sea).
              kla = ka(ib,it+2);
              do k = 1 , ng
                  kla = kla + 1;

                  E_fam_v(it)%blocks(ib)%vec(j) = equ(kla,it+2);

                  do n = nh+1 , nh+nh
                      fg = merge( 'f' , 'g' , n<=nh+nf ); ! (n<=nf) ? 'f' : 'g';

                      nzz = nz( i0f + n-nh );
                      nrr = nr( i0f + n-nh );
                      mll = ml( i0f + n-nh );
                      mx  = 2*( abs(kap) - mll ) - 1; ! abs(kap) = Omega+1/2 = ml+ms+1/2, thus mx=2*ms.

                      ! Eq. (A.18).
                      z = cmplx( +fguv(n,kla,it+2) , 0.d0 , kind=8 );
                      if( fg=='f' .and. mx==+1 ) z = z * cmplx( -1.d0 , +0.d0 , kind=8 );
                      if( fg=='f' .and. mx==-1 ) z = z * cmplx( +0.d0 , -1.d0 , kind=8 );
                      if( fg=='g' .and. mx==+1 ) z = z * cmplx( +0.d0 , -1.d0 , kind=8 );
                      if( fg=='g' .and. mx==-1 ) z = z * cmplx( -1.d0 , +0.d0 , kind=8 );


                      i = index_of_vector( fg , nzz , nrr , mx*mll , ib );

                      v(it)%blocks(ib,ib)%mat(i,j) = z;

                  end do

                  j = j + 1;
              end do

              call assert( j-1 == id_spx(ib) , 'Error in construct_v.' );

          end do
      end do


      return;
      end subroutine construct_v






!======================================================================!

      subroutine construct_qpenergies()

!======================================================================!
      use simplex;
      use fam_energies;
      use u_energy;
      use v_energy;
      implicit none;
      integer          :: it, ib, i;
      double precision :: E1, E2;


      do it = 1 , 2
          do ib = 1 , N_blocks
              do i = 1 , id_spx(ib)

                  E1 = E_fam_u(it)%blocks(ib)%vec(i);
                  E2 = E_fam_v(it)%blocks(ib)%vec(i);
                  call assert( abs(E1-E2)<=1.d-15 , 'E1 =/= E2' );

                  E_fam(it)%blocks(ib)%vec(i) = E1;
              end do
          end do
      end do


      return;
      end subroutine construct_qpenergies
