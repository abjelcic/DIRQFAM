!======================================================================!

      subroutine kpm( Omega_b , N_it )

!======================================================================!
      use simplex;
      use f02f20matrix;
      use KPMdata;
      implicit none;

      double precision , intent(in) :: Omega_b;
      integer          , intent(in) :: N_it;
      integer                       :: tape_mu;
      double precision              :: mu(0:2*N_it);
      integer                       :: it;
      integer                       :: n;


      ! We have X = [ 0 , x ; -x^T , 0 ] and Y = [ 0 , y ; -y^T , 0 ], see Eq. (B.7).
      ! The same holds for F20 = [ 0 , f20 ; -f20^T , 0 ] and F02 = [ 0 , f02 ; -f02^T , 0 ] (see Eq. (B.6)) and
      ! dH20 = [ 0 , dh20 , -dh20^T , 0 ] and dH02 = [ 0 , dh02 , -dh02^T , 0 ] (see. Eq. (B.13)).
      !
      ! If from (X,Y) we calculate (dH20,dH02) then there holds:
      ! ( [ A , B ; B^* , A^* ] * [X;Y] )_ij = [ dH20_ij + (E_i+E_j)X_ij ; dH02_ij +0 (Ei+Ej)Y_ij ].
      ! See e.g. Eq. (9) in N. Hinohara, M. Kortelainen, W. Nazarewicz, Phys. Rev. C 87, 064309 (2013).
      !
      ! This is how we extract the mapping [X,Y] |-> [ A , B ; B^* , A^* ] * [X;Y] needed for KPM,
      ! we simply calculate (dH20,dH02) for a given (X,Y), and add terms (E_i+E_j)X_ij and (E_i+E_j)Y_ij.
      ! Notice that here vectorized matrices use only lower triangular part when vectorizing because
      ! in sums over (mu,nu) indices, we have mu<nu. That is why instead of (X,Y,F20,F02,dH20,dH02),
      ! we can use (x,y,f20,f02,dh20,dh02).


      open( newunit=tape_mu , file='./output/QFAM_output/mu.out' , status='unknown' );
      write(tape_mu,  '(a)'  ) 'Kernel Polynomial Method.';
      write(      6,'(/,a,/)') 'Kernel Polynomial Method.';

      call print_header(tape_mu);

      write(tape_mu,'(a,f10.3,a)') 'Omega_b = ' , Omega_b , ' [MeV].';
      write(      6,'(a,f10.3,a)') 'Omega_b = ' , Omega_b , ' [MeV].';

      write(tape_mu,'(/,8x,a,10x,a,/)') 'n' , 'mu_n';
      write(      6,'(/,8x,a,10x,a,/)') 'n' , 'mu_n';

      call flush(tape_mu);
      call flush(      6);


      ! |alpha_old> = [F20;-F02].
      do it = 1 , 2
          call axpbyBlockMatrix( +1.d0,f20(it) , 0.d0,alpha_old_x(it) );
          call axpbyBlockMatrix( -1.d0,f02(it) , 0.d0,alpha_old_y(it) );
      end do

      ! |alpha_new> = 1/Omega_b * [I,0;0,-I] * [A,B;B^*,A^*] * [F20;-F02].
      call applyMapping_XS( alpha_old_x,alpha_old_y , alpha_new_x,alpha_new_y );
      do it = 1 , 2
          call axpbyBlockMatrix( 0.d0,alpha_new_x(it) , 1/Omega_b,alpha_new_x(it) );
          call axpbyBlockMatrix( 0.d0,alpha_new_y(it) , 1/Omega_b,alpha_new_y(it) );
      end do

      ! mu_0 = 0.5 * [F20;F02]' * |alpha_old>.
      mu(0) = 0.5d0 * real( + sum([( scalarProductBlockMatrix(f20(it),alpha_old_x(it)) , it=1,2 )]) &
                            + sum([( scalarProductBlockMatrix(f02(it),alpha_old_y(it)) , it=1,2 )]) );

      ! mu_1 = [F20;F02]' * |alpha_new>.
      mu(1) = real( + sum([( scalarProductBlockMatrix(f20(it),alpha_new_x(it)) , it=1,2 )]) &
                    + sum([( scalarProductBlockMatrix(f02(it),alpha_new_y(it)) , it=1,2 )]) );

      write(tape_mu,'(i9,5x,e16.9)') 0 , mu(0);
      write(      6,'(i9,5x,e16.9)') 0 , mu(0);

      do n = 1 , N_it

          ! mu_{2*n-1} = 2<alpha_old|[I,0;0,-I]|alpha_new> - mu(1).
          mu(2*n-1) = 2*real( + sum([( scalarProductBlockMatrix(alpha_old_x(it),alpha_new_x(it)) , it=1,2 )])          &
                              - sum([( scalarProductBlockMatrix(alpha_old_y(it),alpha_new_y(it)) , it=1,2 )]) ) - mu(1);

          ! mu_{2*n} = 2<alpha_new|[I,0;0,-I]|alpha_new> - 2*mu(0).
          mu(2*n) = 2*real( + sum([( scalarProductBlockMatrix(alpha_new_x(it),alpha_new_x(it)) , it=1,2 )])            &
                            - sum([( scalarProductBlockMatrix(alpha_new_y(it),alpha_new_y(it)) , it=1,2 )]) ) - 2*mu(0);

          write(tape_mu,'(i9,5x,e16.9)') 2*n-1 , mu(2*n-1);
          write(      6,'(i9,5x,e16.9)') 2*n-1 , mu(2*n-1);

          write(tape_mu,'(i9,5x,e16.9)') 2*n   , mu(2*n);
          write(      6,'(i9,5x,e16.9)') 2*n   , mu(2*n);


          if( n == N_it ) exit;


          ! |alpha_tmp> = |alpha_new>.
          do it = 1 , 2
              call axpbyBlockMatrix( +1.d0,alpha_new_x(it) , 0.d0,alpha_tmp_x(it) );
              call axpbyBlockMatrix( +1.d0,alpha_new_y(it) , 0.d0,alpha_tmp_y(it) );
          end do

          ! |alpha_new> <-- 2/Omega_b * ( [I,0;0,-I] * [A,B;B^*,A^*] * |alpha_new> ) - |alpha_old>.
          call applyMapping_XS( alpha_tmp_x,alpha_tmp_y , alpha_new_x,alpha_new_y );
          do it = 1 , 2
              call axpbyBlockMatrix( -1.d0,alpha_old_x(it) , 2/Omega_b,alpha_new_x(it) );
              call axpbyBlockMatrix( -1.d0,alpha_old_y(it) , 2/Omega_b,alpha_new_y(it) );
          end do

          ! |alpha_old> = |alpha_tmp>.
          do it = 1 , 2
              call axpbyBlockMatrix( +1.d0,alpha_tmp_x(it) , 0.d0,alpha_old_x(it) );
              call axpbyBlockMatrix( +1.d0,alpha_tmp_y(it) , 0.d0,alpha_old_y(it) );
          end do

      end do


      close(tape_mu);


      return;

      contains

          subroutine axpbyBlockMatrix( a , X , b , Y )
              use dataTypes;
              implicit none;
              double precision         , intent(in)    :: a;
              double precision         , intent(in)    :: b;
              type(complexBlockMatrix) , intent(in)    :: X;
              type(complexBlockMatrix) , intent(inout) :: Y;
              integer                                  :: ib1, ib2;
              ! Performs Y <-- a*X + b*Y for block matrices X,Y and scalars a,b.

                  do ib2 = 1 , size(Y%nnzblocks,2)
                      do ib1 = 1 , size(Y%nnzblocks,1)
                          if( Y%nnzblocks(ib1,ib2) ) then

                              Y%blocks(ib1,ib2)%mat(:,:) = b * Y%blocks(ib1,ib2)%mat(:,:);

                              if( X%nnzblocks(ib1,ib2) ) &
                                  Y%blocks(ib1,ib2)%mat(:,:) = a * X%blocks(ib1,ib2)%mat(:,:) + Y%blocks(ib1,ib2)%mat(:,:);

                          end if
                      end do
                  end do

              return;
          end subroutine axpbyBlockMatrix

          double complex function scalarProductBlockMatrix( X , Y ) result(ans)
              use dataTypes;
              implicit none;
              type(complexBlockMatrix) , intent(in)    :: X;
              type(complexBlockMatrix) , intent(inout) :: Y;
              integer                                  :: ib1, ib2, i1, i2;
              ! Returns Tr[X'*Y].

                  ans = cmplx( 0.d0 , 0.d0 , kind=8 );
                  do ib2 = 1 , size(X%nnzblocks,2)
                      do ib1 = 1 , size(X%nnzblocks,1)
                          call assert( X%nnzblocks(ib1,ib2) .eqv. Y%nnzblocks(ib1,ib2) , 'nnzblocks error in scalarProductBlockMatrix.' );
                          if( X%nnzblocks(ib1,ib2) .and. Y%nnzblocks(ib1,ib2) ) then
                              do i2 = 1 , id_spx(ib2)
                                  do i1 = 1 , id_spx(ib1)
                                      ans = ans + conjg(X%blocks(ib1,ib2)%mat(i1,i2)) * Y%blocks(ib1,ib2)%mat(i1,i2);
                                  end do
                              end do
                          end if
                      end do
                  end do

              return;
          end function scalarProductBlockMatrix

          subroutine applyMapping_XS( xin,yin , xout,yout )
          use dataTypes;
          use simplex;
          use fam_energies;
          use dh02dh20matrix;
          use xyfam;
          implicit none;
          ! Calculates [ xout ; yout ] = [I,0;0,-I] * [A,B,B^*,A^*] * [ xin ; yin ].
          type(complexBlockMatrix) , dimension(2) , intent(in)    :: xin;
          type(complexBlockMatrix) , dimension(2) , intent(in)    :: yin;
          type(complexBlockMatrix) , dimension(2) , intent(inout) :: xout;
          type(complexBlockMatrix) , dimension(2) , intent(inout) :: yout;
          integer                                                 :: it, ib1, ib2, i1, i2;
          double precision                                        :: E1, E2;

          ! Inserting [xin;yin] in QFAM [x,y].
          do it = 1 , 2
              call axpbyBlockMatrix( 1.d0,xin(it) , 0.d0,x_fam(it) );
              call axpbyBlockMatrix( 1.d0,yin(it) , 0.d0,y_fam(it) );
          end do

          ! Perform QFAM iteration, i.e. calculate dH20, dH02 using [x,y] as input.
          call fam_drhodkappa();
          call fam_ddensdcurr();
          call fam_dpotentials();
          call fam_dh();
          call fam_ddelta();
          call fam_dh20dh02();

          ! Calculating ( [ A , B ; B^* , A^* ] * [X;Y] )_ij = [ dH20_ij + (E_i+E_j)X_ij ; dH02_ij +0 (Ei+Ej)Y_ij ].
          do it = 1 , 2
              do ib2 = 1 , N_blocks
                  do ib1 = 1 , N_blocks
                      if( dh20(it)%nnzblocks(ib1,ib2) .and. dh02(it)%nnzblocks(ib1,ib2) ) then
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)

                                  E1 = E_fam(it)%blocks(ib1)%vec(i1);
                                  E2 = E_fam(it)%blocks(ib2)%vec(i2);

                                  xout(it)%blocks(ib1,ib2)%mat(i1,i2) = dh20(it)%blocks(ib1,ib2)%mat(i1,i2) + (E1+E2) * xin(it)%blocks(ib1,ib2)%mat(i1,i2);
                                  yout(it)%blocks(ib1,ib2)%mat(i1,i2) = dh02(it)%blocks(ib1,ib2)%mat(i1,i2) + (E1+E2) * yin(it)%blocks(ib1,ib2)%mat(i1,i2);

                              end do
                          end do
                      end if
                  end do
              end do
          end do

          ! Perform multiplication with [I,0;0,-I].
          do it = 1 , 2
              call axpbyBlockMatrix( 0.d0,xout(it) , +1.d0,xout(it) );
              call axpbyBlockMatrix( 0.d0,yout(it) , -1.d0,yout(it) );
          end do

          return;
          end subroutine applyMapping_XS

      end subroutine kpm

