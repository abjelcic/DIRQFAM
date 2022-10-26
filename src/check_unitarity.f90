!======================================================================!

      subroutine check_unitarity()

!======================================================================!
      use dataTypes;
      use simplex;
      use u_matrix;
      use v_matrix;
      implicit none;

      double complex , parameter :: cone  = cmplx( 1.d0 , 0.d0 , kind=8 );
      double complex , parameter :: czero = cmplx( 0.d0 , 0.d0 , kind=8 );
      integer                    :: it, ib, i, j;
      double precision           :: error;
      type(complexBlockMatrix)   :: tmp;


      ! This subroutine checks that the Bogoliubov transformation matrix
      ! W = [ U , V^* ; V , U^* ] is unitary. See Eq. (7.3) from Ring-Schuck.
      ! This is equivalent to Eqs. (7.5) from Ring-Schuck:
      ! U'  * U   + V'  * V   = I ,
      ! U   * U'  + V^* * V^T = I ,
      ! U^T * V   + V^T * U   = 0 ,
      ! U   * V'  + V^* * U^T = 0 .
      !
      ! In simplex-y basis there holds:
      ! U = [ u , 0 ; 0 , u^* ] and V = [ 0 , -v^* ; v , 0 ].


      ! Allocating memory.
      allocate( tmp%blocks(1:N_blocks,1:N_blocks) );
      do ib = 1 , N_blocks
          allocate( tmp%blocks(ib,ib)%mat(1:id_spx(ib),1:id_spx(ib)) );
      end do


      ! U' * U  +  V' * V = I <==> tmp = u' * u + v' * v - I = 0.
      do it = 1 , 2

          do ib = 1 , N_blocks
              call zgemm( 'c' , 'n' , id_spx(ib) , id_spx(ib) , id_spx(ib)               , &
                          cone                                                           , &
                          u(it)%blocks(ib,ib)%mat(1,1) , size(u(it)%blocks(ib,ib)%mat,1) , &
                          u(it)%blocks(ib,ib)%mat(1,1) , size(u(it)%blocks(ib,ib)%mat,1) , &
                          czero                                                          , &
                          tmp%blocks(ib,ib)%mat(1,1)   , size(tmp%blocks(ib,ib)%mat,1)     );

              call zgemm( 'c' , 'n' , id_spx(ib) , id_spx(ib) , id_spx(ib)               , &
                          cone                                                           , &
                          v(it)%blocks(ib,ib)%mat(1,1) , size(v(it)%blocks(ib,ib)%mat,1) , &
                          v(it)%blocks(ib,ib)%mat(1,1) , size(v(it)%blocks(ib,ib)%mat,1) , &
                          cone                                                           , &
                          tmp%blocks(ib,ib)%mat(1,1)   , size(tmp%blocks(ib,ib)%mat,1)     );

              do i = 1 , id_spx(ib)
                  tmp%blocks(ib,ib)%mat(i,i) = tmp%blocks(ib,ib)%mat(i,i) - cone;
              end do
          end do

          error = -1.d0;
          do ib = 1 , N_blocks
              do j = 1 , id_spx(ib)
                  do i = 1 , id_spx(ib)
                      error = max( error , abs(tmp%blocks(ib,ib)%mat(i,j)) );
                  end do
              end do
          end do
          call assert( error < 1.d-8 , 'Bogoliubov transformation is not unitary!' );
#ifdef DEBUG
          write(*,'(a,e15.6)') 'Bogoliubov transformation unitarity error: ' , error;
#endif
      end do


      ! U * U'  +  V * V' = I <==> tmp = u * u' + v * v' - I = 0.
      do it = 1 , 2

          do ib = 1 , N_blocks
              call zgemm( 'n' , 'c' , id_spx(ib) , id_spx(ib) , id_spx(ib)               , &
                          cone                                                           , &
                          u(it)%blocks(ib,ib)%mat(1,1) , size(u(it)%blocks(ib,ib)%mat,1) , &
                          u(it)%blocks(ib,ib)%mat(1,1) , size(u(it)%blocks(ib,ib)%mat,1) , &
                          czero                                                          , &
                          tmp%blocks(ib,ib)%mat(1,1)   , size(tmp%blocks(ib,ib)%mat,1)     );

              call zgemm( 'n' , 'c' , id_spx(ib) , id_spx(ib) , id_spx(ib)               , &
                          cone                                                           , &
                          v(it)%blocks(ib,ib)%mat(1,1) , size(v(it)%blocks(ib,ib)%mat,1) , &
                          v(it)%blocks(ib,ib)%mat(1,1) , size(v(it)%blocks(ib,ib)%mat,1) , &
                          cone                                                           , &
                          tmp%blocks(ib,ib)%mat(1,1)   , size(tmp%blocks(ib,ib)%mat,1)     );

              do i = 1 , id_spx(ib)
                  tmp%blocks(ib,ib)%mat(i,i) = tmp%blocks(ib,ib)%mat(i,i) - cone;
              end do
          end do

          error = -1.d0;
          do ib = 1 , N_blocks
              do j = 1 , id_spx(ib)
                  do i = 1 , id_spx(ib)
                      error = max( error , abs(tmp%blocks(ib,ib)%mat(i,j)) );
                  end do
              end do
          end do
          call assert( error < 1.d-8 , 'Bogoliubov transformation is not unitary!' );
#ifdef DEBUG
          write(*,'(a,e15.6)') 'Bogoliubov transformation unitarity error: ' , error;
#endif
      end do


      ! U^T * V  +  V^T * U = 0 <==> tmp = u' * v - v' * u = 0.
      do it = 1 , 2

          do ib = 1 , N_blocks
              call zgemm( 'c' , 'n' , id_spx(ib) , id_spx(ib) , id_spx(ib)               , &
                          cone                                                           , &
                          u(it)%blocks(ib,ib)%mat(1,1) , size(u(it)%blocks(ib,ib)%mat,1) , &
                          v(it)%blocks(ib,ib)%mat(1,1) , size(v(it)%blocks(ib,ib)%mat,1) , &
                          czero                                                          , &
                          tmp%blocks(ib,ib)%mat(1,1)   , size(tmp%blocks(ib,ib)%mat,1)     );

              call zgemm( 'c' , 'n' , id_spx(ib) , id_spx(ib) , id_spx(ib)               , &
                          -cone                                                          , &
                          v(it)%blocks(ib,ib)%mat(1,1) , size(v(it)%blocks(ib,ib)%mat,1) , &
                          u(it)%blocks(ib,ib)%mat(1,1) , size(u(it)%blocks(ib,ib)%mat,1) , &
                          cone                                                           , &
                          tmp%blocks(ib,ib)%mat(1,1)   , size(tmp%blocks(ib,ib)%mat,1)     );
          end do

          error = -1.d0;
          do ib = 1 , N_blocks
              do j = 1 , id_spx(ib)
                  do i = 1 , id_spx(ib)
                      error = max( error , abs(tmp%blocks(ib,ib)%mat(i,j)) );
                  end do
              end do
          end do
          call assert( error < 1.d-8 , 'Bogoliubov transformation is not unitary!' );
#ifdef DEBUG
          write(*,'(a,e15.6)') 'Bogoliubov transformation unitarity error: ' , error;
#endif
      end do


      ! U * V'  +  V^* * U^T = 0 <==> tmp = u * v' - v * u' = 0.
      do it = 1 , 2

          do ib = 1 , N_blocks
              call zgemm( 'n' , 'c' , id_spx(ib) , id_spx(ib) , id_spx(ib)               , &
                          cone                                                           , &
                          u(it)%blocks(ib,ib)%mat(1,1) , size(u(it)%blocks(ib,ib)%mat,1) , &
                          v(it)%blocks(ib,ib)%mat(1,1) , size(v(it)%blocks(ib,ib)%mat,1) , &
                          czero                                                          , &
                          tmp%blocks(ib,ib)%mat(1,1)   , size(tmp%blocks(ib,ib)%mat,1)     );

              call zgemm( 'n' , 'c' , id_spx(ib) , id_spx(ib) , id_spx(ib)               , &
                          -cone                                                          , &
                          v(it)%blocks(ib,ib)%mat(1,1) , size(v(it)%blocks(ib,ib)%mat,1) , &
                          u(it)%blocks(ib,ib)%mat(1,1) , size(u(it)%blocks(ib,ib)%mat,1) , &
                          cone                                                           , &
                          tmp%blocks(ib,ib)%mat(1,1)   , size(tmp%blocks(ib,ib)%mat,1)     );
          end do

          error = -1.d0;
          do ib = 1 , N_blocks
              do j = 1 , id_spx(ib)
                  do i = 1 , id_spx(ib)
                      error = max( error , abs(tmp%blocks(ib,ib)%mat(i,j)) );
                  end do
              end do
          end do
          call assert( error < 1.d-8 , 'Bogoliubov transformation is not unitary!' );
#ifdef DEBUG
          write(*,'(a,e15.6)') 'Bogoliubov transformation unitarity error: ' , error;
#endif
      end do


      ! Deallocating memory.
      do ib = 1 , N_blocks
          deallocate( tmp%blocks(ib,ib)%mat );
      end do
      deallocate( tmp%blocks );


      return;
      end subroutine check_unitarity
