!======================================================================!

      double precision function fam_gmres( iter , omega_gamma ) result(relResError)

!======================================================================!
      use dh02dh20matrix;
      use dh;
      use dDelta;
      use gmres;
      implicit none;
      integer          ,                intent(in)  :: iter;
      double complex   ,                intent(in)  :: omega_gamma;
      integer                                       :: it, j;
      double complex   , dimension(:) , allocatable :: WORK;
      double complex   , dimension(1)               :: WORKquery;
      integer                                       :: INFO;
      double precision                              :: relResError_check;

      call assert( iter <= maxIter , 'fam_gmres called with too large iter!' )




      ! We are solving x(w) = T(w)*x(w) + b(w), where x(w) stores the non-redundant elements of x(w) = [ dh1(w) ; dDelta1_pl(w) ; conj(dDelta1_mi(w)) ].
      ! In i-th iteration, for given x_i, the full QFAM iteration actually calculates x_{i+1} = T(w)*x_i + b(w).
      ! Vector b(w) is obtained by setting x=0, and performing the full QFAM iteration, i.e. it corresponds to the free response where dh20 = dh02 = 0.
      ! Thus we first calculate the vector b(w), and then solve the linear system A(w)*x(w) = (I-T(w))*x(w) = b(w), where we only have access to the mapping
      ! x_i |--> x_{i+1} = T(w)*x_i + b(w).
      ! In the first step, we need to find the initial residual r_0 = b - A*x_0. Since we are using x_0 = 0 for initial guess, we have r_0 = b(w).
      if( iter == 0 ) then


          ! Calculation of vectors b(w).
          do it = 1 , 2
              call setToZeroBlockMatrix( dh20(it) );
              call setToZeroBlockMatrix( dh02(it) );
          end do
          call fam_xy( omega_gamma );
          call fam_drhodkappa();
          call fam_ddensdcurr();
          call fam_dpotentials();
          call fam_dh();
          call fam_ddelta();
          ! Subroutines insertIntoGMRESvector and extractFromGMRESvector are implemented in gmres module.
          call insertIntoGMRESvector( dh1 , dh2 , dDelta1_pl , dDelta1_mi , b  );


          ! Initialize the GMRES solver and the first Arnoldi vector Q(:,1).
          r0(:)   = b(:);
          Q(:,1)  = r0(:) / norm2cmplx(r0(:));
          H(:,:)  = cmplx( 0.d0 , 0.d0 , kind=8 );
          beta(:) = cmplx( 0.d0 , 0.d0 , kind=8 );
          beta(1) = norm2cmplx(r0(:));


          ! Enter the first Arnoldi vector Q(:,1) into dh1, dh2, dDelta1_pl and dDelta1_mi.
          call extractFromGMRESvector( Q(:,1) , dh1 , dh2 , dDelta1_pl , dDelta1_mi );


          ! Return the relative residual error.
          relResError = norm2cmplx(r0(:)) / norm2cmplx(b(:));


          return;
      end if




      ! If iter>0, we have just run the QFAM iteration for the input vector q = Q(:,iter), and now for the next GMRES step we need
      ! the vector: v = A(w)*q = (I-T(w))*q = q - T(w)*q = q - (T(w)*q+b(w)) + b(w).
      ! As stated in the comment above, notice that the vector T(w)*q+b(w) is contained in the non-redundant elements of [ dh1(w) ; dDelta1_pl(w) ; conj(dDelta1_mi(w)) ]
      ! at the end of the iter-th iteration because we supplied the vector q at the start of the iter-th iteration.
      ! Therefore, we first store the non-redundant elements of [ dh1(w) ; dDelta1_pl(w) ; conj(dDelta1_mi(w)) ] into the vector v, and then set v=q-v+b.
      call insertIntoGMRESvector( dh1 , dh2 , dDelta1_pl , dDelta1_mi , v );
      v(:) = Q(:,iter) - v(:) + b(:);


      ! Arnoldi iteration step.
      do j = 1 , iter
          H(j,iter) = dot_product( Q(:,j) , v(:) );
          v(:)      = v(:) - H(j,iter)*Q(:,j);
      end do
      H( iter+1 , iter   ) = norm2cmplx(v(:));
      Q(   :    , iter+1 ) = v(:) / norm2cmplx(v(:));


      ! Enter the next Arnoldi vector Q(:,iter+1) into dh1, dh2, dDelta1_pl and dDelta2_mi.
      call extractFromGMRESvector( Q(:,iter+1) , dh1 , dh2 , dDelta1_pl , dDelta1_mi );


      ! GMRES step - minimize the residual norm || beta - H*y ||.
      H_tmp   ( 1:iter+1 , 1:iter ) = H   ( 1:iter+1 , 1:iter );
      beta_tmp( 1:iter+1          ) = beta( 1:iter+1          );

      call zgels( 'N'        ,                  &
                   iter+1    , iter           , &
                   1         ,                  &
                   H_tmp     , size(H_tmp,1)  , &
                   beta_tmp  , size(beta_tmp) , &
                   WORKquery , -1             , &
                   INFO                         );
      call assert( INFO==0 , 'zgels workspace query failed!' );
      allocate( WORK( 1:int(real(WORKquery(1)+0.5d0)) ) );

      call zgels( 'N'        ,                  &
                   iter+1    , iter           , &
                   1         ,                  &
                   H_tmp     , size(H_tmp,1)  , &
                   beta_tmp  , size(beta_tmp) , &
                   WORK      , size(WORK)     , &
                   INFO                         );
      call assert( INFO==0 , 'zgels failed!' );
      deallocate( WORK );

      y(1:iter)   = beta_tmp(1:iter);
      relResError = abs(beta_tmp(iter+1)) / norm2cmplx(b(:));

#if DEBUG
      relResError_check = norm2cmplx( beta(1:iter+1) - matmul(H(1:iter+1,1:iter),y(1:iter)) ) / norm2cmplx(b(:));
      call assert( abs(relResError_check-relResError)<1.d-10 , 'Relative residual error check failed!' );
#endif




      return;

      contains

          pure double precision function norm2cmplx(v) result(ans)
          implicit none;
          double complex , dimension(:) , intent(in) :: v;

              ans = sqrt(sum(abs(v(:))**2));

          return;
          end function norm2cmplx

      end function fam_gmres






!======================================================================!

      subroutine fam_gmres_getSolution( iter )

!======================================================================!
      use dh;
      use dDelta;
      use gmres;
      implicit none;
      integer        , intent(in) :: iter;
      double complex , parameter  :: cone  = cmplx( 1.d0 , 0.d0 , kind=8 );
      double complex , parameter  :: czero = cmplx( 0.d0 , 0.d0 , kind=8 );

      ! After the iter-th GMRES step, the solution is given by:
      ! x_solution(:) = x_0(:) + Q(:,1:iter)*y(1:iter) = Q(:,1:iter)*y(1:iter).
      ! Here we calculate x_solution and using it construct the solution of dh1, dh2, dDelta1_pl and dDelta1_mi.
      call zgemv( 'N' , xsize , iter     , &
                  cone       ,             &
                  Q          , size(Q,1) , &
                  y          , 1         , &
                  czero      ,             &
                  x_solution , 1           );

      call extractFromGMRESvector( x_solution , dh1 , dh2 , dDelta1_pl , dDelta1_mi );

      return;
      end subroutine fam_gmres_getSolution
