!======================================================================!

      double complex function fam_strength() result(S)

!======================================================================!
      use simplex;
      use f02f20matrix;
      use xyfam;
      implicit none;
      integer :: it, ib1, ib2, i1, i2;


          ! QRPA strength function is given by:
          ! S = sum_{mu<nu} conjg(F20_{mu,nu})*X_{mu,nu} + conjg(F02_{mu,nu})*Y_{mu,nu}.
          !
          ! When written in simplex-y basis using Eqs. (B.6) and (B.7) one obtains:
          ! S = Tr[f20'*x] + Tr[f02'*y].

          S = cmplx( 0.d0 , 0.d0 , kind=8 );
          do it = 1 , 2
              do ib2 = 1 , N_blocks
                  do ib1 = 1 , N_blocks

                      if( f20(it)%nnzblocks(ib1,ib2) .and. x_fam(it)%nnzblocks(ib1,ib2) ) then
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)
                                  S = S + conjg(f20(it)%blocks(ib1,ib2)%mat(i1,i2)) * x_fam(it)%blocks(ib1,ib2)%mat(i1,i2);
                              end do
                          end do
                      end if

                      if( f02(it)%nnzblocks(ib1,ib2) .and. y_fam(it)%nnzblocks(ib1,ib2) ) then
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)
                                  S = S + conjg(f02(it)%blocks(ib1,ib2)%mat(i1,i2)) * y_fam(it)%blocks(ib1,ib2)%mat(i1,i2);
                              end do
                          end do
                      end if

                  end do
              end do
          end do


      return;
      end function fam_strength
