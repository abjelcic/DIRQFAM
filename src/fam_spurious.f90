!======================================================================!

      subroutine fam_spurious()

!======================================================================!
      use fam_input;
      use simplex;
      use xyfam;
      use spurious;
      implicit none;
      double complex :: lambdaR;
      double complex :: lambdaP;
      integer        :: it, ib1, ib2, i1, i2;


          ! Translational spurious mode removal only if J is odd and K=0 or K=1.
          if( .not.( mod(J_multipole,2)==1 .and. any(K_multipole==[0,1]) ) ) then
              return;
          end if


          ! lambdaR = ( - Tr[p20'*x] + Tr[p20'*y^T] ) / <Phi|[Rcm,Pcm]|Phi>.
          lambdaR = cmplx( 0.d0 , 0.d0 , kind=8 );
          do it = 1 , 2
              do ib2 = 1 , N_blocks
                  do ib1 = 1 , N_blocks

                      if( p20(it)%nnzblocks(ib1,ib2) .and. x_fam(it)%nnzblocks(ib1,ib2) ) then
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)
                                  lambdaR = lambdaR - conjg(p20(it)%blocks(ib1,ib2)%mat(i1,i2)) * x_fam(it)%blocks(ib1,ib2)%mat(i1,i2);
                              end do
                          end do
                      end if

                      if( p20(it)%nnzblocks(ib1,ib2) .and. y_fam(it)%nnzblocks(ib2,ib1) ) then ! Notice here (ib2,ib1).
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)
                                  lambdaR = lambdaR + conjg(p20(it)%blocks(ib1,ib2)%mat(i1,i2)) * y_fam(it)%blocks(ib2,ib1)%mat(i2,i1); ! Notice here (ib2,ib1) and (i2,i1).
                              end do
                          end do
                      end if

                  end do
              end do
          end do
          lambdaR = lambdaR / sum(RcmPcm_commutator(1:2));


          ! lambdaP = ( + Tr[r20'*x] + Tr[r20'*y^T] ) / <Phi|[Rcm,Pcm]|Phi>.
          lambdaP = cmplx( 0.d0 , 0.d0 , kind=8 );
          do it = 1 , 2
              do ib2 = 1 , N_blocks
                  do ib1 = 1 , N_blocks

                      if( r20(it)%nnzblocks(ib1,ib2) .and. x_fam(it)%nnzblocks(ib1,ib2) ) then
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)
                                  lambdaP = lambdaP + conjg(r20(it)%blocks(ib1,ib2)%mat(i1,i2)) * x_fam(it)%blocks(ib1,ib2)%mat(i1,i2);
                              end do
                          end do
                      end if

                      if( r20(it)%nnzblocks(ib1,ib2) .and. y_fam(it)%nnzblocks(ib2,ib1) ) then ! Notice here (ib2,ib1).
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)
                                  lambdaP = lambdaP + conjg(r20(it)%blocks(ib1,ib2)%mat(i1,i2)) * y_fam(it)%blocks(ib2,ib1)%mat(i2,i1); ! Notice here (ib2,ib1) and (i2,i1).
                              end do
                          end do
                      end if

                  end do
              end do
          end do
          lambdaP = lambdaP / sum(RcmPcm_commutator(1:2));


          ! Correction: x <-- x - lambdaR*r20 - lambdaP*p20.
          do it = 1 , 2
              do ib2 = 1 , N_blocks
                  do ib1 = 1 , N_blocks

                      if( r20(it)%nnzblocks(ib1,ib2) ) then
                          call assert( x_fam(it)%nnzblocks(ib1,ib2) , 'x error in fam_spurious.' );
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)
                                  x_fam(it)%blocks(ib1,ib2)%mat(i1,i2) = x_fam(it)%blocks(ib1,ib2)%mat(i1,i2) - lambdaR * r20(it)%blocks(ib1,ib2)%mat(i1,i2);
                              end do
                          end do
                      end if

                      if( p20(it)%nnzblocks(ib1,ib2) ) then
                          call assert( x_fam(it)%nnzblocks(ib1,ib2) , 'x error in fam_spurious.' );
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)
                                  x_fam(it)%blocks(ib1,ib2)%mat(i1,i2) = x_fam(it)%blocks(ib1,ib2)%mat(i1,i2) - lambdaP * p20(it)%blocks(ib1,ib2)%mat(i1,i2);
                              end do
                          end do
                      end if

                  end do
              end do
          end do


          ! Correction: y <-- y + lambdaR*conjg(r20) + lambdaP*conjg(p20).
          do it = 1 , 2
              do ib2 = 1 , N_blocks
                  do ib1 = 1 , N_blocks

                      if( r20(it)%nnzblocks(ib1,ib2) ) then
                          call assert( y_fam(it)%nnzblocks(ib1,ib2) , 'y error in fam_spurious.' );
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)
                                  y_fam(it)%blocks(ib1,ib2)%mat(i1,i2) = y_fam(it)%blocks(ib1,ib2)%mat(i1,i2) + lambdaR * conjg(r20(it)%blocks(ib1,ib2)%mat(i1,i2));
                              end do
                          end do
                      end if

                      if( p20(it)%nnzblocks(ib1,ib2) ) then
                          call assert( y_fam(it)%nnzblocks(ib1,ib2) , 'y error in fam_spurious.' );
                          do i2 = 1 , id_spx(ib2)
                              do i1 = 1 , id_spx(ib1)
                                  y_fam(it)%blocks(ib1,ib2)%mat(i1,i2) = y_fam(it)%blocks(ib1,ib2)%mat(i1,i2) + lambdaP * conjg(p20(it)%blocks(ib1,ib2)%mat(i1,i2));
                              end do
                          end do
                      end if

                  end do
              end do
          end do


      return;
      end subroutine fam_spurious
