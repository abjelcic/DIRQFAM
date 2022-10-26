!======================================================================!

      subroutine fam_xy( omega_gamma )

!======================================================================!
      use fam_input;
      use simplex;
      use fam_energies;
      use f02f20matrix;
      use dh02dh20matrix;
      use xyfam;
      implicit none;

      double complex   , intent(in) :: omega_gamma;
      integer                       :: it, ib1, ib2, i1, i2;
      double precision              :: E1, E2;


          ! Calculation of x(omega), see Eq. (B.16).
          do it = 1 , 2
              do ib2 = 1 , N_blocks
                  do ib1 = 1 , N_blocks
                      if( x_fam(it)%nnzblocks(ib1,ib2) ) then

                          x_fam(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = cmplx(0.d0,0.d0,kind=8);

                          if( f20(it)%nnzblocks(ib1,ib2) ) then
                              do i2 = 1 , id_spx(ib2)
                                  do i1 = 1 , id_spx(ib1)
                                      E1 = E_fam(it)%blocks(ib1)%vec(i1);
                                      E2 = E_fam(it)%blocks(ib2)%vec(i2);

                                      x_fam(it)%blocks(ib1,ib2)%mat(i1,i2) = x_fam(it)%blocks(ib1,ib2)%mat(i1,i2) + ( -f20(it)%blocks(ib1,ib2)%mat(i1,i2) / (E1+E2-omega_gamma) );

                                  end do
                              end do
                          end if

                          if( dh20(it)%nnzblocks(ib1,ib2) ) then
                              do i2 = 1 , id_spx(ib2)
                                  do i1 = 1 , id_spx(ib1)
                                      E1 = E_fam(it)%blocks(ib1)%vec(i1);
                                      E2 = E_fam(it)%blocks(ib2)%vec(i2);

                                      x_fam(it)%blocks(ib1,ib2)%mat(i1,i2) = x_fam(it)%blocks(ib1,ib2)%mat(i1,i2) + ( -dh20(it)%blocks(ib1,ib2)%mat(i1,i2) / (E1+E2-omega_gamma) );

                                  end do
                              end do
                          end if

                      end if
                  end do
              end do
          end do


          ! Calculation of y(omega), see Eq. (B.17).
          do it = 1 , 2
              do ib2 = 1 , N_blocks
                  do ib1 = 1 , N_blocks
                      if( y_fam(it)%nnzblocks(ib1,ib2) ) then

                          y_fam(it)%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = cmplx(0.d0,0.d0,kind=8);

                          if( f02(it)%nnzblocks(ib1,ib2) ) then
                              do i2 = 1 , id_spx(ib2)
                                  do i1 = 1 , id_spx(ib1)
                                      E1 = E_fam(it)%blocks(ib1)%vec(i1);
                                      E2 = E_fam(it)%blocks(ib2)%vec(i2);

                                      y_fam(it)%blocks(ib1,ib2)%mat(i1,i2) = y_fam(it)%blocks(ib1,ib2)%mat(i1,i2) + ( -f02(it)%blocks(ib1,ib2)%mat(i1,i2) / (E1+E2+omega_gamma) );

                                  end do
                              end do
                          end if

                          if( dh02(it)%nnzblocks(ib1,ib2) ) then
                              do i2 = 1 , id_spx(ib2)
                                  do i1 = 1 , id_spx(ib1)
                                      E1 = E_fam(it)%blocks(ib1)%vec(i1);
                                      E2 = E_fam(it)%blocks(ib2)%vec(i2);

                                      y_fam(it)%blocks(ib1,ib2)%mat(i1,i2) = y_fam(it)%blocks(ib1,ib2)%mat(i1,i2) + ( -dh02(it)%blocks(ib1,ib2)%mat(i1,i2) / (E1+E2+omega_gamma) );

                                  end do
                              end do
                          end if

                      end if
                  end do
              end do
          end do


      return;
      end subroutine fam_xy
