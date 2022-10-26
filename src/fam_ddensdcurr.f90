!======================================================================!

      subroutine fam_ddensdcurr()

!======================================================================!
      use fam_input;
      use simplex;
      use PHI;
      use drho;
      use ddens;
      use dcurr;
      use clock_module;
      implicit none;

      double precision , external  :: frobeniousNormOfBlockMatrix;
      double precision , parameter :: pi = 4.d0*atan(1.d0);
      double complex   , parameter :: II = cmplx( 0.d0 , 1.d0 , kind=8 );
      integer                      :: it, ib1, ib2, i1, i2;
      character                    :: fg1, fg2;
      integer                      :: ml1, ml2;
      double precision             :: normrho;
      double complex               :: A1, A2;
      integer                      :: irun, Nruns;
      type(clockClass)             :: clock;




#ifdef DEBUG
      ! Selection rules test.
      ! (drho1+drho2)_{i1,i2} for fg1==fg2 are selected by |ml1-ml2|   = K.
      ! (drho1+drho2)_{i1,i2} for fg1/=fg2 are selected by |ml1+ml2+1| = K.
      ! (drho1-drho2)_{i1,i2} for fg1/=fg2 are selected by |ml1-ml2|   = K.
      do it = 1 , 2
          normrho = max( frobeniousNormOfBlockMatrix(drho1(it)) , &
                         frobeniousNormOfBlockMatrix(drho2(it))   );

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  call assert( drho1(it)%nnzblocks(ib1,ib2) .eqv. drho2(it)%nnzblocks(ib1,ib2) , 'drho1 and drho2 nnzblocks error!' );

                  if( drho1(it)%nnzblocks(ib1,ib2) .and. drho2(it)%nnzblocks(ib1,ib2) ) then
                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , id_spx(ib1)

                              fg1 = fg_spx(ib1)%index(i1);
                              ml1 = ml_spx(ib1)%index(i1);

                              fg2 = fg_spx(ib2)%index(i2);
                              ml2 = ml_spx(ib2)%index(i2);

                              if( fg1==fg2 .and. abs(ml1-ml2)/=K_multipole ) &
                                  call assert( abs( drho1(it)%blocks(ib1,ib2)%mat(i1,i2) + drho2(it)%blocks(ib1,ib2)%mat(i1,i2) ) / normrho < 1.d-12 , 'drho1+drho2 selection rules failed!' );

                              if( fg1/=fg2 .and. abs(ml1+ml2+1)/=K_multipole ) &
                                  call assert( abs( drho1(it)%blocks(ib1,ib2)%mat(i1,i2) + drho2(it)%blocks(ib1,ib2)%mat(i1,i2) ) / normrho < 1.d-12 , 'drho1+drho2 selection rules failed!' );

                              if( fg1/=fg2 .and. abs(ml1-ml2)/=K_multipole ) &
                                  call assert( abs( drho1(it)%blocks(ib1,ib2)%mat(i1,i2) - drho2(it)%blocks(ib1,ib2)%mat(i1,i2) ) / normrho < 1.d-12 , 'drho1-drho2 selection rules failed!' );

                          end do
                      end do
                  end if

              end do
          end do

      end do
#endif




      ! Calculation of drho_v.
      drho_v( -NGH:+NGH , 1:NGL , 1:2 ) = cmplx( 0.d0 , 0.d0 , kind=8 );
      do it = 1 , 2

          ! Construction of upper triangular block matrices Ar and Ai.
          do ib2 = 1 , N_blocks
              do ib1 = 1 , ib2
                  if( Ar_tmp%nnzblocks(ib1,ib2) .and. Ai_tmp%nnzblocks(ib1,ib2) ) then

                      Ar_tmp%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = 0.d0;
                      Ai_tmp%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = 0.d0;

                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , merge( i2 , id_spx(ib1) , ib1==ib2 )

                              fg1 = fg_spx(ib1)%index(i1);
                              ml1 = ml_spx(ib1)%index(i1);

                              fg2 = fg_spx(ib2)%index(i2);
                              ml2 = ml_spx(ib2)%index(i2);

                              if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then

                                  A1 = drho1(it)%blocks(ib1,ib2)%mat(i1,i2) + drho2(it)%blocks(ib1,ib2)%mat(i1,i2);
                                  A2 = drho1(it)%blocks(ib2,ib1)%mat(i2,i1) + drho2(it)%blocks(ib2,ib1)%mat(i2,i1);

                                  Ar_tmp%blocks(ib1,ib2)%mat(i1,i2) = merge( 0.5d0 , 1.d0 , all([ib1,i1]==[ib2,i2]) ) * real(A1+A2);
                                  Ai_tmp%blocks(ib1,ib2)%mat(i1,i2) = merge( 0.5d0 , 1.d0 , all([ib1,i1]==[ib2,i2]) ) * imag(A1+A2);

                              end if

                          end do
                      end do

                  end if
              end do
          end do

          ! Nruns = 1000;
          ! call clock%tic();
          ! do irun = 1 , Nruns
                call diagVtUtAUV( Ar_tmp , Ai_tmp , drho_v(:,:,it) );
          ! end do
          ! call clock%toc();
          ! write(6,'(a,f12.7,a)') 'Elapsed time per diagVtUtAUV: ' , clock%getTime()/Nruns , ' s.';

          drho_v(:,:,it) = 1/(2*pi) * drho_v(:,:,it);

      end do




      ! Calculation of drho_s. Since we don't need proton/neutron individually, we calculate the sum.
      drho_s( -NGH:+NGH , 1:NGL ) = cmplx( 0.d0 , 0.d0 , kind=8 );
      do it = 3 , 3

          ! Construction of upper triangular block matrices Ar and Ai.
          do ib2 = 1 , N_blocks
              do ib1 = 1 , ib2
                  if( Ar_tmp%nnzblocks(ib1,ib2) .and. Ai_tmp%nnzblocks(ib1,ib2) ) then

                      Ar_tmp%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = 0.d0;
                      Ai_tmp%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = 0.d0;

                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , merge( i2 , id_spx(ib1) , ib1==ib2 )

                              fg1 = fg_spx(ib1)%index(i1);
                              ml1 = ml_spx(ib1)%index(i1);

                              fg2 = fg_spx(ib2)%index(i2);
                              ml2 = ml_spx(ib2)%index(i2);

                              if( fg1==fg2 .and. abs(ml1-ml2)==K_multipole ) then

                                  A1 = + ( drho1(1)%blocks(ib1,ib2)%mat(i1,i2) + drho1(2)%blocks(ib1,ib2)%mat(i1,i2) ) &
                                       + ( drho2(1)%blocks(ib1,ib2)%mat(i1,i2) + drho2(2)%blocks(ib1,ib2)%mat(i1,i2) );

                                  A2 = + ( drho1(1)%blocks(ib2,ib1)%mat(i2,i1) + drho1(2)%blocks(ib2,ib1)%mat(i2,i1) ) &
                                       + ( drho2(1)%blocks(ib2,ib1)%mat(i2,i1) + drho2(2)%blocks(ib2,ib1)%mat(i2,i1) );

                                  A1 = merge( -1.d0 , +1.d0 , all([fg1,fg2]==['g','g']) ) * A1;
                                  A2 = merge( -1.d0 , +1.d0 , all([fg1,fg2]==['g','g']) ) * A2;

                                  Ar_tmp%blocks(ib1,ib2)%mat(i1,i2) = merge( 0.5d0 , 1.d0 , all([ib1,i1]==[ib2,i2]) ) * real(A1+A2);
                                  Ai_tmp%blocks(ib1,ib2)%mat(i1,i2) = merge( 0.5d0 , 1.d0 , all([ib1,i1]==[ib2,i2]) ) * imag(A1+A2);

                              end if

                          end do
                      end do

                  end if
              end do
          end do

          ! Nruns = 1000;
          ! call clock%tic();
          ! do irun = 1 , Nruns
                call diagVtUtAUV( Ar_tmp , Ai_tmp , drho_s(:,:) );
          ! end do
          ! call clock%toc();
          ! write(6,'(a,f12.7,a)') 'Elapsed time per diagVtUtAUV: ' , clock%getTime()/Nruns , ' s.';

          drho_s(:,:) = 1/(2*pi) * drho_s(:,:);

      end do




      ! Calculation of dj_1.
      dj_1( -NGH:+NGH , 1:NGL , 1:2 ) = cmplx( 0.d0 , 0.d0 , kind=8 );
      do it = 1 , 2

          ! Construction of upper triangular block matrices Ar and Ai.
          do ib2 = 1 , N_blocks
              do ib1 = 1 , ib2
                  if( Ar_tmp%nnzblocks(ib1,ib2) .and. Ai_tmp%nnzblocks(ib1,ib2) ) then

                      Ar_tmp%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = 0.d0;
                      Ai_tmp%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = 0.d0;

                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , merge( i2 , id_spx(ib1) , ib1==ib2 )

                              fg1 = fg_spx(ib1)%index(i1);
                              ml1 = ml_spx(ib1)%index(i1);

                              fg2 = fg_spx(ib2)%index(i2);
                              ml2 = ml_spx(ib2)%index(i2);

                              if( fg1/=fg2 .and. ml1+ml2+1==+K_multipole ) then

                                  A1 = drho1(it)%blocks(ib1,ib2)%mat(i1,i2) + drho2(it)%blocks(ib1,ib2)%mat(i1,i2);
                                  A2 = drho1(it)%blocks(ib2,ib1)%mat(i2,i1) + drho2(it)%blocks(ib2,ib1)%mat(i2,i1);

                                  A1 = merge( +1.d0 , -1.d0 , all([fg1,fg2]==['f','g']) ) * A1;
                                  A2 = merge( +1.d0 , -1.d0 , all([fg1,fg2]==['f','g']) ) * A2;

                                  Ar_tmp%blocks(ib1,ib2)%mat(i1,i2) = merge( 0.5d0 , 1.d0 , all([ib1,i1]==[ib2,i2]) ) * real(A1-A2);
                                  Ai_tmp%blocks(ib1,ib2)%mat(i1,i2) = merge( 0.5d0 , 1.d0 , all([ib1,i1]==[ib2,i2]) ) * imag(A1-A2);

                              end if

                          end do
                      end do

                  end if
              end do
          end do

          ! Nruns = 1000;
          ! call clock%tic();
          ! do irun = 1 , Nruns
                call diagVtUtAUV( Ar_tmp , Ai_tmp , dj_1(:,:,it) );
          ! end do
          ! call clock%toc();
          ! write(6,'(a,f12.7,a)') 'Elapsed time per diagVtUtAUV: ' , clock%getTime()/Nruns , ' s.';

          dj_1(:,:,it) = sqrt(2.d0)/(1+merge(1,0,K_multipole==0)) * 1/(2*pi*II) * dj_1(:,:,it);

      end do




      ! Calculation of dj_2.
      dj_2( -NGH:+NGH , 1:NGL , 1:2 ) = cmplx( 0.d0 , 0.d0 , kind=8 );
      do it = 1 , 2

          ! Construction of upper triangular block matrices Ar and Ai.
          do ib2 = 1 , N_blocks
              do ib1 = 1 , ib2
                  if( Ar_tmp%nnzblocks(ib1,ib2) .and. Ai_tmp%nnzblocks(ib1,ib2) ) then

                      Ar_tmp%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = 0.d0;
                      Ai_tmp%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = 0.d0;

                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , merge( i2 , id_spx(ib1) , ib1==ib2 )

                              fg1 = fg_spx(ib1)%index(i1);
                              ml1 = ml_spx(ib1)%index(i1);

                              fg2 = fg_spx(ib2)%index(i2);
                              ml2 = ml_spx(ib2)%index(i2);

                              if( fg1/=fg2 .and. ml1+ml2+1==-K_multipole ) then

                                  A1 = drho1(it)%blocks(ib1,ib2)%mat(i1,i2) + drho2(it)%blocks(ib1,ib2)%mat(i1,i2);
                                  A2 = drho1(it)%blocks(ib2,ib1)%mat(i2,i1) + drho2(it)%blocks(ib2,ib1)%mat(i2,i1);

                                  A1 = merge( +1.d0 , -1.d0 , all([fg1,fg2]==['f','g']) ) * A1;
                                  A2 = merge( +1.d0 , -1.d0 , all([fg1,fg2]==['f','g']) ) * A2;

                                  Ar_tmp%blocks(ib1,ib2)%mat(i1,i2) = merge( 0.5d0 , 1.d0 , all([ib1,i1]==[ib2,i2]) ) * real(A1-A2);
                                  Ai_tmp%blocks(ib1,ib2)%mat(i1,i2) = merge( 0.5d0 , 1.d0 , all([ib1,i1]==[ib2,i2]) ) * imag(A1-A2);

                              end if

                          end do
                      end do

                  end if
              end do
          end do

          ! Nruns = 1000;
          ! call clock%tic();
          ! do irun = 1 , Nruns
                call diagVtUtAUV( Ar_tmp , Ai_tmp , dj_2(:,:,it) );
          ! end do
          ! call clock%toc();
          ! write(6,'(a,f12.7,a)') 'Elapsed time per diagVtUtAUV: ' , clock%getTime()/Nruns , ' s.';

          dj_2(:,:,it) = sqrt(2.d0)/(1+merge(1,0,K_multipole==0)) * 1/(2*pi*II) * dj_2(:,:,it);

      end do




      ! Calculation of dj_3.
      dj_3( -NGH:+NGH , 1:NGL , 1:2 ) = cmplx( 0.d0 , 0.d0 , kind=8 );
      do it = 1 , 2

          ! Construction of upper triangular block matrices Ar and Ai.
          do ib2 = 1 , N_blocks
              do ib1 = 1 , ib2
                  if( Ar_tmp%nnzblocks(ib1,ib2) .and. Ai_tmp%nnzblocks(ib1,ib2) ) then

                      Ar_tmp%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = 0.d0;
                      Ai_tmp%blocks(ib1,ib2)%mat( 1:id_spx(ib1) , 1:id_spx(ib2) ) = 0.d0;

                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , merge( i2 , id_spx(ib1) , ib1==ib2 )

                              fg1 = fg_spx(ib1)%index(i1);
                              ml1 = ml_spx(ib1)%index(i1);

                              fg2 = fg_spx(ib2)%index(i2);
                              ml2 = ml_spx(ib2)%index(i2);

                              if( fg1/=fg2 .and. abs(ml1-ml2)==K_multipole ) then

                                  A1 = drho1(it)%blocks(ib1,ib2)%mat(i1,i2) - drho2(it)%blocks(ib1,ib2)%mat(i1,i2);
                                  A2 = drho1(it)%blocks(ib2,ib1)%mat(i2,i1) - drho2(it)%blocks(ib2,ib1)%mat(i2,i1);

                                  Ar_tmp%blocks(ib1,ib2)%mat(i1,i2) = merge( 0.5d0 , 1.d0 , all([ib1,i1]==[ib2,i2]) ) * real(A1+A2);
                                  Ai_tmp%blocks(ib1,ib2)%mat(i1,i2) = merge( 0.5d0 , 1.d0 , all([ib1,i1]==[ib2,i2]) ) * imag(A1+A2);

                              end if

                          end do
                      end do

                  end if
              end do
          end do

          ! Nruns = 1000;
          ! call clock%tic();
          ! do irun = 1 , Nruns
                call diagVtUtAUV( Ar_tmp , Ai_tmp , dj_3(:,:,it) );
          ! end do
          ! call clock%toc();
          ! write(6,'(a,f12.7,a)') 'Elapsed time per diagVtUtAUV: ' , clock%getTime()/Nruns , ' s.';

          dj_3(:,:,it) = 1/(2*pi) * dj_3(:,:,it);

      end do




      ! Calculation of dj_r, dj_p and dj_z.
      dj_r(:,:,:) = ( + dj_1(:,:,:) + dj_2(:,:,:) ) / sqrt(2.d0);
      dj_p(:,:,:) = ( - dj_1(:,:,:) + dj_2(:,:,:) ) / sqrt(2.d0);
      dj_z(:,:,:) =   + dj_3(:,:,:);




      ! Calculating Laplacians of the induced currents and densities.
      call fam_laplacian();




      return;
      end subroutine fam_ddensdcurr






!======================================================================!

      subroutine fam_laplacian()

!======================================================================!
      use fam_input;
      use quadrature;
      use basis;
      use ddens;
      use dcurr;
      use dlaplace;
      use famlaplacianmod;
      implicit none;

      logical          , save     :: first_call = .true.;
      integer                     :: i;
      integer                     :: Nz, Nr, K;
      integer                     :: ih, il;
      double precision            :: z, r;
      double precision , external :: phi_nz;
      double precision , external :: phi_nr_ml
      double complex              :: integral;
      double precision            :: bzz, bpp, fac;




      if( first_call ) then
          first_call = .false.;

          ! Initializing NZZ and NRR such that NZZ(i)+2*NRR(i) <= Nshells.
          i = 0;
          do Nr = 0 , Nshells/2
              do Nz = 0 , Nshells-2*Nr
                  i = i + 1;
                  NRR(i) = Nr;
                  NZZ(i) = Nz;
              end do
          end do
          call assert( i==size(NRR) .and. i==size(NZZ) .and. i==size(c) , 'Basis error in fam_laplacian.' );

          ! Initializing phiz.
          do Nz = 0 , Nshells
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;

                  z = sign(1,ih) * zb_fam(abs(ih));

                  phiz(ih,Nz) = phi_nz( Nz , bz/sqrt(2.d0) , z );

              end do
          end do

          ! Initializing phirK.
          do K = 0 , J_MAX+1
              do Nr = 0 , Nshells/2
                  do il = 1 , NGL

                      r = rb_fam(il);

                      phirK(il,Nr,K) = phi_nr_ml( Nr , K , bp/sqrt(2.d0) , r );

                  end do
              end do
          end do

      end if




      ! Calculation of Delta_{z,r,K} drho_v(z,r) for protons.
      c(:) = cmplx( 0.d0 , 0.d0 , kind=8 );
      do i = 1 , size(c)
          Nz = NZZ(i);
          Nr = NRR(i);

          integral = cmplx( 0.d0 , 0.d0 , kind=8 );
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;
                  integral = integral + wzwr(abs(ih),il) * drho_v(ih,il,2) * phiz(ih,Nz) * phirK(il,Nr,abs(K_multipole));
              end do
          end do
          c(i) = integral;

      end do

      ldrho_vprot = cmplx( 0.d0 , 0.d0 , kind=8 );
      do i = 1 , size(c)
          Nz = NZZ(i);
          Nr = NRR(i);

          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;

                  z   = sign(1,ih) * zb_fam(abs(ih));
                  r   = rb_fam(il);
                  bzz = bz/sqrt(2.d0);
                  bpp = bp/sqrt(2.d0);

                  fac = z**2/bzz**4 + r**2/bpp**4 - (2*Nz+1)/bzz**2 - 2*(2*Nr+abs(K_multipole)+1)/bpp**2;

                  ldrho_vprot(ih,il) = ldrho_vprot(ih,il) + c(i) * fac * phiz(ih,Nz) * phirK(il,Nr,abs(K_multipole));
              end do
          end do

      end do




      ! Calculation of Delta_{z,r,K} drho_s(z,r).
      c(:) = cmplx( 0.d0 , 0.d0 , kind=8 );
      do i = 1 , size(c)
          Nz = NZZ(i);
          Nr = NRR(i);

          integral = cmplx( 0.d0 , 0.d0 , kind=8 );
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;
                  integral = integral + wzwr(abs(ih),il) * drho_s(ih,il) * phiz(ih,Nz) * phirK(il,Nr,abs(K_multipole));
              end do
          end do
          c(i) = integral;

      end do

      ldrho_s = cmplx( 0.d0 , 0.d0 , kind=8 );
      do i = 1 , size(c)
          Nz = NZZ(i);
          Nr = NRR(i);

          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;

                  z   = sign(1,ih) * zb_fam(abs(ih));
                  r   = rb_fam(il);
                  bzz = bz/sqrt(2.d0);
                  bpp = bp/sqrt(2.d0);

                  fac = z**2/bzz**4 + r**2/bpp**4 - (2*Nz+1)/bzz**2 - 2*(2*Nr+abs(K_multipole)+1)/bpp**2;

                  ldrho_s(ih,il) = ldrho_s(ih,il) + c(i) * fac * phiz(ih,Nz) * phirK(il,Nr,abs(K_multipole));
              end do
          end do

      end do




      ! Calculation of Delta_{z,r,|K-1|} dj_1(z,r) for protons.
      c(:) = cmplx( 0.d0 , 0.d0 , kind=8 );
      do i = 1 , size(c)
          Nz = NZZ(i);
          Nr = NRR(i);

          integral = cmplx( 0.d0 , 0.d0 , kind=8 );
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;
                  integral = integral + wzwr(abs(ih),il) * dj_1(ih,il,2) * phiz(ih,Nz) * phirK(il,Nr,abs(K_multipole-1));
              end do
          end do
          c(i) = integral;

      end do

      ldj_1prot = cmplx( 0.d0 , 0.d0 , kind=8 );
      do i = 1 , size(c)
          Nz = NZZ(i);
          Nr = NRR(i);

          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;

                  z   = sign(1,ih) * zb_fam(abs(ih));
                  r   = rb_fam(il);
                  bzz = bz/sqrt(2.d0);
                  bpp = bp/sqrt(2.d0);

                  fac = z**2/bzz**4 + r**2/bpp**4 - (2*Nz+1)/bzz**2 - 2*(2*Nr+abs(K_multipole-1)+1)/bpp**2;

                  ldj_1prot(ih,il) = ldj_1prot(ih,il) + c(i) * fac * phiz(ih,Nz) * phirK(il,Nr,abs(K_multipole-1));
              end do
          end do

      end do




      ! Calculation of Delta_{z,r,|K+1|} dj_2(z,r) for protons.
      c(:) = cmplx( 0.d0 , 0.d0 , kind=8 );
      do i = 1 , size(c)
          Nz = NZZ(i);
          Nr = NRR(i);

          integral = cmplx( 0.d0 , 0.d0 , kind=8 );
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;
                  integral = integral + wzwr(abs(ih),il) * dj_2(ih,il,2) * phiz(ih,Nz) * phirK(il,Nr,abs(K_multipole+1));
              end do
          end do
          c(i) = integral;

      end do

      ldj_2prot = cmplx( 0.d0 , 0.d0 , kind=8 );
      do i = 1 , size(c)
          Nz = NZZ(i);
          Nr = NRR(i);

          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;

                  z   = sign(1,ih) * zb_fam(abs(ih));
                  r   = rb_fam(il);
                  bzz = bz/sqrt(2.d0);
                  bpp = bp/sqrt(2.d0);

                  fac = z**2/bzz**4 + r**2/bpp**4 - (2*Nz+1)/bzz**2 - 2*(2*Nr+abs(K_multipole+1)+1)/bpp**2;

                  ldj_2prot(ih,il) = ldj_2prot(ih,il) + c(i) * fac * phiz(ih,Nz) * phirK(il,Nr,abs(K_multipole+1));
              end do
          end do

      end do




      ! Calculation of Delta_{z,r,|K|} dj_3(z,r) for protons.
      c(:) = cmplx( 0.d0 , 0.d0 , kind=8 );
      do i = 1 , size(c)
          Nz = NZZ(i);
          Nr = NRR(i);

          integral = cmplx( 0.d0 , 0.d0 , kind=8 );
          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;
                  integral = integral + wzwr(abs(ih),il) * dj_3(ih,il,2) * phiz(ih,Nz) * phirK(il,Nr,abs(K_multipole));
              end do
          end do
          c(i) = integral;

      end do

      ldj_3prot = cmplx( 0.d0 , 0.d0 , kind=8 );
      do i = 1 , size(c)
          Nz = NZZ(i);
          Nr = NRR(i);

          do il = 1 , NGL
              do ih = -NGH , +NGH
                  if( ih == 0 ) cycle;

                  z   = sign(1,ih) * zb_fam(abs(ih));
                  r   = rb_fam(il);
                  bzz = bz/sqrt(2.d0);
                  bpp = bp/sqrt(2.d0);

                  fac = z**2/bzz**4 + r**2/bpp**4 - (2*Nz+1)/bzz**2 - 2*(2*Nr+abs(K_multipole)+1)/bpp**2;

                  ldj_3prot(ih,il) = ldj_3prot(ih,il) + c(i) * fac * phiz(ih,Nz) * phirK(il,Nr,abs(K_multipole));
              end do
          end do

      end do




      return;
      end subroutine fam_laplacian






!======================================================================!

      subroutine diagVtUtAUV( Ar , Ai , u )

!======================================================================!
      use dataTypes;
      use fam_input;
      use simplex;
      use PHI;
      implicit none;

      type(realBlockMatrix) ,                                  intent(in)    :: Ar;
      type(realBlockMatrix) ,                                  intent(in)    :: Ai;
      double complex        , dimension( -NGH:+NGH , 1:NGL ) , intent(inout) :: u;

      double precision      ,                                  parameter     :: one  = 1.d0;
      double precision      ,                                  parameter     :: zero = 0.d0;
      integer                                                                :: ib1, ib2;
      integer                                                                :: iG0, iG;
      integer                                                                :: Nresidual;
      integer                                                                :: ihl, ih, il;
      double precision      ,                                  external      :: ddot;




      u(:,:) = cmplx( 0.d0 , 0.d0 , kind=8 );




      ! Real part of u.
      ! Calculation of AU = A*U_PHI.
      AU(:,:) = 0.d0;
      do ib1 = 1 , N_blocks
          do ib2 = ib1 , N_blocks
              if( Ar%nnzblocks(ib1,ib2) ) then

                  call dgemm( 'n' , 'n' , id_spx(ib1) , size(U_PHI,2) , id_spx(ib2)        , &
                              one                                                          , &
                              Ar%blocks(ib1,ib2)%mat(1,1) , size(Ar%blocks(ib1,ib2)%mat,1) , &
                              U_PHI(ia_spx(ib2),1)        , size(U_PHI,1)                  , &
                              one                                                          , &
                              AU(ia_spx(ib1),1)           , size(AU,1)                       );

              end if
          end do
      end do

      ! Calculation of UtAU = U_PHI^T * A * U_PHI = U_PHI^T * (AU).
      call dgemm( 't' , 'n' , size(U_PHI,2) , size(AU,2) , N_total , &
                  one                                              , &
                  U_PHI(1,1) , size(U_PHI,1)                       , &
                  AU(1,1)    , size(AU,1)                          , &
                  zero                                             , &
                  UtAU(1,1)  , size(UtAU,1)                          );

      ! Calculation of VtUtAUV(iG) = V_PHI(:,iG)^T * ( U_PHI^T * A * U_PHI ) * V_PHI(:,iG) in batches.
      do iG0 = 1 , (size(V_PHI,2)/nbsize)*nbsize , nbsize

          call dgemm( 'n' , 'n' , size(UtAU,1) , nbsize , size(UtAU,2) , &
                      one                                              , &
                      UtAU(1,1)    , size(UtAU,1)                      , &
                      V_PHI(1,iG0) , size(V_PHI,1)                     , &
                      zero                                             , &
                      UtAUV(1,1)   , size(UtAUV,1)                       );

          do iG = iG0 , iG0+nbsize-1
              VtUtAUV(iG) = ddot( size(V_PHI,1) , UtAUV(1,iG-iG0+1),1 , V_PHI(1,iG),1  );
          end do

      end do

      Nresidual = size(V_PHI,2) - iG0 + 1;
      if( Nresidual > 0 ) then

          call dgemm( 'n' , 'n' , size(UtAU,1) , Nresidual , size(UtAU,2) , &
                      one                                                 , &
                      UtAU(1,1)    , size(UtAU,1)                         , &
                      V_PHI(1,iG0) , size(V_PHI,1)                        , &
                      zero                                                , &
                      UtAUV(1,1)   , size(UtAUV,1)                          );

          do iG = iG0 , size(V_PHI,2)
              VtUtAUV(iG) = ddot( size(V_PHI,1) , UtAUV(1,iG-iG0+1),1 , V_PHI(1,iG),1  );
          end do

      end if

      ihl = 0;
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih == 0 ) cycle;
              ihl = ihl + 1;

              u(ih,il) = u(ih,il) + cmplx( VtUtAUV(ihl) , 0.d0 , kind=8 );

          end do
      end do




      ! Imaginary part of u.
      ! Calculation of AU = A*U_PHI.
      AU(:,:) = 0.d0;
      do ib1 = 1 , N_blocks
          do ib2 = ib1 , N_blocks
              if( Ai%nnzblocks(ib1,ib2) ) then

                  call dgemm( 'n' , 'n' , id_spx(ib1) , size(U_PHI,2) , id_spx(ib2)        , &
                              one                                                          , &
                              Ai%blocks(ib1,ib2)%mat(1,1) , size(Ai%blocks(ib1,ib2)%mat,1) , &
                              U_PHI(ia_spx(ib2),1)        , size(U_PHI,1)                  , &
                              one                                                          , &
                              AU(ia_spx(ib1),1)           , size(AU,1)                       );

              end if
          end do
      end do

      ! Calculation of UtAU = U_PHI^T * A * U_PHI = U_PHI^T * (AU).
      call dgemm( 't' , 'n' , size(U_PHI,2) , size(AU,2) , N_total , &
                  one                                              , &
                  U_PHI(1,1) , size(U_PHI,1)                       , &
                  AU(1,1)    , size(AU,1)                          , &
                  zero                                             , &
                  UtAU(1,1)  , size(UtAU,1)                          );

      ! Calculation of VtUtAUV(iG) = V_PHI(:,iG)^T * ( U_PHI^T * A * U_PHI ) * V_PHI(:,iG) in batches.
      do iG0 = 1 , (size(V_PHI,2)/nbsize)*nbsize , nbsize

          call dgemm( 'n' , 'n' , size(UtAU,1) , nbsize , size(UtAU,2) , &
                      one                                              , &
                      UtAU(1,1)    , size(UtAU,1)                      , &
                      V_PHI(1,iG0) , size(V_PHI,1)                     , &
                      zero                                             , &
                      UtAUV(1,1)   , size(UtAUV,1)                       );

          do iG = iG0 , iG0+nbsize-1
              VtUtAUV(iG) = ddot( size(V_PHI,1) , UtAUV(1,iG-iG0+1),1 , V_PHI(1,iG),1  );
          end do

      end do

      Nresidual = size(V_PHI,2) - iG0 + 1;
      if( Nresidual > 0 ) then

          call dgemm( 'n' , 'n' , size(UtAU,1) , Nresidual , size(UtAU,2) , &
                      one                                                 , &
                      UtAU(1,1)    , size(UtAU,1)                         , &
                      V_PHI(1,iG0) , size(V_PHI,1)                        , &
                      zero                                                , &
                      UtAUV(1,1)   , size(UtAUV,1)                          );

          do iG = iG0 , size(V_PHI,2)
              VtUtAUV(iG) = ddot( size(V_PHI,1) , UtAUV(1,iG-iG0+1),1 , V_PHI(1,iG),1  );
          end do

      end if

      ihl = 0;
      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih == 0 ) cycle;
              ihl = ihl + 1;

              u(ih,il) = u(ih,il) + cmplx( 0.d0 , VtUtAUV(ihl) , kind=8 );

          end do
      end do




      return;
      end subroutine diagVtUtAUV
