!======================================================================!

      subroutine init_coulomb()

!======================================================================!
      use fam_input;
      use quadrature;
      use fam_green;
      implicit none;

      integer                      :: ih1, il1, ih2, il2;
      double precision             :: z1, r1, z2, r2, a;
      double precision , parameter :: hbc   = 197.328284d0;
      double precision , parameter :: alpha = 1.d0 / 137.03602d0; ! e^2/(4*pi*eps_0*hbar*c).
      double precision , external  :: I_K;


      G1( -NGH:+NGH , 1:NGL , -NGH:+NGH , 1:NGL ) = 0.d0;
      G2( -NGH:+NGH , 1:NGL , -NGH:+NGH , 1:NGL ) = 0.d0;
      G3( -NGH:+NGH , 1:NGL , -NGH:+NGH , 1:NGL ) = 0.d0;

      if( include_coulomb == 0 ) then
          return;
      end if

      do il1 = 1 , NGL
          do ih1 = -NGH , +NGH
              if( ih1 == 0 ) cycle;

              do il2 = 1 , NGL
                  do ih2 = -NGH , +NGH
                      if( ih2 == 0 ) cycle;

                      z1 = sign(1,ih1) * zb_fam(abs(ih1)); ! z
                      z2 = sign(1,ih2) * zb_fam(abs(ih2)); ! z'

                      r1 = rb_fam(il1);                    ! r
                      r2 = rb_fam(il2);                    ! r'

                      a = 4*r1*r2 / ( (r1+r2)**2 + (z1-z2)**2 );

                    ! Gi( z', r',  z,  r) = quadrature weight  *  2*e^2/(4*pi*eps_0) * [see Eq. (D.15)].
                      G1(ih2,il2,ih1,il1) = wzwr(abs(ih2),il2) * (2.d0*hbc*alpha)    * sqrt( (r1+r2)**2 + (z1-z2)**2 ) * I_K( a , abs(K_multipole-1) );
                      G2(ih2,il2,ih1,il1) = wzwr(abs(ih2),il2) * (2.d0*hbc*alpha)    * sqrt( (r1+r2)**2 + (z1-z2)**2 ) * I_K( a , abs(K_multipole+1) );
                      G3(ih2,il2,ih1,il1) = wzwr(abs(ih2),il2) * (2.d0*hbc*alpha)    * sqrt( (r1+r2)**2 + (z1-z2)**2 ) * I_K( a , abs(K_multipole+0) );
                  end do
              end do

          end do
      end do


      return;
      end subroutine init_coulomb
