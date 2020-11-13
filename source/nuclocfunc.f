c======================================================================c

      subroutine nuclocfunc( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE simplex;
      USE v_matrix;
      USE basis;
      USE nnz_blocks;
      USE drho;
      USE quadrature;
      USE nuclearlocfunc;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      DOUBLE PRECISION f0r( -NGH:NGH , 1:NGL , 2 );
      DOUBLE COMPLEX   dfr( -NGH:NGH , 1:NGL , 2 );

      DOUBLE PRECISION f0z( -NGH:NGH , 1:NGL , 2 );
      DOUBLE COMPLEX   dfz( -NGH:NGH , 1:NGL , 2 );

      DOUBLE PRECISION g0 ( -NGH:NGH , 1:NGL , 2 );
      DOUBLE COMPLEX   dg ( -NGH:NGH , 1:NGL , 2 );

      DOUBLE PRECISION h0 ( -NGH:NGH , 1:NGL , 2 );
      DOUBLE COMPLEX   dh ( -NGH:NGH , 1:NGL , 2 );

      DOUBLE PRECISION F0 ( -NGH:NGH , 1:NGL , 2 );
      DOUBLE COMPLEX   dF ( -NGH:NGH , 1:NGL , 2 );

      DOUBLE PRECISION D0 ( -NGH:NGH , 1:NGL , 2 );
      DOUBLE COMPLEX   dD ( -NGH:NGH , 1:NGL , 2 );

      CHARACTER fg1, fg2;
      DOUBLE COMPLEX rho0( NBSX , NBSX , NBX , 2 );
      DOUBLE COMPLEX z;
      pi = 3.14159265358979324D0;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN nuclocfunc() **************************';
      write(6,*) '';
      endif






c-----Calculation of rho0 matrix
      do it = 1 , 2
        do ib = 1 , N_blocks
          call zgemm( 'n'        ,   'c'        ,
     &                id_spx(ib) ,   id_spx(ib) , id_spx(ib) ,
     &                COMPLEX( 1.D0 , 0.D0 )    ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                v(1,1,ib,it)              , NBSX       ,
     &                COMPLEX( 0.D0 , 0.D0 )    ,
     &                rho0(1,1,ib,it)           , NBSX         );
        enddo
      enddo






c-----Calculation of f0r, f0z, g0 and h0
      f0r = 0.D0;
      f0z = 0.D0;
      g0  = 0.D0;
      h0  = 0.D0;
      do it = 1 , 2
        do il = 1 , NGL
          do ih = -NGH , +NGH
            if( ih .eq. 0 ) CYCLE;


            do ib = 1 , N_blocks
              do i = 1 , id_spx(ib)
                fg1 = fg_spx(i,ib);
                ml1 = ml_spx(i,ib);
                do j = 1 , id_spx(ib)
                  fg2 = fg_spx(j,ib);
                  ml2 = ml_spx(j,ib);

                  if(         fg1  .ne. fg2 ) CYCLE;
                  if( abs(ml1-ml2) .ne.  0  ) CYCLE;

                  ii = i-1+ia_spx(ib);
                  jj = j-1+ia_spx(ib);


                  r = rb_fam(il);
                  call assert( r.gt.1.D-6 , 'r close to zero' );

                   phiz1 =  phi_z(ih,ii)
                  dphiz1 = dphi_z(ih,ii);
                   phir1 =  phi_r(il,ii);
                  dphir1 = dphi_r(il,ii);

                   phiz2 =  phi_z(ih,jj)
                  dphiz2 = dphi_z(ih,jj);
                   phir2 =  phi_r(il,jj);
                  dphir2 = dphi_r(il,jj);

                  x = DREAL(rho0(i,j,ib,it));
                  x = x * phiz1*phiz2;
                  x = x * 0.5D0 * ( dphir1*phir2 + phir1*dphir2 );
                  x = x / (2.D0*pi);
                  f0r(ih,il,it) = f0r(ih,il,it) + x;

                  x = DREAL(rho0(i,j,ib,it));
                  x = x * 0.5D0 * ( dphiz1*phiz2 + phiz1*dphiz2 );
                  x = x * phir1*phir2;
                  x = x / (2.D0*pi);
                  f0z(ih,il,it) = f0z(ih,il,it) + x;

                  x = DREAL(rho0(i,j,ib,it));
                  x = x * phiz1*phiz2;
                  x = x * phir1*phir2;
                  x = x / (2.D0*pi);
                  g0(ih,il,it) = g0(ih,il,it) + x;

                  x = DBLE(ml1*ml2)/(r*r);
                  x = x *  phiz1*phiz2  *  phir1*phir2;
                  x = x + dphiz1*dphiz2 *  phir1*phir2;
                  x = x +  phiz1*phiz2  * dphir1*dphir2;
                  x = x / (2.D0*pi);
                  x = x * DREAL(rho0(i,j,ib,it));
                  h0(ih,il,it) = h0(ih,il,it) + x;


                enddo
              enddo
            enddo


          enddo
        enddo
      enddo






c-----Calculation of dfr, dfz, dg and dh
      dfr = COMPLEX( 0.D0 , 0.D0 );
      dfz = COMPLEX( 0.D0 , 0.D0 );
      dg  = COMPLEX( 0.D0 , 0.D0 );
      dh  = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2
        do il = 1 , NGL
          do ih = -NGH , +NGH
            if( ih .eq. 0 ) CYCLE;


            do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                ! Recall that drho_1/2 have nnz pattern as dh_nnz
                if( dh_nnz(ib1,ib2) .eqv. .false. ) CYCLE;

                j0 = ia_spx(ib2);
                j1 = ia_spx(ib2)+id_spx(ib2)-1;
                do j = j0 , j1
                  fg2 = fg_spx(j-j0+1,ib2);
                  ml2 = ml_spx(j-j0+1,ib2);

                  i0 = ia_spx(ib1);
                  i1 = ia_spx(ib1)+id_spx(ib1)-1;
                  do i = i0 , i1
                    fg1 = fg_spx(i-i0+1,ib1);
                    ml1 = ml_spx(i-i0+1,ib1);

                    if(         fg1  .ne.         fg2 ) CYCLE;
                    if( abs(ml1-ml2) .ne. K_multipole ) CYCLE;


                    r = rb_fam(il);
                    call assert( r.gt.1.D-6 , 'r close to zero' );

                     phiz1 =  phi_z(ih,i)
                    dphiz1 = dphi_z(ih,i);
                     phir1 =  phi_r(il,i);
                    dphir1 = dphi_r(il,i);

                     phiz2 =  phi_z(ih,j)
                    dphiz2 = dphi_z(ih,j);
                     phir2 =  phi_r(il,j);
                    dphir2 = dphi_r(il,j);

                    x = phiz1*phiz2;
                    x = x * 0.5D0 * ( dphir1*phir2 + phir1*dphir2 );
                    x = x / (2.D0*pi);
                    z = x * 0.5D0 * ( drho_1(i,j,it) + drho_2(i,j,it) );
                    dfr(ih,il,it) = dfr(ih,il,it) + z;

                    x = 0.5D0 * ( dphiz1*phiz2 + phiz1*dphiz2 );
                    x = x * phir1*phir2;
                    x = x / (2.D0*pi);
                    z = x * 0.5D0 * ( drho_1(i,j,it) + drho_2(i,j,it) );
                    dfz(ih,il,it) = dfz(ih,il,it) + z;

                    x = phiz1*phiz2;
                    x = x * phir1*phir2;
                    x = x / (2.D0*pi);
                    z = x * 0.5D0 * ( drho_1(i,j,it) + drho_2(i,j,it) );
                    dg(ih,il,it) = dg(ih,il,it) + z;

                    x = DBLE(ml1*ml2)/(r*r);
                    x = x *  phiz1*phiz2  *  phir1*phir2;
                    x = x + dphiz1*dphiz2 *  phir1*phir2;
                    x = x +  phiz1*phiz2  * dphir1*dphir2;
                    x = x / (2.D0*pi);
                    z = x * 0.5D0 * ( drho_1(i,j,it) + drho_2(i,j,it) );
                    dh(ih,il,it) = dh(ih,il,it) + z;


                  enddo
                enddo
              enddo
            enddo


          enddo
        enddo
      enddo






c-----Calculation of F0 and dF
      do it = 1 , 2
        do il = 1 , NGL
          do ih = -NGH , +NGH
            if( ih .eq. 0 ) CYCLE;


            x            = f0r(ih,il,it)**2.D0 + f0z(ih,il,it)**2.D0;
            x            = x / g0(ih,il,it);
            x            = h0(ih,il,it) - x;
            F0(ih,il,it) = x;


            x            = f0r(ih,il,it)**2.D0 + f0z(ih,il,it)**2.D0;
            x            = x / ( g0(ih,il,it)**2.D0 );
            z            = COMPLEX( 0.D0 , 0.D0 );
            z            = z + f0r(ih,il,it) * dfr(ih,il,it)
            z            = z + f0z(ih,il,it) * dfz(ih,il,it);
            z            = z * ( -2.D0 / g0(ih,il,it) );
            z            = z + dh(ih,il,it);
            z            = z + x * dg(ih,il,it);
            dF(ih,il,it) = z;


          enddo
        enddo
      enddo






c-----Calculation of D0 and dD
      do it = 1 , 2
        do il = 1 , NGL
          do ih = -NGH , +NGH
            if( ih .eq. 0 ) CYCLE;


            fac          = (3.D0/5.D0) * (6.D0*pi*pi)**(2.D0/3.D0);
            tauTF        = fac * ( g0(ih,il,it) )**(5.D0/3.D0);
            D0(ih,il,it) = F0(ih,il,it) / tauTF;


            z            = dg(ih,il,it) / g0(ih,il,it);
            z            = z * (-5.D0/3.D0) * F0(ih,il,it) / tauTF;
            z            = z + dF(ih,il,it) / tauTF;
            dD(ih,il,it) = z;


          enddo
        enddo
      enddo






c-----Calculation of C0 and dC
      do it = 1 , 2
        do il = 1 , NGL
          do ih = -NGH , +NGH
            if( ih .eq. 0 ) CYCLE;


            C0(ih,il,it) = 1.D0 / ( 1.D0 + D0(ih,il,it)**2.D0 );


            x            = -2.D0 * D0(ih,il,it);
            x            = x / ( 1.D0 + D0(ih,il,it)**2.D0 )**2.D0;
            dC(ih,il,it) = x * dD(ih,il,it);


          enddo
        enddo
      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END nuclocfunc() ****************************';
      write(6,*) '';
      endif

      return;
      end;
