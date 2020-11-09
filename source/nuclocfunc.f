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



      DOUBLE PRECISION rhoAVG0            ( -NGH:NGH , 1:NGL , 2 );
      DOUBLE COMPLEX   drhoAVG            ( -NGH:NGH , 1:NGL , 2 );

      DOUBLE PRECISION tauAVG0            ( -NGH:NGH , 1:NGL , 2 );
      DOUBLE COMPLEX   dtauAVG            ( -NGH:NGH , 1:NGL , 2 );

      DOUBLE PRECISION grad2_rhoAVG0      ( -NGH:NGH , 1:NGL , 2 );
      DOUBLE COMPLEX   dot_rhoAVG0_drhoAVG( -NGH:NGH , 1:NGL , 2 );

      DOUBLE PRECISION dr_rhoAVG0         ( -NGH:NGH , 1:NGL , 2 );
      DOUBLE COMPLEX   dr_drhoAVG         ( -NGH:NGH , 1:NGL , 2 );

      DOUBLE PRECISION dz_rhoAVG0         ( -NGH:NGH , 1:NGL , 2 );
      DOUBLE COMPLEX   dz_drhoAVG         ( -NGH:NGH , 1:NGL , 2 );

      DOUBLE PRECISION f0                 ( -NGH:NGH , 1:NGL , 2 );
      DOUBLE COMPLEX   df                 ( -NGH:NGH , 1:NGL , 2 );

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






c-----Calculation of rhoAVG0
      rhoAVG0 = 0.D0;
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

                  phi1 = phi_z(ih,ii) * phi_r(il,ii);
                  phi2 = phi_z(ih,jj) * phi_r(il,jj);

                  x = DREAL(rho0(i,j,ib,it)) * phi1*phi2/(2.D0*pi);

                  rhoAVG0(ih,il,it) = rhoAVG0(ih,il,it) + x;


                enddo
              enddo
            enddo


            !call assert(rhoAVG0(ih,il,it).ge.0.D0,'rhoAVG0 negative');

          enddo
        enddo
      enddo






c-----Calculation of drhoAVG
      drhoAVG = COMPLEX( 0.D0 , 0.D0 );
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


                    phi1 = phi_z(ih,i) * phi_r(il,i);
                    phi2 = phi_z(ih,j) * phi_r(il,j);

                    z = 0.5D0 * ( drho_1(i,j,it) + drho_2(i,j,it) );
                    z = z * phi1 * phi2 / (2.D0*pi);

                    drhoAVG(ih,il,it) = drhoAVG(ih,il,it) + z;


                  enddo
                enddo
              enddo
            enddo


          enddo
        enddo
      enddo






c-----Calculation of tauAVG0
      tauAVG0 = 0.D0;
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

                  if( fg1 .eq. fg2 ) CYCLE;

                  ii = i-1+ia_spx(ib);
                  jj = j-1+ia_spx(ib);


                  r      = rb_fam(il);

                   phiz1 =  phi_z(ih,ii)
                  dphiz1 = dphi_z(ih,ii);
                   phir1 =  phi_r(il,ii);
                  dphir1 = dphi_r(il,ii);

                   phiz2 =  phi_z(ih,jj)
                  dphiz2 = dphi_z(ih,jj);
                   phir2 =  phi_r(il,jj);
                  dphir2 = dphi_r(il,jj);


                  if( abs(ml1-ml2) .eq. 0 ) then
                    x = DIMAG(rho0(i,j,ib,it));
                    x = x * (-0.5D0) * ( phiz1*dphiz2 - phiz2*dphiz1 );
                    x = x * phir1*phir2;
                    x = x / (2.D0*pi);

                    tauAVG0(ih,il,it) = tauAVG0(ih,il,it) + x;
                  endif

                  if( abs(ml1+ml2+1) .eq. 0 ) then
                    call assert( r.gt.1.D-6 , 'r close to zero' );

                    x = DBLE(ml2-ml1)/r * phir1*phir2;
                    x = x + dphir1*phir2 - phir1*dphir2;
                    x = x * (-0.5D0) * phiz1*phiz2;
                    x = x / (2.D0*pi);
                    x = x * DREAL(rho0(i,j,ib,it));
                    if( fg1.eq.'g' .and. fg2.eq.'f' ) then
                      x = - x;
                    endif

                    tauAVG0(ih,il,it) = tauAVG0(ih,il,it) + x;
                  endif


                enddo
              enddo
            enddo


            !call assert(tauAVG0(ih,il,it).ge.0.D0,'tauAVG0 negative');

          enddo
        enddo
      enddo






c-----Calculation of dtauAVG
      dtauAVG = COMPLEX( 0.D0 , 0.D0 );
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

                    if( fg1 .eq. fg2 ) CYCLE;


                    r      = rb_fam(il);

                     phiz1 =  phi_z(ih,i)
                    dphiz1 = dphi_z(ih,i);
                     phir1 =  phi_r(il,i);
                    dphir1 = dphi_r(il,i);

                     phiz2 =  phi_z(ih,j)
                    dphiz2 = dphi_z(ih,j);
                     phir2 =  phi_r(il,j);
                    dphir2 = dphi_r(il,j);

                    if( abs(ml1-ml2) .eq. K_multipole ) then
                      z = 0.5D0 * (drho_1(i,j,it)-drho_2(i,j,it));
                      z = z * 0.5D0 * ( phiz1*dphiz2 - phiz2*dphiz1 );
                      z = z * COMPLEX( 0.D0 , 1.D0 );
                      z = z * phir1*phir2;
                      z = z / (2.D0*pi);

                      dtauAVG(ih,il,it) = dtauAVG(ih,il,it) + z;
                    endif

                    if( abs(ml1+ml2+1) .eq. K_multipole ) then
                      call assert( r.gt.1.D-6 , 'r close to zero' );

                      z = DBLE(ml2-ml1)/r * phir1*phir2;
                      z = z + dphir1*phir2 - phir1*dphir2;
                      z = z * (-0.5D0) * phiz1*phiz2;
                      z = z / (2.D0*pi);
                      z = z * 0.5D0 * (drho_1(i,j,it)+drho_2(i,j,it));
                      if( fg1.eq.'g' .and. fg2.eq.'f' ) then
                          z = - z;
                      endif

                       dtauAVG(ih,il,it) = dtauAVG(ih,il,it) + z;
                    endif


                  enddo
                enddo
              enddo
            enddo


          enddo
        enddo
      enddo






c-----Calculation of dr_rhoAVG0 and dz_rhoAVG0
      dr_rhoAVG0 = 0.D0;
      dz_rhoAVG0 = 0.D0;
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


                   phiz1 =  phi_z(ih,ii)
                  dphiz1 = dphi_z(ih,ii);
                   phir1 =  phi_r(il,ii);
                  dphir1 = dphi_r(il,ii);

                   phiz2 =  phi_z(ih,jj)
                  dphiz2 = dphi_z(ih,jj);
                   phir2 =  phi_r(il,jj);
                  dphir2 = dphi_r(il,jj);

                  x = DREAL(rho0(i,j,ib,it));
                  x = x * ( dphir1*phir2 + phir1*dphir2 );
                  x = x * phiz1*phiz2;
                  x = x / (2.D0*pi);
                  dr_rhoAVG0(ih,il,it) = dr_rhoAVG0(ih,il,it) + x;

                  x = DREAL(rho0(i,j,ib,it));
                  x = x * ( dphiz1*phiz2 + phiz1*dphiz2 );
                  x = x * phir1*phir2;
                  x = x / (2.D0*pi);
                  dz_rhoAVG0(ih,il,it) = dz_rhoAVG0(ih,il,it) + x;


                enddo
              enddo
            enddo


          enddo
        enddo
      enddo






c-----Calculation of dr_drhoAVG and dz_drhoAVG
      dr_drhoAVG = COMPLEX( 0.D0 , 0.D0 );
      dz_drhoAVG = COMPLEX( 0.D0 , 0.D0 );
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


                     phiz1 =  phi_z(ih,i)
                    dphiz1 = dphi_z(ih,i);
                     phir1 =  phi_r(il,i);
                    dphir1 = dphi_r(il,i);

                     phiz2 =  phi_z(ih,j)
                    dphiz2 = dphi_z(ih,j);
                     phir2 =  phi_r(il,j);
                    dphir2 = dphi_r(il,j);

                    z = 0.5D0 * ( drho_1(i,j,it) + drho_2(i,j,it) );
                    z = z * ( dphir1*phir2 + phir1*dphir2 );
                    z = z * phiz1*phiz2;
                    z = z / (2.D0*pi);
                    dr_drhoAVG(ih,il,it) = dr_drhoAVG(ih,il,it) + z;

                    z = 0.5D0 * ( drho_1(i,j,it) + drho_2(i,j,it) );
                    z = z * ( dphiz1*phiz2 + phiz1*dphiz2 );
                    z = z * phir1*phir2;
                    z = z / (2.D0*pi);
                    dz_drhoAVG(ih,il,it) = dz_drhoAVG(ih,il,it) + z;


                  enddo
                enddo
              enddo
            enddo


          enddo
        enddo
      enddo






c-----Calculation of grad2_rhoAVG0 and dot_rhoAVG0_drhoAVG
      do it = 1 , 2
        do il = 1 , NGL
          do ih = -NGH , +NGH
            if( ih .eq. 0 ) CYCLE;

            x1 = dr_rhoAVG0(ih,il,it);
            x2 = dz_rhoAVG0(ih,il,it)

            grad2_rhoAVG0(ih,il,it) = x1**2.D0 + x2**2.D0;

            dot_rhoAVG0_drhoAVG(ih,il,it) = + x1 * dr_drhoAVG(ih,il,it)
     &                                      + x2 * dz_drhoAVG(ih,il,it);

          enddo
        enddo
      enddo






c-----Calculation of f0 and df
      do it = 1 , 2
        do il = 1 , NGL
          do ih = -NGH , +NGH
            if( ih .eq. 0 ) CYCLE;

            fac    = (3.D0/5.D0) * (6.D0*pi*pi)**(2.D0/3.D0);
            tau0TF = fac * rhoAVG0(ih,il,it)**(5.D0/3.D0);

            x = rhoAVG0(ih,il,it) * tauAVG0(ih,il,it);
            x = x - 0.25D0 * grad2_rhoAVG0(ih,il,it);
            x = x / ( rhoAVG0(ih,il,it) * tau0TF );

            f0(ih,il,it) = x;


            z =   + drhoAVG(ih,il,it)*tauAVG0(ih,il,it);
            z = z + dtauAVG(ih,il,it)*rhoAVG0(ih,il,it);
            z = z * 2.D0;
            z = z - dot_rhoAVG0_drhoAVG(ih,il,it);
            z = z / ( 2.D0 * rhoAVG0(ih,il,it) * tau0TF );
            z = z - (8.D0/3.D0)*x*drhoAVG(ih,il,it)/rhoAVG0(ih,il,it);

            df(ih,il,it) = z;


          enddo
        enddo
      enddo






c-----Calculation of C0 and dC
      do it = 1 , 2
        do il = 1 , NGL
          do ih = -NGH , +NGH
            if( ih .eq. 0 ) CYCLE;

            C0(ih,il,it) = 1.D0 / ( 1.D0 + f0(ih,il,it)**2.D0 );


            dC(ih,il,it) = - 2.D0 * C0(ih,il,it)**2.D0
     &                            * f0(ih,il,it)
     &                            * df(ih,il,it);


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
