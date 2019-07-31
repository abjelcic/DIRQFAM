c======================================================================c

      subroutine init_basis( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /baspar/ hom, hb0, b0;
      common /defbas/ beta0, q, bp, bz;
      common /gaussh/ xh(0:NGH), wh(0:NGH), zb(0:NGH);
      common /gaussl/ xl(0:NGL), wl(0:NGL), sxl(0:NGL), rb(0:NGL);

      CHARACTER fg_spx;
      common /simplex/ N_total         , N_blocks        ,
     &                 ia_spx(NBX)     , id_spx(NBX)     ,
     &                 nf_size(NBX)    , ng_size(NBX)    ,
     &                 nz_spx(NBSX,NBX), nr_spx(NBSX,NBX),
     &                 ml_spx(NBSX,NBX), fg_spx(NBSX,NBX);

      common /wavefunc/ qh ( -NGH:NGH , NTX ),
     &                  qh1( -NGH:NGH , NTX ),
     &                  ql (    1:NGL , NTX ),
     &                  ql1(    1:NGL , NTX );

      common /basis/ qhql ( -NGH:NGH , 1:NGL , NTX ),
     &               qh1ql( -NGH:NGH , 1:NGL , NTX ),
     &               qhql1( -NGH:NGH , 1:NGL , NTX ),
     &               wqhql(    1:NGH , 1:NGL , NTX );

      common /k_PHI/ k_PHI, k_PHIr, k_dzPHI, k_drPHI;
      common /PHI/ PHI_U    (    NTX , KTRUNC ),
     &             PHI_SVt  ( KTRUNC , NCOORD ),
     &
     &             PHIr_U   (    NTX , KTRUNC ),
     &             PHIr_SVt ( KTRUNC , NCOORD ),
     &
     &             dzPHI_U  (    NTX , KTRUNC ),
     &             dzPHI_SVt( KTRUNC , NCOORD ),
     &
     &             drPHI_U  (    NTX , KTRUNC ),
     &             drPHI_SVt( KTRUNC , NCOORD );



      REAL*8 PHI  ( NTX , NCOORD );
      REAL*8 PHIr ( NTX , NCOORD );
      REAL*8 dzPHI( NTX , NCOORD );
      REAL*8 drPHI( NTX , NCOORD );

      parameter( THRESHOLD = 1.D-12 );
      REAL*8 S (                   min(NTX,NCOORD) );
      REAL*8 U (             NTX , min(NTX,NCOORD) );
      REAL*8 VT( min(NTX,NCOORD) ,         NCOORD  );

      parameter( NWORK = 10000*NTX );
      REAL*8     WORK( NWORK );
      INTEGER*4 IWORK( 12*min(NTX,NCOORD) );



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN init_basis() **************************';
      write(6,*) '';
      endif






c-----Calculation of qh
      do ib = 1 , N_blocks
          do i = 1 , id_spx(ib)
              nz = nz_spx(i,ib);
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  z = DBLE(isign(1,ih)) * zb(abs(ih));
                  qh(ih,i-1+ia_spx(ib)) = phi_nz(nz,z);
              enddo
          enddo
      enddo






c-----Calculation of ql
      do ib = 1 , N_blocks
          do i = 1 , id_spx(ib)
              nr = nr_spx(i,ib);
              ml = ml_spx(i,ib);
              do il = 1 , NGL
                  r = rb(il);
                  ql(il,i-1+ia_spx(ib)) = phi_nr_ml(nr,abs(ml),r);
              enddo
          enddo
      enddo






c-----Calculation of qh1
      do ib = 1 , N_blocks
          do i = 1 , id_spx(ib)
              nz = nz_spx(i,ib);
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  z = DBLE(isign(1,ih)) * zb(abs(ih));
                  qh1(ih,i-1+ia_spx(ib)) = d_phi_nz(nz,z);
              enddo
          enddo
      enddo






c-----Calculation of ql1
      do ib = 1 , N_blocks
          do i = 1 , id_spx(ib)
              nr = nr_spx(i,ib);
              ml = ml_spx(i,ib);
              do il = 1 , NGL
                  r = rb(il);
                  ql1(il,i-1+ia_spx(ib)) = d_phi_nr_ml(nr,abs(ml),r);
              enddo
          enddo
      enddo






c-----Calculation of qhql, qh1ql and qhql1
      do i = 1 , N_total
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  qhql (ih,il,i) = qh (ih,i) * ql (il,i);
                  qh1ql(ih,il,i) = qh1(ih,i) * ql (il,i);
                  qhql1(ih,il,i) = qh (ih,i) * ql1(il,i);
              enddo
          enddo
      enddo






c-----Calculation of wqhql
      do i = 1 , N_total
          do il = 1 , NGL
              do ih = 1 , NGH
                  ! Notice that in Ground State code, Gauss-Hermite
                  ! weights wh() are multiplied by a factor 2
                  !                                  |
                  !                                  ˇ
                  w = (0.5D0*b0*b0*b0*bp*bp*bz)*(0.5D0*wh(ih))*wl(il);
                  wqhql(ih,il,i) = DSQRT(w) * qhql(ih,il,i);
              enddo
          enddo
      enddo






c-----Calculation of PHI_U, PHI_SVt
      do i = 1 , N_total
          j = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  j = j+1;
                  PHI(i,j) = qh(ih,i)*ql(il,i);
              enddo
          enddo
      enddo

      call dgesvdx(       'V' ,             'V' ,
     &                    'V' ,
     &                N_total ,          NCOORD ,
     &                    PHI ,             NTX ,
     &              THRESHOLD ,         1.D+100 ,
     &                      0 ,               0 ,
     &                  k_PHI ,               S ,
     &                      U ,             NTX ,
     &                      VT, min(NTX,NCOORD) ,
     &                   WORK ,              -1 ,
     &                  IWORK ,            INFO   );
      if( WORK(1) .gt. DBLE(NWORK) ) then
          write(6,'(a,f10.1)') 'WORK(1) = ' , WORK(1);
          write(6,'(a,i10)'  ) 'NWORK   = ' ,   NWORK;
          stop 'Error: NWORK too small, please increase it!';
      endif
      call dgesvdx(       'V' ,             'V' ,
     &                    'V' ,
     &                N_total ,          NCOORD ,
     &                    PHI ,             NTX ,
     &              THRESHOLD ,         1.D+100 ,
     &                      0 ,               0 ,
     &                  k_PHI ,               S ,
     &                      U ,             NTX ,
     &                      VT, min(NTX,NCOORD) ,
     &                   WORK ,           NWORK ,
     &                  IWORK ,            INFO   );
      if( INFO .ne. 0 ) then
          stop 'Error: dgesvdx()!';
      endif
      if( k_PHI .gt. KTRUNC ) then
          write(6,'(a,i4)') 'k_PHI  = ' ,  k_PHI;
          write(6,'(a,i4)') 'KTRUNC = ' , KTRUNC;
          stop 'Error: KTRUNC too small, please incrase it';
      endif

      do j = 1 , k_PHI
          do i = 1 , N_total
              PHI_U(i,j) = U(i,j);
          enddo
      enddo
      do j = 1 , NCOORD
          do i = 1 , k_PHI
              PHI_SVt(i,j) = S(i) * VT(i,j);
          enddo
      enddo

      if( lpr ) then
          write(6,'(a,i4)') 'k_PHI = ' , k_PHI;

          do i = 1 , N_total
              j = 0;
              do il = 1 , NGL
                  do ih = -NGH , NGH
                      if( ih .eq. 0 ) CYCLE;
                      j = j+1;
                      PHI(i,j) = qh(ih,i)*ql(il,i);
                  enddo
              enddo
          enddo

          acc1 = 0.D0;
          acc2 = 0.D0;
          do i = 1 , N_total
              do j = 1 , NCOORD

                  acc1 = acc1 + PHI(i,j)*PHI(i,j);

                  x = 0.D0;
                  do k = 1 , k_PHI
                      x = x + PHI_U(i,k)*PHI_SVt(k,j);
                  enddo
                  acc2 = acc2 + (x-PHI(i,j))*(x-PHI(i,j));

              enddo
          enddo

          write(6,'(a,E20.10)')
     &    '|| PHI - PHI(k_PHI) ||_F / || PHI ||_F = ',
     &    DSQRT(acc2) / DSQRT(acc1);

      endif






c-----Calculation of PHIr_U, PHIr_SVt
      do i = 1 , N_total
          j = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  j = j+1;
                  PHIr(i,j) = qh(ih,i)*ql(il,i) / rb(il);
              enddo
          enddo
      enddo

      call dgesvdx(       'V' ,             'V' ,
     &                    'V' ,
     &                N_total ,          NCOORD ,
     &                   PHIr ,             NTX ,
     &              THRESHOLD ,         1.D+100 ,
     &                      0 ,               0 ,
     &                 k_PHIr ,               S ,
     &                      U ,             NTX ,
     &                      VT, min(NTX,NCOORD) ,
     &                   WORK ,              -1 ,
     &                  IWORK ,            INFO   );
      if( WORK(1) .gt. DBLE(NWORK) ) then
          write(6,'(a,f10.1)') 'WORK(1) = ' , WORK(1);
          write(6,'(a,i10)'  ) 'NWORK   = ' ,   NWORK;
          stop 'Error: NWORK too small, please increase it!';
      endif
      call dgesvdx(       'V' ,             'V' ,
     &                    'V' ,
     &                N_total ,          NCOORD ,
     &                   PHIr ,             NTX ,
     &              THRESHOLD ,         1.D+100 ,
     &                      0 ,               0 ,
     &                 k_PHIr ,               S ,
     &                      U ,             NTX ,
     &                      VT, min(NTX,NCOORD) ,
     &                   WORK ,           NWORK ,
     &                  IWORK ,            INFO   );
      if( INFO .ne. 0 ) then
          stop 'Error: dgesvdx()!';
      endif
      if( k_PHIr .gt. KTRUNC ) then
          write(6,'(a,i4)') 'k_PHIr = ' , k_PHIr;
          write(6,'(a,i4)') 'KTRUNC = ' , KTRUNC;
          stop 'Error: KTRUNC too small, please incrase it';
      endif

      do j = 1 , k_PHIr
          do i = 1 , N_total
              PHIr_U(i,j) = U(i,j);
          enddo
      enddo
      do j = 1 , NCOORD
          do i = 1 , k_PHIr
              PHIr_SVt(i,j) = S(i) * VT(i,j);
          enddo
      enddo

      if( lpr ) then
          write(6,'(a,i4)') 'k_PHIr = ' , k_PHIr;

          do i = 1 , N_total
              j = 0;
              do il = 1 , NGL
                  do ih = -NGH , NGH
                      if( ih .eq. 0 ) CYCLE;
                      j = j+1;
                      PHIr(i,j) = qh(ih,i)*ql(il,i) / rb(il);
                  enddo
              enddo
          enddo

          acc1 = 0.D0;
          acc2 = 0.D0;
          do i = 1 , N_total
              do j = 1 , NCOORD

                  acc1 = acc1 + PHIr(i,j)*PHIr(i,j);

                  x = 0.D0;
                  do k = 1 , k_PHIr
                      x = x + PHIr_U(i,k)*PHIr_SVt(k,j);
                  enddo
                  acc2 = acc2 + (x-PHIr(i,j))*(x-PHIr(i,j));

              enddo
          enddo

          write(6,'(a,E20.10)')
     &    '|| PHIr - PHIr(k_PHIr) ||_F / || PHIr ||_F = ',
     &    DSQRT(acc2) / DSQRT(acc1);

      endif






c-----Calculation of dzPHI_U, dzPHI_SVt
      do i = 1 , N_total
          j = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  j = j+1;
                  dzPHI(i,j) = qh1(ih,i)*ql(il,i);
              enddo
          enddo
      enddo

      call dgesvdx(       'V' ,             'V' ,
     &                    'V' ,
     &                N_total ,          NCOORD ,
     &                  dzPHI ,             NTX ,
     &              THRESHOLD ,         1.D+100 ,
     &                      0 ,               0 ,
     &                k_dzPHI ,               S ,
     &                      U ,             NTX ,
     &                      VT, min(NTX,NCOORD) ,
     &                   WORK ,              -1 ,
     &                  IWORK ,            INFO   );
      if( WORK(1) .gt. DBLE(NWORK) ) then
          write(6,'(a,f10.1)') 'WORK(1) = ' , WORK(1);
          write(6,'(a,i10)'  ) 'NWORK   = ' ,   NWORK;
          stop 'Error: NWORK too small, please increase it!';
      endif
      call dgesvdx(       'V' ,             'V' ,
     &                    'V' ,
     &                N_total ,          NCOORD ,
     &                  dzPHI ,             NTX ,
     &              THRESHOLD ,         1.D+100 ,
     &                      0 ,               0 ,
     &                k_dzPHI ,               S ,
     &                      U ,             NTX ,
     &                      VT, min(NTX,NCOORD) ,
     &                   WORK ,           NWORK ,
     &                  IWORK ,            INFO   );
      if( INFO .ne. 0 ) then
          stop 'Error: dgesvdx()!';
      endif
      if( k_dzPHI .gt. KTRUNC ) then
          write(6,'(a,i4)') 'k_dzPHI = ' , k_dzPHI;
          write(6,'(a,i4)') 'KTRUNC  = ' ,  KTRUNC;
          stop 'Error: KTRUNC too small, please incrase it';
      endif

      do j = 1 , k_dzPHI
          do i = 1 , N_total
              dzPHI_U(i,j) = U(i,j);
          enddo
      enddo
      do j = 1 , NCOORD
          do i = 1 , k_dzPHI
              dzPHI_SVt(i,j) = S(i) * VT(i,j);
          enddo
      enddo

      if( lpr ) then
          write(6,'(a,i4)') 'k_dzPHI = ' , k_dzPHI;

          do i = 1 , N_total
              j = 0;
              do il = 1 , NGL
                  do ih = -NGH , NGH
                      if( ih .eq. 0 ) CYCLE;
                      j = j+1;
                      dzPHI(i,j) = qh1(ih,i)*ql(il,i);
                  enddo
              enddo
          enddo

          acc1 = 0.D0;
          acc2 = 0.D0;
          do i = 1 , N_total
              do j = 1 , NCOORD

                  acc1 = acc1 + dzPHI(i,j)*dzPHI(i,j);

                  x = 0.D0;
                  do k = 1 , k_dzPHI
                      x = x + dzPHI_U(i,k)*dzPHI_SVt(k,j);
                  enddo
                  acc2 = acc2 + (x-dzPHI(i,j))*(x-dzPHI(i,j));

              enddo
          enddo

          write(6,'(a,E20.10)')
     &    '|| dzPHI - dzPHI(k_dzPHI) ||_F / || dzPHI ||_F = ',
     &    DSQRT(acc2) / DSQRT(acc1);

      endif






c-----Calculation of drPHI_U, drPHI_SVt
      do i = 1 , N_total
          j = 0;
          do il = 1 , NGL
              do ih = -NGH , NGH
                  if( ih .eq. 0 ) CYCLE;
                  j = j+1;
                  drPHI(i,j) = qh(ih,i)*ql1(il,i);
              enddo
          enddo
      enddo

      call dgesvdx(       'V' ,             'V' ,
     &                    'V' ,
     &                N_total ,          NCOORD ,
     &                  drPHI ,             NTX ,
     &              THRESHOLD ,         1.D+100 ,
     &                      0 ,               0 ,
     &                k_drPHI ,               S ,
     &                      U ,             NTX ,
     &                      VT, min(NTX,NCOORD) ,
     &                   WORK ,              -1 ,
     &                  IWORK ,            INFO   );
      if( WORK(1) .gt. DBLE(NWORK) ) then
          write(6,'(a,f10.1)') 'WORK(1) = ' , WORK(1);
          write(6,'(a,i10)'  ) 'NWORK   = ' ,   NWORK;
          stop 'Error: NWORK too small, please increase it!';
      endif
      call dgesvdx(       'V' ,             'V' ,
     &                    'V' ,
     &                N_total ,          NCOORD ,
     &                  drPHI ,             NTX ,
     &              THRESHOLD ,         1.D+100 ,
     &                      0 ,               0 ,
     &                k_drPHI ,               S ,
     &                      U ,             NTX ,
     &                      VT, min(NTX,NCOORD) ,
     &                   WORK ,           NWORK ,
     &                  IWORK ,            INFO   );
      if( INFO .ne. 0 ) then
          stop 'Error: dgesvdx()!';
      endif
      if( k_drPHI .gt. KTRUNC ) then
          write(6,'(a,i4)') 'k_drPHI = ' , k_drPHI;
          write(6,'(a,i4)') 'KTRUNC  = ' ,  KTRUNC;
          stop 'Error: KTRUNC too small, please incrase it';
      endif

      do j = 1 , k_drPHI
          do i = 1 , N_total
              drPHI_U(i,j) = U(i,j);
          enddo
      enddo
      do j = 1 , NCOORD
          do i = 1 , k_drPHI
              drPHI_SVt(i,j) = S(i) * VT(i,j);
          enddo
      enddo

      if( lpr ) then
          write(6,'(a,i4)') 'k_drPHI = ' , k_drPHI;

          do i = 1 , N_total
              j = 0;
              do il = 1 , NGL
                  do ih = -NGH , NGH
                      if( ih .eq. 0 ) CYCLE;
                      j = j+1;
                      drPHI(i,j) = qh(ih,i)*ql1(il,i);
                  enddo
              enddo
          enddo

          acc1 = 0.D0;
          acc2 = 0.D0;
          do i = 1 , N_total
              do j = 1 , NCOORD

                  acc1 = acc1 + drPHI(i,j)*drPHI(i,j);

                  x = 0.D0;
                  do k = 1 , k_drPHI
                      x = x + drPHI_U(i,k)*drPHI_SVt(k,j);
                  enddo
                  acc2 = acc2 + (x-drPHI(i,j))*(x-drPHI(i,j));

              enddo
          enddo

          write(6,'(a,E20.10)')
     &    '|| drPHI - drPHI(k_drPHI) ||_F / || drPHI ||_F = ',
     &    DSQRT(acc2) / DSQRT(acc1);

      endif






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END init_basis() ****************************';
      write(6,*) '';
      endif

      return;
      end;
