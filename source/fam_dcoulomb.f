c======================================================================c

      subroutine fam_dcoulomb( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE fam_green;
      USE dlaplace;
      USE dcoul;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      DOUBLE COMPLEX acc0, acc1, acc2, acc3;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_dcoulomb() ************************';
      write(6,*) '';
      endif


      dVCou_0 = COMPLEX( 0.D0 , 0.D0 );
      dVCou_r = COMPLEX( 0.D0 , 0.D0 );
      dVCou_p = COMPLEX( 0.D0 , 0.D0 );
      dVCou_z = COMPLEX( 0.D0 , 0.D0 );


      if( i_coulomb .eq. 0 ) then
          return;
      endif


      do il = 1 , NGL
          do ih = -NGH , +NGH
              if( ih .eq. 0 ) CYCLE;

              acc0 = COMPLEX( 0.D0 , 0.D0 );
              acc1 = COMPLEX( 0.D0 , 0.D0 );
              acc2 = COMPLEX( 0.D0 , 0.D0 );
              acc3 = COMPLEX( 0.D0 , 0.D0 );
              do il1 = 1 , NGL
                  do ih1 = -NGH , +NGH
                      if( ih1 .eq. 0 ) CYCLE;
                      acc0 = acc0 + G3(ih1,il1,ih,il)*ldrho_vp(ih1,il1);
                      acc1 = acc1 + G1(ih1,il1,ih,il)*  ldj_1p(ih1,il1);
                      acc2 = acc2 + G2(ih1,il1,ih,il)*  ldj_2p(ih1,il1);
                      acc3 = acc3 + G3(ih1,il1,ih,il)*  ldj_3p(ih1,il1);
                  enddo
              enddo

              dVCou_0(ih,il) = acc0;
              dVCou_r(ih,il) = ( + acc1 + acc2 ) / DSQRT(2.D0);
              dVCou_p(ih,il) = ( - acc1 + acc2 ) / DSQRT(2.D0);
              dVCou_z(ih,il) = acc3;

          enddo
      enddo



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_dcoulomb() **************************';
      write(6,*) '';
      endif

      return;
      end;
