c======================================================================c

      subroutine fam_spurious( lpr )

c======================================================================c

      USE dirqfampar;
      USE fam;
      USE simplex;
      USE xyfam;
      USE spurious;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_spurious() ************************';
      write(6,*) '';
      endif


      if( J_multipole.ne.1 .and. J_multipole.ne.3 ) then
          return;
      endif
      if( K_multipole.ne.0 .and. K_multipole.ne.1 ) then
          return;
      endif


      lamR = COMPLEX( 0.D0 , 0.D0 );
      lamP = COMPLEX( 0.D0 , 0.D0 );
      do it = 1 , 2
          do j = 1 , N_total
              do i = 1 , N_total
                  lamR = lamR + ( - x_fam(i,j,it) + y_fam(i,j,it) )
     &                        * DCONJG( p20(i,j,it) );
                  lamP = lamP + ( + x_fam(i,j,it) + y_fam(i,j,it) )
     &                        * DCONJG( r20(i,j,it) );
              enddo
          enddo
      enddo
      lamR = lamR / ( RcmPcm_commutator(1)+RcmPcm_commutator(2) );
      lamP = lamP / ( RcmPcm_commutator(1)+RcmPcm_commutator(2) );


      do it = 1 , 2
          do j = 1 , N_total
              do i = 1 , N_total
                  x_fam(i,j,it) = x_fam(i,j,it) - lamR*r20(i,j,it)
     &                                          - lamP*p20(i,j,it);
                  y_fam(i,j,it) = y_fam(i,j,it) + lamR*r20(i,j,it)
     &                                          - lamP*p20(i,j,it);
              enddo
          enddo
      enddo






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_spurious() **************************';
      write(6,*) '';
      endif

      return;
      end;
