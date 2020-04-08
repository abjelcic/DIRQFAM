c======================================================================c

      REAL*8 function fam_strength( lpr , it )

c======================================================================c

      USE dirqfampar;
      USE simplex;
      USE fmatrix;
      USE drho;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;
      INTEGER it;



      DOUBLE COMPLEX trace;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_strength() ************************';
      write(6,*) '';
      endif






c-----Calculates S(f,omega) = -1/pi * Im[Tr[hermconj(f)*drho(omega)]]
      trace = COMPLEX( 0.D0 , 0.D0 );

      do j = 1 , N_total
          do i = 1 , N_total
              trace = trace + DCONJG(f1_JK(i,j,it))*drho_1(i,j,it)
     &                      + DCONJG(f2_JK(i,j,it))*drho_2(i,j,it);
          enddo
      enddo

      pi = 3.14159265358979324D0;
      fam_strength = - 1.D0/pi * DIMAG( trace );






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_strength() **************************';
      write(6,*) '';
      endif

      return;
      end;
