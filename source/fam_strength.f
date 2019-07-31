c======================================================================c

      REAL*8 function fam_strength( lpr , it )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /mathco/ zero, one, two, half, third, pi;

      CHARACTER fg_spx;
      common /simplex/ N_total         , N_blocks        ,
     &                 ia_spx(NBX)     , id_spx(NBX)     ,
     &                 nf_size(NBX)    , ng_size(NBX)    ,
     &                 nz_spx(NBSX,NBX), nr_spx(NBSX,NBX),
     &                 ml_spx(NBSX,NBX), fg_spx(NBSX,NBX);

      COMPLEX*16 f1_JK, f2_JK;
      common /f_matrix/ f1_JK( NTX , NTX , 2 ),
     &                  f2_JK( NTX , NTX , 2 );

      COMPLEX*16 drho_1, drho_2;
      common /delta_rho/ drho_1( NTX , NTX , 2 ),
     &                   drho_2( NTX , NTX , 2 );



      COMPLEX*16 trace;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN fam_strength() ************************';
      write(6,*) '';
      endif






c-----Calculates S(f,omega) = -1/pi * Im[Tr[hermconj(f)*drho(omega)]]
      trace = COMPLEX( 0.D0 , 0.D0 );

      do i = 1 , N_total
          do j = 1 , N_total
              trace = trace + DCONJG(f1_JK(i,j,it))*drho_1(i,j,it)
     &                      + DCONJG(f2_JK(i,j,it))*drho_2(i,j,it);
          enddo
      enddo


      fam_strength = - 1.D0/pi * DIMAG( trace );






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END fam_strength() **************************';
      write(6,*) '';
      endif

      return;
      end;
