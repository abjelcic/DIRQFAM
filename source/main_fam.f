c======================================================================c

      subroutine main_fam( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z);
      include 'dirqfam.par'
      LOGICAL lpr;

      CHARACTER parname*10;
      common /partyp/ parname;
      common /basnnn/ n0f, n0b;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN main_fam() *****************************';
      write(6,*) '';
      endif


      call assert( MOD(n0f,2).eq.0 , 'n0f has to be even number' );
      call assert( NGH.ge.30  , 'NGH has to be at least 30' );
      call assert( NGL.ge.30  , 'NGL has to be at least 30' );
      call assert( NGH.lt.100 , 'NGH has to be less than 100' );
      call assert( NGL.lt.100 , 'NGL has to be less than 100' );
      call assert( parname.eq.'DD-PC1' , 'Only DD-PC1 supports QFAM' );


      write(6,'(a)') '                                               ';
      write(6,'(a)') '                                               ';
      write(6,'(a)') '                                               ';
      write(6,'(a)') '  *******************************************  ';
      write(6,'(a)') '  *      Relativistic Quasiparticle RPA     *  ';
      write(6,'(a)') '  *           Simplex-y Basis               *  ';
      write(6,'(a)') '  *                                         *  ';
      write(6,'(a)') '  *   Finite Amplitude Method Calculation   *  ';
      write(6,'(a)') '  *         with Density-Dependent          *  ';
      write(6,'(a)') '  *          Point-Coupling Force           *  ';
      write(6,'(a)') '  *         and Separable Pairing           *  ';
      write(6,'(a)') '  * --------------------------------------- *  ';
      write(6,'(a)') '  *                                         *  ';
      write(6,'(a)') '  *           A.Bjelcic, T.Niksic           *  ';
      write(6,'(a)') '  *                                         *  ';
      write(6,'(a)') '  *******************************************  ';
      write(6,'(a)') '                                               ';
      call flush(6);






c-----Constructs simplex-y basis quantum numbers
      call base_simplex( .false. );

c-----Constructs U and V matrices in simplex-y basis
      do it = 1 , 2
          call construct_v( .false. , it );
          call construct_u( .false. , it );
      enddo

c-----Consistency check of transformation to simplex-y basis
      if( .true. ) then
          call check_gs_dens( .false. );
          do it = 1 , 2
              call check_unitarity( .false. , it );
          enddo
      endif

c-----Initializes QFAM submodule
      call init_fam( .false. );

c-----Executes QFAM submodule
      call start_fam( .false. );






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END main_fam() *******************************';
      write(6,*) '';
      endif

      return;
      end;
