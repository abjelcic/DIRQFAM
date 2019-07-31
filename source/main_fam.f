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


      if( MOD(n0f,2) .ne. 0 ) then
          stop 'Error: n0f has to be even number!';
      endif
      if( parname .ne. 'DD-PC1' ) then
          stop 'Error: FAM implemented only for DD-PC1!';
      endif


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

c-----Initializes FAM submodule
      call init_fam( .false. );

c-----Executes FAM submodule
      call start_fam( .false. );






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END main_fam() *******************************';
      write(6,*) '';
      endif

      return;
      end;
