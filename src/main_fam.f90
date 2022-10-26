!======================================================================!

      subroutine main_fam()

!======================================================================!
      implicit none;


      write(6,'(a)') '                                                ';
      write(6,'(a)') '                                                ';
      write(6,'(a)') '                                                ';
      write(6,'(a)') '  ********************************************  ';
      write(6,'(a)') '  *              DIRQFAM v2.0.0              *  ';
      write(6,'(a)') '  *                                          *  ';
      write(6,'(a)') '  *      Relativistic Quasiparticle RPA      *  ';
      write(6,'(a)') '  *             Simplex-y Basis              *  ';
      write(6,'(a)') '  *                                          *  ';
      write(6,'(a)') '  *   Finite Amplitude Method Calculation    *  ';
      write(6,'(a)') '  *          with Density-Dependent          *  ';
      write(6,'(a)') '  *  Meson-Exchange or Point-Coupling Force  *  ';
      write(6,'(a)') '  *          and Separable Pairing           *  ';
      write(6,'(a)') '  * ---------------------------------------- *  ';
      write(6,'(a)') '  *                                          *  ';
      write(6,'(a)') '  *            A.Bjelcic, T.Niksic           *  ';
      write(6,'(a)') '  *                                          *  ';
      write(6,'(a)') '  ********************************************  ';
      write(6,'(a)') '                                                ';
      call flush(6);


          ! Reading QFAM input data.
          call read_input();

          ! Constructing the simplex-y basis quantum numbers.
          call base_simplex();

          ! Allocating the memory used in QFAM submodule.
          call allocqfam();

          ! Constructing the U and V matrices in simplex-y basis.
          call construct_u();
          call construct_v();
          call construct_qpenergies();

          ! Consistency check of transformation to simplex-y basis.
          if( .true. ) then
              call check_gs_dens();
              call check_unitarity();
          end if

          ! Initializing the QFAM submodule.
          call init_fam();

          ! Executing the QFAM submodule.
          call start_fam();

          ! Deallocating the memory used in QFAM submodule.
          call deallocqfam();


      return;
      end subroutine main_fam
