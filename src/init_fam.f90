!======================================================================!

      subroutine init_fam()

!======================================================================!
      implicit none;

          write(6,'(x,a)') 'Initializing QFAM submodule.';
          call flush(6);

          ! Calculation of the Gaussian quadrature weights+nodes and basis functions.
          call init_basis();

          ! Calculation of the ground state densities and meson fields.
          call init_gs();

          ! Calculation of the multipole excitation operator.
          call init_multipole();

          ! Calculation of the Rcm and Pcm operators.
          call init_spurious();

          ! Preparation of mesons.
          call init_mesons();

          ! Preparation of Green's functions for electromagnetic interaction.
          call init_coulomb();

          ! Preparation of pairing (W coefficients).
          call init_pairing();

          write(6,'(x,a)') 'Initialization complete.';
          call flush(6);

      return;
      end subroutine init_fam
