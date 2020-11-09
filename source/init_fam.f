c======================================================================c

      subroutine init_fam( lpr )

c======================================================================c

      USE dirqfampar;
      USE simplex;
      USE fam;
      USE ddpc1ddme2;
      USE pairparams;
      USE u_energy;
      USE v_energy;
      USE fam_energies;
      USE nnz_blocks;
      IMPLICIT DOUBLE PRECISION(a-h,o-z)
      IMPLICIT INTEGER(i-n)
      LOGICAL lpr;



      CHARACTER fg1, fg2;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN init_fam() ****************************';
      write(6,*) '';
      endif

      write(6,'(a)') 'Initializing QFAM submodule.';
      call flush(6);






c-----Reading QFAM parameters
      infam = 123;
      open( infam , file = 'dirqfam.dat' , status = 'old' );
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,'(18x,i9)')    i_calculation_type;
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,'(18x,i9)')    i_coulomb;
      read(infam,'(18x,i9)')    i_pairing;
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,'(18x,i9)')    J_multipole;
      read(infam,'(18x,i9)')    K_multipole;
      read(infam,'(18x,i9)')    ISO;
      read(infam,'(18x,1f9.3)') gamma_smear;
      read(infam,*);
      read(infam,'(18x,1f9.3)') omega_start;
      read(infam,'(18x,1f9.3)') omega_end;
      read(infam,'(18x,1f9.3)') delta_omega;
      read(infam,*);
      read(infam,'(18x,1f9.3)') omega_print;
      close(infam);

      tape_strength   = 100;
      tape_rhov       = 110;
      tape_nuclocfunc = 120;





c-----DD-PC1 parameters
      a_s     = - 10.0462D0; ![fm^2]
      b_s     = -  9.1504D0; ![fm^2]
      c_s     = -  6.4273D0; ![fm^2]
      d_s     = +  1.3724D0; ![    ]

      a_v     = +  5.9195D0; ![fm^2]
      b_v     = +  8.8637D0; ![fm^2]
      c_v     = +  0.0000D0; ![fm^2]
      d_v     = +  0.6584D0; ![    ]

      a_tv    = +  0.0000D0; ![fm^2]
      b_tv    = +  1.8360D0; ![fm^2]
      c_tv    = +  0.0000D0; ![fm^2]
      d_tv    = +  0.6403D0; ![    ]

      del_s   = -  0.8149D0; ![fm^4]

      rho_sat = +  0.1520D0; ![fm^-3]

c-----DD-ME2 parameters
      a_sig   = +   1.388148180042941D0; ![       ]
      b_sig   = +   1.094300000000000D0; ![       ]
      c_sig   = +   1.705700000000000D0; ![       ]
      d_sig   = +   0.442066950716287D0; ![       ]

      a_ome   = +   1.389291065981963D0; ![       ]
      b_ome   = +   0.923974841420762D0; ![       ]
      c_ome   = +   1.462000000000000D0; ![       ]
      d_ome   = +   0.477491545490170D0; ![       ]

      a_rho   = +   0.564700000000000D0; ![       ]

      g0_sig  = +  10.539600000000000D0; ![       ]
      g0_ome  = +  13.018900000000000D0; ![       ]
      g0_rho  = +   3.683600000000000D0; ![       ]

      m_sig   = + 550.123800000000000D0; ![MeV/c^2]
      m_ome   = + 783.000000000000000D0; ![MeV/c^2]
      m_rho   = + 763.000000000000000D0; ![MeV/c^2]

      rho_sat = +   0.152000000000000D0; ![ fm^-3 ]

c-----TMR separable pairing parameters
      G_pairing = + 728.000000000D0; ![MeV/fm^3]
      a_pairing = +   0.644204936D0; ![   fm   ]






c-----Calculation of QFAM energies
      do it = 1 , 2
          ie = 0;
          do ib = 1 , N_blocks
              do i = 1 , id_spx(ib)
                  Eu = E_fam_u(i,ib,it);
                  Ev = E_fam_v(i,ib,it);
                  call assert( DABS(Eu-Ev).lt.1.D-15 , 'Eu =/= Ev' );

                  ie = ie + 1;
                  E_fam(ie,it) = Eu;
              enddo
          enddo
          call assert( ie .eq. N_total , 'ie =/= N_total' );
      enddo






c-----Calculation of nnz_blocks
      dh_nnz     = .false.;
      dDelta_nnz = .false.;
      dkappa_nnz = .false.;
      f_nnz      = .false.;
      do ib1 = 1 , N_blocks
          do i = 1 , id_spx(ib1);
              fg1 = fg_spx(i,ib1);
              ml1 = ml_spx(i,ib1);
              do ib2 = 1 , N_blocks
                  do j = 1 , id_spx(ib2)
                      fg2 = fg_spx(j,ib2);
                      ml2 = ml_spx(j,ib2);

                      if( fg1.eq.'f' .and. fg2.eq.'f' ) then
                          if( abs(ml1-ml2).eq.K_multipole ) then
                              dh_nnz    (ib1,ib2) = .true.;
                              dDelta_nnz(ib1,ib2) = .true.;
                              dkappa_nnz(ib1,ib2) = .true.;
                              f_nnz     (ib1,ib2) = .true.;
                          endif
                      endif
                      if( fg1.eq.'g' .and. fg2.eq.'g' ) then
                          if( abs(ml1-ml2).eq.K_multipole ) then
                              dh_nnz(ib1,ib2) = .true.;
                              f_nnz (ib1,ib2) = .true.;
                          endif
                      endif
                      if( fg1.ne.fg2 ) then
                          if( abs(ml1+ml2+1).eq.K_multipole  .or.
     &                        abs(ml1-ml2)  .eq.K_multipole      ) then
                              dh_nnz(ib1,ib2) = .true.;
                          endif
                      endif

                  enddo
              enddo
          enddo
      enddo






c-----Calculation of the quadrature weights+nodes and basis functions
      call init_basis( .false. );






c-----Calculation of the ground state densities and meson fields
      call init_gs( .false. );






c-----Calculation of the multipole excitation operator
      call init_multipole( .false. );






c-----Calculation of the Rcm and Pcm operators
      call init_spurious( .false. );






c-----Preparation of mesons
      call init_mesons( .false. );






c-----Preparation of Coulomb (Green's function)
      call init_coulomb( .false. );






c-----Preparation of pairing (W coefficients)
      call init_pairing( .false. );






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END init_fam() ******************************';
      write(6,*) '';
      endif

      return;
      end;
