c======================================================================c

      subroutine init_fam( lpr )

c======================================================================c

      implicit REAL*8 (a-h,o-z)
      include 'dirqfam.par'
      LOGICAL lpr;

      common /fam/ omega_start, omega_end, delta_omega, omega_print,
     &             omega, gamma_smear,
     &             i_calculation_type, i_coulomb, i_pairing,
     &             J_multipole, K_multipole, ISO;

      CHARACTER fg_spx;
      common /simplex/ N_total         , N_blocks        ,
     &                 ia_spx(NBX)     , id_spx(NBX)     ,
     &                 nf_size(NBX)    , ng_size(NBX)    ,
     &                 nz_spx(NBSX,NBX), nr_spx(NBSX,NBX),
     &                 ml_spx(NBSX,NBX), fg_spx(NBSX,NBX);

      common /DDPC1/  a_s , b_s , c_s , d_s ,
     &                a_v , b_v , c_v , d_v ,
     &                a_tv, b_tv, c_tv, d_tv,
     &                del_s,
     &                rho_sat;

      common /TMR_param/ G_pairing, a_pairing;

      common /v_energy/ E_fam_v( NBSX , NBX , 2 );

      common /u_energy/ E_fam_u( NBSX , NBX , 2 );

      common /fam_energies/ E_fam( NTX , 2 );

      LOGICAL dh_nnz, dDelta_nnz, dkappa_nnz, f_nnz;
      common /nnz_blocks/ dh_nnz    ( NBX , NBX ),
     &                    dDelta_nnz( NBX , NBX ),
     &                    dkappa_nnz( NBX , NBX ),
     &                    f_nnz     ( NBX , NBX );



      CHARACTER fg1, fg2;



      if(lpr) then
      write(6,*) '';
      write(6,*) '****** BEGIN init_fam() *****************************';
      write(6,*) '';
      endif

      write(6,'(a)') 'Initializing FAM.';
      call flush(6);






c-----Reading FAM parameters
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
      read(infam,'(18x,i9)') i_calculation_type;
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,'(18x,i9)') i_coulomb;
      read(infam,'(18x,i9)') i_pairing;
      read(infam,*);
      read(infam,*);
      read(infam,*);
      read(infam,'(18x,i9)') J_multipole;
      read(infam,'(18x,i9)') K_multipole;
      read(infam,'(18x,i9)') ISO;
      read(infam,'(18x,1f9.3)') gamma_smear;
      read(infam,*);
      read(infam,'(18x,1f9.3)') omega_start;
      read(infam,'(18x,1f9.3)') omega_end;
      read(infam,'(18x,1f9.3)') delta_omega;
      read(infam,*);
      read(infam,'(18x,1f9.3)') omega_print;
      close(infam);






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

c-----TMR separable pairing parameters
      G_pairing = 728.D0;        ![MeV/fm^3]
      a_pairing = 0.644204936D0; ![fm]






c-----Calculation of FAM energies
      do it = 1 , 2
          ie = 1;
          do ib = 1 , N_blocks
              do i = 1 , id_spx(ib)
                  if( E_fam_v(i,ib,it) .ne. E_fam_u(i,ib,it) ) then
                      stop 'Error: E_fam_v and E_fam_u wrong!';
                  endif
                  E_fam(ie,it) = E_fam_v(i,ib,it);
                  ie = ie + 1;
              enddo
          enddo
          if( ie-1 .ne. N_total ) then
               stop 'Error: ie-1 =/= N_total!';
          endif
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






c-----Calculation of basis wave-functions
      call init_basis( .false. );






c-----Calculation of multipole excitation operator
      call init_multipole( .false. );






c-----Calculation of Rcm and Pcm operators
      call init_spurious( .false. );






c-----Preparation of Coulomb (Green's function)
      call init_coulomb( .false. );






c-----Calculation of W coefficients
      call init_pairing( .false. );






      if(lpr) then
      write(6,*) '';
      write(6,*) '****** END init_fam() *******************************';
      write(6,*) '';
      endif

      return;
      end;
