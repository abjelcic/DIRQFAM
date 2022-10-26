!======================================================================!

      subroutine print_header( tape  )

!======================================================================!
      use fam_input;
      use basis;
      implicit none;
      integer , intent(in) :: tape;

      double precision :: ea, rms, betg, gamg; ! From the ground state part of the code.
      common /erwar /     ea, rms, betg, gamg; ! betg contains the ground state self-consistent deformation parameter.

          write(tape,'(a)') '';
          write(tape,'(a,a,a,i0,a)') 'Nucleus: ' , nucleusName , ' ' , nucleusZ+nucleusN , '.';
          write(tape,'(a,a,a)') 'Parametrization of the Lagrangian: ' , LagrangianModelName , '.';
          write(tape,'(a,sp,f6.3,a)') 'Ground state deformation beta = ' , betg , '.';
          write(tape,'(a,a,i0,a,i0,a)') merge('Isoscalar','Isovector',Isospin==0) , ' J = ' , J_multipole , ', K = ' , K_multipole , ' excitation.';
          if( any(calculation_type==[3,4]) ) then
              write(tape,'(a)') 'Gamma smearing = N/A.';
          else
              write(tape,'(a,f5.3,a)') 'Gamma smearing = ' , gamma_smear , ' [MeV].';
          end if
          write(tape,'(a)') merge( 'Free response.                 ' ,                     &
                                   'Fully self-consistent response.' , calculation_type==0 );
          write(tape,'(a)') '';


          write(tape,'(a,i0,a,i0,a,a)') 'n0f = ' , n0f , '(',n0f+1,')' , ' oscillator shells for large(small) spinor component.';
          if( LagrangianModelName=='DD-ME2' ) write(tape,'(a,i0,a)') 'n0b = ' , n0b , ' oscillator shells for meson fields.';
          write(tape,'(a,i0,a,i0,a)') 'Gaussian quadrature mesh: NGH = ' , NGH , ', NGL = ' , NGL , '.';
          write(tape,'(a,sp,f6.3,a)') 'Oscillator basis deformation beta0 = ' , beta0 , '.';
          write(tape,'(a)') merge( 'Coulomb included in calculation.  ' ,                    &
                                   'Coulomb excluded from calculation.' , include_coulomb==1 );
          write(tape,'(a)') merge( 'Pairing included in calculation.  ' ,                    &
                                   'Pairing excluded from calculation.' , include_pairing==1 );

          write(tape,'(a)') '';
          call flush(tape);

      return;
      end subroutine print_header






!======================================================================!

      subroutine print_basis()

!======================================================================!
      use fam_input;
      use simplex;
      use basis;
      implicit none;

      integer   :: tape;
      integer   :: ib, i, il;
      character :: fg;
      integer   :: nz;
      integer   :: nr;
      integer   :: ml;

          open( newunit=tape , file='./output/QFAM_output/basis.out' , action='write' );

          write(tape,'(a)') 'Define a QHO 2D spinor: |nz,nr,ml,ms>(z,r,phi) = phi_{nz}(z,bz) * phi_{nr}^{|ml|}(r,br) * exp(i*ml*phi)/sqrt(2*pi) * X_{ms}(1/2)   ,';
          write(tape,'(a)') 'where phi_{nz}(z,bz)        = 1/sqrt(bz) * 1/sqrt(2^nz*nz!*sqrt(pi)) *                   H_nz(z/bz)          * exp(-0.5*(z/bz)^2)  ,';
          write(tape,'(a)') '      phi_{nr}^{|ml|}(r,br) = 1/br       * sqrt(2*nr!/(nr+|ml|)!)    * (r/br)^{|ml|}   * L_{nr}^{|ml|}(r/br) * exp(-0.5*(r/br)^2)  ,';
          write(tape,'(a)') '      X_{ms=+1/2}(1/2) = [ 1 ; 0 ]                                                                                                 ,';
          write(tape,'(a)') '      X_{ms=-1/2}(1/2) = [ 0 ; 1 ]                                                                                                 .';
          write(tape,'(a)') '                                                                                                                                    ';
          write(tape,'(a)') 'Define simplex-y vectors: |nz,nr,ml,s=+i> = i/sqrt(2)*|nz,nr,+ml,ms=+1/2> + 1/sqrt(2)*|nz,nr,-ml,ms=-1/2>                          ,';
          write(tape,'(a)') '                          |nz,nr,ml,s=-i> = 1/sqrt(2)*|nz,nr,+ml,ms=+1/2> + i/sqrt(2)*|nz,nr,-ml,ms=-1/2>                          .';
          write(tape,'(a)') '                                                                                                                                    ';
          write(tape,'(a)') 'We define the basis B=[B1,B2] used in DIRQFAM code:                                          ';
          write(tape,'(a)') ' B1 = [ all vectors |fg,nz,nr,ml,s=+i> such that fg={f,g}, and nz+2*nr+|ml|<=n0f if fg=f, and nz+2*nr+|ml|<=n0f+1 if fg=g ]        ,';
          write(tape,'(a)') ' B2 = [ all vectors |fg,nz,nr,ml,s=-i> such that fg={f,g}, and nz+2*nr+|ml|<=n0f if fg=f, and nz+2*nr+|ml|<=n0f+1 if fg=g ]        ,';
          write(tape,'(a)') 'where the order in B1 and B2 are the same.                                                                                          ';
          write(tape,'(a)') '                                                                                                                                    ';
          write(tape,'(a)') '4-dimensional spinors |fg,nz,nr,ml,s> are given by:                                                                                 ';
          write(tape,'(a)') '   |fg=f,nz,nr,ml,s=+i> = [ |nz,nr,ml,+i> ; [0;0]          ]                                                                       ,';
          write(tape,'(a)') '   |fg=g,nz,nr,ml,s=+i> = [ [0;0]         ; i|nz,nr,ml,-i> ]                                                                       ,';
          write(tape,'(a)') '   |fg=f,nz,nr,ml,s=-i> = [ |nz,nr,ml,-i> ; [0;0]          ]                                                                       ,';
          write(tape,'(a)') '   |fg=g,nz,nr,ml,s=-i> = [ [0;0]         ; i|nz,nr,ml,+i> ]                                                                       .';
          write(tape,'(a)') '                                                                                                                                    ';
          write(tape,'(a)') 'Below we print the order of vectors (fg,nz,nr,ml) in block B1. Order is the same in block B2.                                       ';
          write(tape,'(a)') 'If fg=f, then we include only nz+2*nr+|ml|<=n0f, while if fg=g we include only nz+2*nr+|ml|<=n0f+1.                                 ';
          write(tape,'(a)') '                                                                                                                                    ';

          write(tape,'(a,i0,a,i0,a)') 'Total of ' , N_total , ' basis vectors are in B1 for number of shells n0f = ' , n0f , '.';
          write(tape,'(a,f9.5,a,f9.5,a,/)') 'bz = ' , bz , ' [fm^3], br = ' , bp , ' [fm^3].';


          write(tape,'(5x,a,5x,a,4x,a,4x,a,4x,a,/)') 'k' , 'fg' , 'nz' , 'nr' , 'ml';
          do il = 1 , N_total

              ib = getBlock(il);
              i  = getIndexWithinBlock(il);

              fg = fg_spx(ib)%index(i);
              nz = nz_spx(ib)%index(i);
              nr = nr_spx(ib)%index(i);
              ml = ml_spx(ib)%index(i);

              write(tape,'(i6,5x,a,i6,i6,i6)') il , fg , nz , nr , ml;

          end do

          close(tape);

      return;
      end subroutine print_basis






!======================================================================!

      subroutine print_qpEnergies()

!======================================================================!
      use fam_input;
      use simplex;
      use fam_energies;
      implicit none;

      integer :: tape;
      integer :: it, il, ib, i;

      do it = 1 , 2

          if( it == 1 ) open( newunit=tape , file='./output/QFAM_output/qpenergy_neut.out', action='write' );
          if( it == 2 ) open( newunit=tape , file='./output/QFAM_output/qpenergy_prot.out', action='write' );

          write(tape,'(a,i0,a)') 'Quasiparticle energies are given in two simplex blocks: [ E , E ]. There are ' , N_total , ' energies in E.';
          write(tape,'(a)')      'In simplex-y basis, there holds U=[u,0;0,conjg(u)] and V=[0,-conjg(v);v,0], where u_{k,mu} and v_{k,mu} are printed in u.out and v.out.';
          write(tape,'(a,/)')    'Basis vectors k are given in basis.out and here we give quasiparticle energies E_mu [MeV].';

          write(tape,'(5x,a,5x,a,/)') 'mu' , 'E_mu';

          do il = 1 , N_total
              ib = getBlock(il);
              i  = getIndexWithinBlock(il);

              write(tape,'(i6,e17.7E3)') il , E_fam(it)%blocks(ib)%vec(i);
          end do

          close(tape);

      enddo

      return;
      end subroutine print_qpEnergies






!======================================================================!

      subroutine print_uv()

!======================================================================!
      use fam_input;
      use simplex;
      use u_matrix;
      use v_matrix;
      implicit none;

      integer :: tape;
      integer :: it, ib1, ib2, i1, i2;

      ! Printing u.
      do it = 1 , 2

          if( it == 1 ) open( newunit=tape , file='./output/QFAM_output/U_neut.out' , action='write' );
          if( it == 2 ) open( newunit=tape , file='./output/QFAM_output/U_prot.out' , action='write' );

          write(tape,'(a,i0,a,i0,a)') 'U = [ u , 0 ; 0 , conjg(u) ]. Here ',N_total,'x',N_total,' matrix u_{k,mu} is printed. Other elements of u_{k,mu} not printed here are zero.';

          call print_header( tape );

          write(tape,'(/,/,4x,a,7x,a,11x,a,/)') 'k' , 'mu' , 'u_{k,mu}';

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( u(it)%nnzblocks(ib1,ib2) ) then
                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , id_spx(ib1)

                              write(tape,'(i6,2x,i6,2x,e17.7E3,a,e17.7E3,a)') ia_spx(ib1)+i1-1 , ia_spx(ib2)+i2-1 , real( u(it)%blocks(ib1,ib2)%mat(i1,i2) ) , ' + ' , imag( u(it)%blocks(ib1,ib2)%mat(i1,i2) ) , ' i';

                          end do
                      end do
                  end if
              end do
          end do

          close(tape);

      end do

      ! Printing v.
      do it = 1 , 2

          if( it == 1 ) open( newunit=tape , file='./output/QFAM_output/V_neut.out' , action='write' );
          if( it == 2 ) open( newunit=tape , file='./output/QFAM_output/V_prot.out' , action='write' );

          write(tape,'(a,i0,a,i0,a)') 'V = [ 0 , -conjg(v) ; v , 0 ]. Here ',N_total,'x',N_total,' matrix v_{k,mu} is printed. Other elements of v_{k,mu} not printed here are zero.';

          call print_header( tape );

          write(tape,'(/,/,4x,a,7x,a,11x,a,/)') 'k' , 'mu' , 'v_{k,mu}';

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( v(it)%nnzblocks(ib1,ib2) ) then
                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , id_spx(ib1)

                              write(tape,'(i6,2x,i6,2x,e17.7E3,a,e17.7E3,a)') ia_spx(ib1)+i1-1 , ia_spx(ib2)+i2-1 , real( v(it)%blocks(ib1,ib2)%mat(i1,i2) ) , ' + ' , imag( v(it)%blocks(ib1,ib2)%mat(i1,i2) ) , ' i';

                          end do
                      end do
                  end if
              end do
          end do

          close(tape);

      end do

      return;
      end subroutine print_uv






!======================================================================!

      subroutine print_densities( omega )

!======================================================================!
      use fam_input;
      use quadrature;
      use gs_dens;
      use ddens;
      implicit none;

      double precision , intent(in) :: omega;
      integer                       :: tape;
      integer                       :: it, ih, il;
      double precision              :: z, r;

          do it = 1 , 2
              if( it == 1 ) open( newunit=tape , file='./output/QFAM_output/rhov_neut.out' , action='write' );
              if( it == 2 ) open( newunit=tape , file='./output/QFAM_output/rhov_prot.out' , action='write' );

              write(tape,'(a)') 'rho_v(z,r,phi,t) = rho0_v(z,r) + 2*eta*Re[ exp(-i*omega*t) * drho_v(z,r,phi) ] + O(eta^2).';
              write(tape,'(a,i0,a)') 'drho_v(z,r,phi) = drho_v(z,r) * cos(' , K_multipole , '*phi).';
              call print_header(tape);
              write(tape,'(a,e16.9,a)') 'omega = ' , omega , ' [MeV].';

              write(tape,'(/,/,3x,a,6x,a,7x,a,6x,a,/)') 'r[fm]' , 'z[fm]' , 'rho0_v(z,r)[fm^-3]' , 'drho_v(z,r)[fm^-3]';
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;

                      z = sign(1,ih) * zb_fam(abs(ih));
                      r = rb_fam(il);

                      write(tape,'(f10.5,2x,f10.5,2x,e16.7E3,8x,e16.7E3,a,e16.7E3,a)')  r , z , rhov_GS(ih,il,it) , real(drho_v(ih,il,it)) , '  +' , imag(drho_v(ih,il,it)) , ' i';

                  end do
              end do

              close(tape);
          end do

      return;
      end subroutine print_densities






!======================================================================!

      subroutine print_currents( omega )

!======================================================================!
      use fam_input;
      use quadrature;
      use dcurr;
      implicit none;

      double precision , intent(in) :: omega;
      integer                       :: tape;
      integer                       :: it, ih, il;
      double precision              :: z, r;

          ! Printing j_z.
          do it = 1 , 2
              if( it == 1 ) open( newunit=tape , file='./output/QFAM_output/jz_neut.out' , action='write' );
              if( it == 2 ) open( newunit=tape , file='./output/QFAM_output/jz_prot.out' , action='write' );

              write(tape,'(a)') 'j_z(z,r,phi,t) = 2*eta*Re[ exp(-i*omega*t) * dj_z(z,r,phi) ] + O(eta^2).';
              write(tape,'(a,i0,a)') 'dj_z(z,r,phi) = dj_z(z,r) * cos(' , K_multipole , '*phi).';
              call print_header(tape);
              write(tape,'(a,e16.9,a)') 'omega = ' , omega , ' [MeV].';

              write(tape,'(/,/,3x,a,6x,a,7x,a,/)') 'r[fm]' , 'z[fm]' , 'dj_z(z,r)[cfm^-3]';
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;

                      z = sign(1,ih) * zb_fam(abs(ih));
                      r = rb_fam(il);

                      write(tape,'(f10.5,2x,f10.5,2x,e16.7E3,a,e16.7E3,a)')  r , z , real(dj_z(ih,il,it)) , '  +' , imag(dj_z(ih,il,it)) , ' i';

                  end do
              end do

              close(tape);
          end do

          ! Printing j_r.
          do it = 1 , 2
              if( it == 1 ) open( newunit=tape , file='./output/QFAM_output/jr_neut.out' , action='write' );
              if( it == 2 ) open( newunit=tape , file='./output/QFAM_output/jr_prot.out' , action='write' );

              write(tape,'(a)') 'j_r(z,r,phi,t) = 2*eta*Re[ exp(-i*omega*t) * dj_r(z,r,phi) ] + O(eta^2).';
              write(tape,'(a,i0,a)') 'dj_r(z,r,phi) = dj_r(z,r) * cos(' , K_multipole , '*phi).';
              call print_header(tape);
              write(tape,'(a,e16.9,a)') 'omega = ' , omega , ' [MeV].';

              write(tape,'(/,/,3x,a,6x,a,7x,a,/)') 'r[fm]' , 'z[fm]' , 'dj_r(z,r)[cfm^-3]';
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;

                      z = sign(1,ih) * zb_fam(abs(ih));
                      r = rb_fam(il);

                      write(tape,'(f10.5,2x,f10.5,2x,e16.7E3,a,e16.7E3,a)')  r , z , real(dj_r(ih,il,it)) , '  +' , imag(dj_r(ih,il,it)) , ' i';

                  end do
              end do

              close(tape);
          end do

          ! Printing j_phi.
          do it = 1 , 2
              if( it == 1 ) open( newunit=tape , file='./output/QFAM_output/jphi_neut.out' , action='write' );
              if( it == 2 ) open( newunit=tape , file='./output/QFAM_output/jphi_prot.out' , action='write' );

              write(tape,'(a)') 'j_phi(z,r,phi,t) = 2*eta*Re[ exp(-i*omega*t) * dj_phi(z,r,phi) ] + O(eta^2).';
              write(tape,'(a,i0,a)') 'dj_phi(z,r,phi) = dj_phi(z,r) * sin(' , K_multipole , '*phi).';
              call print_header(tape);
              write(tape,'(a,e16.9,a)') 'omega = ' , omega , ' [MeV].';

              write(tape,'(/,/,3x,a,6x,a,7x,a,/)') 'r[fm]' , 'z[fm]' , 'dj_phi(z,r)[cfm^-3]';
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;

                      z = sign(1,ih) * zb_fam(abs(ih));
                      r = rb_fam(il);

                      write(tape,'(f10.5,2x,f10.5,2x,e16.7E3,a,e16.7E3,a)')  r , z , real(dj_p(ih,il,it)) , '  +' , imag(dj_p(ih,il,it)) , ' i';

                  end do
              end do

              close(tape);
          end do

      return;
      end subroutine print_currents






!======================================================================!

      subroutine print_nuclocfunc( omega )

!======================================================================!
      use fam_input;
      use quadrature;
      use nuclearlocfunc;
      implicit none;

      double precision , intent(in) :: omega;
      integer                       :: tape;
      integer                       :: it, ih, il;
      double precision              :: z, r;

          do it = 1 , 2
              if( it == 1 ) open( newunit=tape , file='./output/QFAM_output/nuclocfunc_neut.out' , action='write' );
              if( it == 2 ) open( newunit=tape , file='./output/QFAM_output/nuclocfunc_prot.out' , action='write' );

              write(tape,'(a)') 'C(z,r,phi,t) = C0(z,r) + 2*eta*Re[ exp(-i*omega*t) * dC(z,r,phi) ] + O(eta^2).';
              write(tape,'(a,i0,a)') 'dC(z,r,phi) = dC(z,r) * cos(' , K_multipole , '*phi).';
              call print_header(tape);
              write(tape,'(a,e16.9,a)') 'omega = ' , omega , ' [MeV].';

              write(tape,'(/,/,3x,a,6x,a,7x,a,17x,a,/)') 'r[fm]' , 'z[fm]' , 'C0(z,r)' , 'dC(z,r)';
              do il = 1 , NGL
                  do ih = -NGH , +NGH
                      if( ih == 0 ) cycle;

                      z = sign(1,ih) * zb_fam(abs(ih));
                      r = rb_fam(il);

                      write(tape,'(f10.5,2x,f10.5,2x,e16.7E3,8x,e16.7E3,a,e16.7E3,a)')  r , z , C0(ih,il,it) , real(dC(ih,il,it)) , '  +' , imag(dC(ih,il,it)) , ' i';

                  end do
              end do

              close(tape);
          end do

      return;
      end subroutine print_nuclocfunc






!======================================================================!

      subroutine print_xy( omega )

!======================================================================!
      use fam_input;
      use simplex;
      use xyfam;
      implicit none;

      double precision , intent(in) :: omega;
      integer                       :: tape;
      integer                       :: it, ib1, ib2, i1, i2;

      ! Printing x.
      do it = 1 , 2

          if( it == 1 ) open( newunit=tape , file='./output/QFAM_output/X_neut.out' , action='write' );
          if( it == 2 ) open( newunit=tape , file='./output/QFAM_output/X_prot.out' , action='write' );

          write(tape,'(a,i0,a,i0,a)') 'X(omega) = [ 0 , x(omega) ; -x(omega)^T , 0 ]. Here ',N_total,'x',N_total,' matrix x_{mu,nu}(omega) is printed. Other elements of x(omega) not printed here are zero.';

          call print_header( tape );

          write(tape,'(a,e16.9,a)') 'omega = ' , omega , ' [MeV].';

          write(tape,'(/,/,4x,a,6x,a,11x,a,/)') 'mu' , 'nu' , 'x_{mu,nu}(omega)';

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( x_fam(it)%nnzblocks(ib1,ib2) ) then
                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , id_spx(ib1)

                              write(tape,'(i6,2x,i6,2x,e17.7E3,a,e17.7E3,a)') ia_spx(ib1)+i1-1 , ia_spx(ib2)+i2-1 , real( x_fam(it)%blocks(ib1,ib2)%mat(i1,i2) ) , ' + ' , imag( x_fam(it)%blocks(ib1,ib2)%mat(i1,i2) ) , ' i';

                          end do
                      end do
                  end if
              end do
          end do

          close(tape);

      end do

      ! Printing y.
      do it = 1 , 2

          if( it == 1 ) open( newunit=tape , file='./output/QFAM_output/Y_neut.out' , action='write' );
          if( it == 2 ) open( newunit=tape , file='./output/QFAM_output/Y_prot.out' , action='write' );

          write(tape,'(a,i0,a,i0,a)') 'Y(omega) = [ 0 , y(omega) ; -y(omega)^T , 0 ]. Here ',N_total,'x',N_total,' matrix y_{mu,nu}(omega) is printed. Other elements of y(omega) not printed here are zero.';

          call print_header( tape );

          write(tape,'(a,e16.9,a)') 'omega = ' , omega , ' [MeV].';

          write(tape,'(/,/,4x,a,6x,a,11x,a,/)') 'mu' , 'nu' , 'y_{mu,nu}(omega)';

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( y_fam(it)%nnzblocks(ib1,ib2) ) then
                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , id_spx(ib1)

                              write(tape,'(i6,2x,i6,2x,e17.7E3,a,e17.7E3,a)') ia_spx(ib1)+i1-1 , ia_spx(ib2)+i2-1 , real( y_fam(it)%blocks(ib1,ib2)%mat(i1,i2) ) , ' + ' , imag( y_fam(it)%blocks(ib1,ib2)%mat(i1,i2) ) , ' i';

                          end do
                      end do
                  end if
              end do
          end do

          close(tape);

      end do

      return;
      end subroutine print_xy






!======================================================================!

      subroutine print_dh20dh02( omega )

!======================================================================!
      use fam_input;
      use simplex;
      use dh02dh20matrix;
      implicit none;

      double precision , intent(in) :: omega;
      integer                       :: tape;
      integer                       :: it, ib1, ib2, i1, i2;

      ! Printing dh20.
      do it = 1 , 2

          if( it == 1 ) open( newunit=tape , file='./output/QFAM_output/dH20_neut.out' , action='write' );
          if( it == 2 ) open( newunit=tape , file='./output/QFAM_output/dH20_prot.out' , action='write' );

          write(tape,'(a,i0,a,i0,a)') 'dH20(omega) = [ 0 , dh20(omega) ; -dh20(omega)^T , 0 ]. Here ',N_total,'x',N_total,' matrix dh20_{mu,nu}(omega) is printed. Other elements of dh20(omega) not printed here are zero.';

          call print_header( tape );

          write(tape,'(a,e16.9,a)') 'omega = ' , omega , ' [MeV].';

          write(tape,'(/,/,4x,a,6x,a,11x,a,/)') 'mu' , 'nu' , 'dh20_{mu,nu}(omega)';

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh20(it)%nnzblocks(ib1,ib2) ) then
                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , id_spx(ib1)

                              write(tape,'(i6,2x,i6,2x,e17.7E3,a,e17.7E3,a)') ia_spx(ib1)+i1-1 , ia_spx(ib2)+i2-1 , real( dh20(it)%blocks(ib1,ib2)%mat(i1,i2) ) , ' + ' , imag( dh20(it)%blocks(ib1,ib2)%mat(i1,i2) ) , ' i';

                          end do
                      end do
                  end if
              end do
          end do

          close(tape);

      end do

      ! Printing dh02.
      do it = 1 , 2

          if( it == 1 ) open( newunit=tape , file='./output/QFAM_output/dH02_neut.out' , action='write' );
          if( it == 2 ) open( newunit=tape , file='./output/QFAM_output/dH02_prot.out' , action='write' );

          write(tape,'(a,i0,a,i0,a)') 'dH02(omega) = [ 0 , dh02(omega) ; -dh02(omega)^T , 0 ]. Here ',N_total,'x',N_total,' matrix dh02_{mu,nu}(omega) is printed. Other elements of dh02(omega) not printed here are zero.';

          call print_header( tape );

          write(tape,'(a,e16.9,a)') 'omega = ' , omega , ' [MeV].';

          write(tape,'(/,/,4x,a,6x,a,11x,a,/)') 'mu' , 'nu' , 'dh02_{mu,nu}(omega)';

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  if( dh02(it)%nnzblocks(ib1,ib2) ) then
                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , id_spx(ib1)

                              write(tape,'(i6,2x,i6,2x,e17.7E3,a,e17.7E3,a)') ia_spx(ib1)+i1-1 , ia_spx(ib2)+i2-1 , real( dh02(it)%blocks(ib1,ib2)%mat(i1,i2) ) , ' + ' , imag( dh02(it)%blocks(ib1,ib2)%mat(i1,i2) ) , ' i';

                          end do
                      end do
                  end if
              end do
          end do

          close(tape);

      end do

      return;
      end subroutine print_dh20dh02






!======================================================================!

      subroutine print_dh11( omega )

!======================================================================!
      use fam_input;
      use simplex;
      use dh11matrix;
      implicit none;

      double precision , intent(in) :: omega;
      integer                       :: tape;
      integer                       :: it, ib1, ib2, i1, i2;

      ! Printing dh11_1 and dh11_2.
      do it = 1 , 2

          if( it == 1 ) open( newunit=tape , file='./output/QFAM_output/dH11_neut.out' , action='write' );
          if( it == 2 ) open( newunit=tape , file='./output/QFAM_output/dH11_prot.out' , action='write' );

          write(tape,'(a,i0,a,i0,a)') 'dH11(omega) = [ dh11_1(omega) , 0 ; 0 , dh11_2(omega) ]. Here ',N_total,'x',N_total,' matrices (dh11_1)_{mu,nu}(omega) and (dh11_2)_{mu,nu}(omega) are printed. Other elements of dh11_1(omega) and dh11_2(omega) not printed here are zero.';

          call print_header( tape );

          write(tape,'(a,e16.9,a)') 'omega = ' , omega , ' [MeV].';

          write(tape,'(/,/,4x,a,6x,a,11x,a,26x,a/)') 'mu' , 'nu' , '(dh11_1)_{mu,nu}(omega)' , '(dh11_2)_{mu,nu}(omega)';

          do ib2 = 1 , N_blocks
              do ib1 = 1 , N_blocks
                  call assert( dh11_1(it)%nnzblocks(ib1,ib2) .eqv. dh11_2(it)%nnzblocks(ib1,ib2) , 'dh11_1 and dh11_2 nnzblocks error.' );
                  if( dh11_1(it)%nnzblocks(ib1,ib2) .and. dh11_2(it)%nnzblocks(ib1,ib2) ) then
                      do i2 = 1 , id_spx(ib2)
                          do i1 = 1 , id_spx(ib1)

                              write(tape,'(i6,2x,i6,2x,e17.7E3,a,e17.7E3,a,10x,e17.7E3,a,e17.7E3,a)') ia_spx(ib1)+i1-1 , ia_spx(ib2)+i2-1 , real( dh11_1(it)%blocks(ib1,ib2)%mat(i1,i2) ) , ' + ' , imag( dh11_1(it)%blocks(ib1,ib2)%mat(i1,i2) ) , ' i' &
                                                                                                                                          , real( dh11_2(it)%blocks(ib1,ib2)%mat(i1,i2) ) , ' + ' , imag( dh11_2(it)%blocks(ib1,ib2)%mat(i1,i2) ) , ' i' ;
                          end do
                      end do
                  end if
              end do
          end do

          close(tape);

      end do

      return;
      end subroutine print_dh11
