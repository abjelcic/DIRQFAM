FC1 = gfortran ./source/prep.f -o run && ./run && rm run

FC2 = gfortran -O3 -mcmodel=medium -funroll-loops

OBJ =                        \
./source/dirhbz.f            \
./source/main_fam.f          \
./source/base_simplex.f      \
./source/construct_uv.f      \
./source/check_gs_dens.f     \
./source/check_unitarity.f   \
./source/init_fam.f          \
./source/init_basis.f        \
./source/init_gs.f           \
./source/init_multipole.f    \
./source/init_spurious.f     \
./source/init_mesons.f       \
./source/init_coulomb.f      \
./source/init_pairing.f      \
./source/start_fam.f         \
./source/iter_fam.f          \
./source/fam_drhodkappa.f    \
./source/fam_ddensdcurr.f    \
./source/fam_dpotentials.f   \
./source/fam_dmesons.f       \
./source/fam_dcoulomb.f      \
./source/fam_dh1.f           \
./source/fam_ddelta.f        \
./source/fam_broyden.f       \
./source/fam_strength.f      \
./source/printout.f          \
./source/utility_functions.f \

prep:
	$(FC1)

run: $(OBJ)
	$(FC2) $(OBJ) -lblas -llapack -o run
