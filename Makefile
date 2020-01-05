FC = gfortran

PREPROC = -cpp

MCMODEL = -mcmodel=large

LDFLAGS = -lblas -llapack

OPT = $(PREPROC) -O3 $(MCMODEL) -march=native

DEBUG = $(PREPROC) -Og -g $(MCMODEL) -Wall                              \
                                     -Wextra                            \
                                     -Wconversion                       \
                                     -Wno-unused-parameter              \
                                     -Warray-temporaries                \
                                     -Wsurprising                       \
                                     -Wline-truncation                  \
                                     -Wcharacter-truncation             \
                                     -fcheck=all                        \
                                     -fbacktrace                        \
                                     -ffpe-trap=zero,overflow,underflow \
                                     -fimplicit-none                    \



SRCGS =                      \
./source/dirhbz.f            \

SRCFAM =                     \
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
./source/fam_h20h02.f        \
./source/fam_xy.f            \
./source/fam_drho.f          \
./source/fam_dkappa.f        \
./source/fam_ddensdcurr.f    \
./source/fam_dpotentials.f   \
./source/fam_dmesons.f       \
./source/fam_dcoulomb.f      \
./source/fam_dh1.f           \
./source/fam_ddelta.f        \
./source/fam_broyden.f       \
./source/fam_spurious.f      \
./source/fam_strength.f      \
./source/printout.f          \
./source/utility_functions.f \



run: prep $(SRCGS) $(SRCFAM)
	$(FC) $(OPT) $(SRCGS) $(SRCFAM) $(LDFLAGS) -o run

prep: ./dirqfam.dat
	$(FC) $(DEBUG) ./source/prep.f -o run_prep && ./run_prep && rm run_prep

debug: prep $(SRCGS) $(SRCFAM)
	$(FC) -c        $(OPT)  -fcheck=all $(SRCGS)       $(LDFLAGS) -o gs.o && \
	$(FC) -D DEBUG  $(DEBUG)            $(SRCFAM) gs.o $(LDFLAGS) -o run  && \
	rm -f *.o

clear:
	rm -f *.o run
