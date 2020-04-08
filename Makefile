FC = gfortran

PREPROC = -cpp

LDFLAGS = -lblas -llapack

OPT = $(PREPROC) -O3 -march=native

DEBUG = $(PREPROC) -Og -g -Wall                              \
                          -Wno-tabs                          \
                          -Wno-compare-reals                 \
                          -Wno-unused-label                  \
                          -Wno-unused-dummy-argument         \
                          -Wno-unused-parameter              \
                          -Wextra                            \
                          -Wconversion                       \
                          -Warray-temporaries                \
                          -Wsurprising                       \
                          -Wline-truncation                  \
                          -Wcharacter-truncation             \
                          -fcheck=all                        \
                          -fbacktrace                        \
#                         -ffpe-trap=zero,overflow,underflow \
#                         -fimplicit-none                    \
#                         -Wpedantic                         \

SRCPAR =                       \
./source/modules/dirqfampar.f  \

SRCMOD1 =                      \
./source/modules/modules_grs.f \

SRCMOD2 =					   \
./source/modules/modules_qfam.f \

SRCGRS =                       \
./source/dirhbz.f              \

SRCFAM =                       \
./source/main_fam.f            \
./source/base_simplex.f        \
./source/construct_uv.f        \
./source/check_gs_dens.f       \
./source/check_unitarity.f     \
./source/init_fam.f            \
./source/init_basis.f          \
./source/init_gs.f             \
./source/init_multipole.f      \
./source/init_spurious.f       \
./source/init_mesons.f         \
./source/init_coulomb.f        \
./source/init_pairing.f        \
./source/start_fam.f           \
./source/iter_fam.f            \
./source/fam_h20h02.f          \
./source/fam_xy.f              \
./source/fam_drho.f            \
./source/fam_dkappa.f          \
./source/fam_ddensdcurr.f      \
./source/fam_dpotentials.f     \
./source/fam_dmesons.f         \
./source/fam_dcoulomb.f        \
./source/fam_dh1.f             \
./source/fam_ddelta.f          \
./source/fam_broyden.f         \
./source/fam_spurious.f        \
./source/fam_strength.f        \
./source/printout.f            \
./source/utility_functions.f   \

MODLOC =				       \
-I./source/modules/mods        \
-J./source/modules/mods        \

run: $(SRCPAR) $(SRCMOD1) $(SRCMOD2) $(SRCGRS) $(SRCFAM)
	$(FC) $(MODLOC) $(OPT) $(SRCPAR) $(SRCMOD1) $(SRCMOD2) $(SRCGRS) $(SRCFAM) $(LDFLAGS) -o run

debug: $(SRCPAR) $(SRCMOD1) $(SRCMOD2) $(SRCGRS) $(SRCFAM)
	$(FC) -D DEBUG $(MODLOC) $(DEBUG) $(SRCPAR) $(SRCMOD1) $(SRCMOD2) $(SRCGRS) $(SRCFAM) $(LDFLAGS) -o dbg

clear:
	rm -f *.o run dbg ./source/modules/mods/*.mod
