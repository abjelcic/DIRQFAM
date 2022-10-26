#If it fails to link with Intel oneAPI MKL: export LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/latest/lib/intel64
#Set number of Intel oneAPI MKL threads   : export MKL_NUM_THREADS=2
#Set number of OpenBLAS         threads   : export OPENBLAS_NUM_THREADS=2

FC = gfortran

GENFLGS = -cpp -ffree-line-length-none

LDFLAGS = -lblas -llapack
#If Intel oneAPI MKL is used put this instead:
#LDFLAGS = -I/opt/intel/oneapi/mkl/latest/include -L/opt/intel/oneapi/mkl/latest/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl

OPT = $(GENFLGS) -O3

DBG = $(GENFLGS) -D DEBUG -Og -g -Wall                              \
                                 -Wno-tabs                          \
                                 -Wno-compare-reals                 \
                                 -Wno-unused-label                  \
                                 -Wno-unused-dummy-argument         \
                                 -Wno-unused-parameter              \
                                 -Wextra                            \
                                 -Wconversion                       \
                                 -Wsurprising                       \
                                 -Wline-truncation                  \
                                 -Wcharacter-truncation             \
                                 -Wpedantic                         \
                                 -fcheck=all                        \
                                 -fbacktrace                        \
                                 -finit-real=nan                    \
                                 -finit-integer=-999999             \
                                 -finit-logical=true                \
                                 -ftrapv                            \
                                 -fno-fast-math                     \
                                 -fno-omit-frame-pointer            \
                                 -fsanitize=address                 \
                                 -fstack-protector                  \

SRCMOD =                       \
./src/modules/dirhbpar.f90     \
./src/modules/modules_grs.f90  \
./src/modules/modules_qfam.f90 \

SRCGRS =                       \
./src/dirhbz.f                 \

SRCFAM =                       \
./src/utility.f90              \
./src/main_fam.f90             \
./src/read_input.f90           \
./src/base_simplex.f90         \
./src/construct_uv.f90         \
./src/check_gs_dens.f90        \
./src/check_unitarity.f90      \
./src/init_fam.f90             \
./src/init_basis.f90           \
./src/init_gs.f90              \
./src/init_multipole.f90       \
./src/init_spurious.f90        \
./src/init_mesons.f90          \
./src/init_coulomb.f90         \
./src/init_pairing.f90         \
./src/start_fam.f90            \
./src/iter_fam.f90             \
./src/fam_dh20dh02dh11.f90     \
./src/fam_xy.f90               \
./src/fam_drhodkappa.f90       \
./src/fam_ddensdcurr.f90       \
./src/fam_dpotentials.f90      \
./src/fam_dmesons.f90          \
./src/fam_dcoulomb.f90         \
./src/fam_dh.f90               \
./src/fam_ddelta.f90           \
./src/fam_gmres.f90            \
./src/fam_spurious.f90         \
./src/fam_strength.f90	       \
./src/printout.f90             \
./src/nuclocfunc.f90           \
./src/kpm.f90                  \

MODLOC =                       \
-I./src/modules/mods           \
-J./src/modules/mods           \

run: $(SRCMOD) $(SRCGRS) $(SRCFAM)
	$(FC) $(MODLOC) $(OPT) $(SRCMOD) $(SRCGRS) $(SRCFAM) $(LDFLAGS) -o run

dbg: $(SRCMOD) $(SRCGRS) $(SRCFAM)
	$(FC) $(MODLOC) $(DBG) $(SRCMOD) $(SRCGRS) $(SRCFAM) $(LDFLAGS) -o dbg

clear:
	rm -f run dbg ./*.out ./src/modules/mods/*.mod ./output/GS_output/* ./output/QFAM_output/*
