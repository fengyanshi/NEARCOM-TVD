# ------------------------------------------------------------------------------
#                      Makefile for building SWAN program
# ------------------------------------------------------------------------------
#
# Before compilation, type "make config" first!
#
# To compile the serial executable type "make ser"
# To compile the OpenMP executable type "make omp"
# To compile the MPI executable type "make mpi"
#
# To remove compiled objects and modules: type "make clean"
#
# To remove compiled objects, modules, and executable: type "make allclean"
#
# Please do not change anything below, unless you're very sure what you're doing
# ------------------------------------------------------------------------------

include macros.inc

SWAN_EXE = master.exe

SWAN_INCS = \
timecomm.inc \
swcomm1.inc \
swcomm2.inc \
swcomm3.inc \
swcomm4.inc \
ocpcomm1.inc \
ocpcomm2.inc \
ocpcomm3.inc \
ocpcomm4.inc

SWAN_OBJS = \
swmod1.$(EXTO) \
sizes.$(EXTO) \
pass.$(EXTO) \
interp.$(EXTO) \
master.$(EXTO) \
sadi.$(EXTO) \
sed.$(EXTO) \
sadicirc.$(EXTO) \
master_time_ctr.$(EXTO) \
swanmain.$(EXTO) \
swanpre1.$(EXTO) \
swanpre2.$(EXTO) \
swancom1.$(EXTO) \
swancom2.$(EXTO) \
swancom3.$(EXTO) \
swancom4.$(EXTO) \
swancom5.$(EXTO) \
swancomi.$(EXTO) \
pass_in_swan.$(EXTO) \
pass_out_swan.$(EXTO) \
swanout1.$(EXTO) \
swanout2.$(EXTO) \
swanser.$(EXTO) \
swanparll.$(EXTO) \
ocpids.$(EXTO) \
ocpcre.$(EXTO) \
ocpmix.$(EXTO)

.SUFFIXES: .f .F .for

.PHONEY: help

help:
	@echo "This Makefile supports the following:"
	@echo "make config    -- makes machine-dependent macros include file"
	@echo "make ser       -- makes the Serial $(SWAN_EXE) executable"
	@echo "make omp       -- makes the OpenMP $(SWAN_EXE) executable"
	@echo "make mpi       -- makes the    MPI $(SWAN_EXE) executable"
	@echo "make clean     -- removes compiled objects and modules"
	@echo "make allclean  -- removes compiled objects, modules and $(SWAN_EXE)"

config:
	@perl platform.pl

install:
	@perl platform.pl

ser:
	@perl switch.pl $(swch) *.ftn
	$(MAKE) FOR=$(F90_SER) FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_SER)" \
                INCS="$(INCS_SER)" LIBS="$(LIBS_SER)" $(SWAN_EXE)

omp:
	@perl switch.pl $(swch) -omp *.ftn
	$(MAKE) FOR=$(F90_OMP) FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_OMP)" \
                INCS="$(INCS_OMP)" LIBS="$(LIBS_OMP)" $(SWAN_EXE)

mpi:
	@perl switch.pl $(swch) -mpi *.ftn
	$(MAKE) FOR=$(F90_MPI) FFLAGS="$(FLAGS_OPT) $(FLAGS_MSC) $(FLAGS_MPI)" \
                INCS="$(INCS_MPI)" LIBS="$(LIBS_MPI)" $(SWAN_EXE)

$(SWAN_EXE): $(SWAN_OBJS) $(SWAN_INCS)
	$(FOR) $(SWAN_OBJS) $(FFLAGS) $(OUT)$(SWAN_EXE) $(INCS) $(LIBS)

.f.o:
	$(FOR) $< -c $(FFLAGS) $(INCS)

.F.o:
	$(FOR) $< -c $(FFLAGS) $(INCS)

.for.o:
	$(FOR) $< -c $(FFLAGS) $(INCS)

.for.obj:
	$(FOR) $< -c $(FFLAGS) $(INCS)

clean:
	$(RM) *.$(EXTO) *.mod

allclean:
	$(RM) *.$(EXTO) *.mod *.f ocp*.F sw*.F *.for $(SWAN_EXE)
