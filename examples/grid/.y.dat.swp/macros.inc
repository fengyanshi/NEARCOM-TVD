##############################################################################
# IA32_Intel:		Intel Pentium with Linux using Intel compiler.
# Note: -tpp6 is for Pentium III & -tpp7 is for Pentium IV
##############################################################################
F90_SER = ifc
F90_OMP = ifc
F90_MPI = mpif90
FLAGS_OPT = -O2 -tpp7
FLAGS_MSC = -W0 -auto
FLAGS_SER =
FLAGS_OMP = -openmp -fpp2
FLAGS_MPI =
INCS_SER =
INCS_OMP =
INCS_MPI =
LIBS_SER =
LIBS_OMP =
LIBS_MPI =
OUT = -o 
EXTO = o
MAKE = make
RM = rm -f
swch = -unix -f95
