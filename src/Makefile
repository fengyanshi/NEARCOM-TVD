#-----------BEGIN MAKEFILE-------------------------------------------------
            SHELL         = /bin/sh
            DEF_FLAGS     = -P -C -traditional 
            EXEC          = nearcom
#==========================================================================
#--------------------------------------------------------------------------
#        PRECISION          DEFAULT PRECISION: SINGLE                     
#                           UNCOMMENT TO SELECT DOUBLE PRECISION
#--------------------------------------------------------------------------

#            FLAG_1 = -DDOUBLE_PRECISION 
#            FLAG_2 = -DPARALLEL
            FLAG_3 = -DSAMPLES
#             FLAG_4 = -DSEDIMENT
#             FLAG_6 = -DINTEL
#             FLAG_7 = -DRESIDUAL
             FLAG_8 = -DCURVILINEAR
#             FLAG_9 = -DSTATIONARY
#            FLAG_9 = -DSTATIONARY_RE_INITIAL !have to use with stationary
#             FLAG_10 = -DCYCLEUPDATE
#             FLAG_11 = -DONLINE
#             FLAG_12 = -DMULTISIZE
#              FLAG_13 = -DSUSPENDED_LOAD
#              FLAG_14 = -DDYE
#               FLAG_15 = -DCOUPLING
#              FLAG_16 = -DONED_IN_X
#              FLAG_17 = -DONED_IN_Y
              FLAG_18 = -DMOMENTUM_BALANCE
#               FLAG_19 = -DTRACKING
#               FLAG_20 = -DONLINE_RESIDUAL
#               FLAG_21 = -DDEBUG
#               FLAG_22 = -DCOUPLING_NO_UV
#--------------------------------------------------------------------------
#  mpi defs 
#--------------------------------------------------------------------------
         CPP      = /usr/bin/cpp 
         CPPFLAGS = $(DEF_FLAGS)
#         FC       = gfortran
          FC       = mpif90
         DEBFLGS  = 
         OPT      = 
#         OPT      = -O3  -qstrict -qtune=auto -qcache=auto -qalign=4k -w -qfixed
         CLIB     = 
#==========================================================================

         FFLAGS = $(DEBFLGS) $(OPT) 
         MDEPFLAGS = --cpp --fext=f90 --file=-
         RANLIB = ranlib
#--------------------------------------------------------------------------
#  CAT Preprocessing Flags
#--------------------------------------------------------------------------
           CPPARGS = $(CPPFLAGS) $(DEF_FLAGS) $(FLAG_1) $(FLAG_2) \
			$(FLAG_3) $(FLAG_4) $(FLAG_5) $(FLAG_6) \
			$(FLAG_7) $(FLAG_8) $(FLAG_9) $(FLAG_10)  \
			$(FLAG_11) $(FLAG_12) $(FLAG_13) $(FLAG_14) \
			$(FLAG_15) $(FLAG_16) $(FLAG_17) $(FLAG_18) \
			$(FLAG_19) $(FLAG_20) $(FLAG_21) $(FLAG_22) \
			$(FLAG_23) $(FLAG_24)
#--------------------------------------------------------------------------
#  Libraries           
#--------------------------------------------------------------------------

#            LIBS  = $(PV3LIB) $(CLIB)  $(PARLIB) $(IOLIBS) $(MPILIB) $(GOTMLIB)
#            INCS  = $(IOINCS) $(GOTMINCS)


#--------------------------------------------------------------------------
#  Preprocessing and Compilation Directives
#--------------------------------------------------------------------------
.SUFFIXES: .fsw .fsw90 .fsc .fsd

.fsw.o:
	$(CPP) $(CPPARGS) $*.fsw > $*.f
	$(FC)  $*.f -c $(FFLAGS) $(INCS)

.fsw90.o:
	$(CPP) $(CPPARGS) $*.fsw90 > $*.f90
	$(FC)  $*.f90 -c $(FFLAGS) $(INCS)
.fsc.o:
	$(CPP) $(CPPARGS) $*.fsc > $*.f90
	$(FC)  $*.f90 -c $(FFLAGS) $(INCS)
.fsd.o:
	$(CPP) $(CPPARGS) $*.fsd > $*.f90
	$(FC)  $*.f90 -c $(FFLAGS) $(INCS)
#	\rm $*.f90
#--------------------------------------------------------------------------
#  FUNWAVE-TVD Source Code.
#--------------------------------------------------------------------------

MODS  = sc_mod_param.F sc_mod_global.F sw_swmod1.F sw_swmod2.F sw_m_fileio.F sw_m_constants.F sw_serv_xnl4v5.F sw_mod_xnl4v5.F  sw_ocpcre.F sw_ocpmix.F mod_pass.F sw_common.F sc_mod_util.F sd_mod.F tracer.F\

MAIN  = coupler.F main_pass.F sw_init.F sw_swanmain.F sw_cycle.F sw_swancom1.F sw_swancom2.F sw_swancom3.F sw_swancom4.F sw_swancom5.F sw_swanout1.F sw_swanout2.F sw_swanparll.F sw_swanpre1.F sw_swanpre2.F sw_swanser.F sw_ocpids.F sc_init.F sc_main.F sc_shorecircmain.F sc_fluxes.F sc_bc.F sc_io.F sc_trid.F sd_main.F\

#MODS  = mod_param.F	mod_global.F	mod_util.F \
#MAIN  = main.F bc.F fluxes.F init.F io.F trid.F       \

SRCS = $(MODS)  $(MAIN)

OBJS = $(SRCS:.F=.o)

#--------------------------------------------------------------------------
#  Linking Directives               
#--------------------------------------------------------------------------

$(EXEC):	$(OBJS)
		$(FC) $(FFLAGS) $(LDFLAGS) -o $(EXEC) $(OBJS) $(LIBS)
#	mv $(EXEC) ../work/.
#--------------------------------------------------------------------------

#--------------------------------------------------------------------------
#  Tar Up Code                           
#--------------------------------------------------------------------------

tarfile:
	tar cvf funwave_tvd.tar *.F  Makefile

#--------------------------------------------------------------------------
#  Cleaning targets.
#--------------------------------------------------------------------------

clean:
		/bin/rm -f *.o *.mod *.f *.f90

clobber:	clean
		/bin/rm -f *.f90 *.o mytvd







