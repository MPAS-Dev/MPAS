.SUFFIXES: .F .o

OBJS = mpas_ocn_analysis_driver.o \
       mpas_ocn_global_stats.o \
       mpas_ocn_zonal_mean.o \
       mpas_ocn_paraview_catalyst.o

all: $(OBJS)

mpas_ocn_analysis_driver.o: mpas_ocn_global_stats.o mpas_ocn_zonal_mean.o mpas_ocn_paraview_catalyst.o

mpas_ocn_global_stats.o:

mpas_ocn_zonal_mean.o:

mpas_ocn_paraview_catalyst.o:

clean:
	$(RM) *.o *.i *.mod *.f90 lib*.a 

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) 
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) 
endif
