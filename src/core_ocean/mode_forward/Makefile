.SUFFIXES: .F .o

OBJS = mpas_ocn_mpas_core.o \
       mpas_ocn_time_integration.o \
       mpas_ocn_time_integration_rk4.o \
       mpas_ocn_time_integration_split.o

all: core_ocean

core_ocean: $(OBJS) 

mpas_ocn_time_integration.o: mpas_ocn_time_integration_rk4.o mpas_ocn_time_integration_split.o

mpas_ocn_time_integration_rk4.o: 

mpas_ocn_time_integration_split.o: 

mpas_ocn_mpas_core.o: mpas_ocn_time_integration.o \
                      mpas_ocn_time_integration_rk4.o \
                      mpas_ocn_time_integration_split.o

clean:
	$(RM) *.o *.mod *.f90 libdycore.a

ifeq "$(USE_PARAVIEW)" "true"
PARAVIEW_FLAGS := -I../../paraview_catalyst
endif

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) $(PARAVIEW_FLAGS)
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) $(PARAVIEW_FLAGS)
endif
