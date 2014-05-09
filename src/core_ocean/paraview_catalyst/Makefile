.SUFFIXES: .F .o

OBJS = mpas_paraview_catalyst.o

all: paraview

paraview: $(OBJS)
	( cd Adaptor-build; cmake -DParaView_DIR="$(PARAVIEW)/lib/cmake/paraview-4.0" ../Adaptor; make)
	ar -ru libparaview.a $(OBJS)

mpas_paraview_catalyst.o: 

clean:
	$(RM) *.o *.mod *.f90 libparaview.a
	( cd Adaptor-build; make clean)

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) -I../../inc $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../../framework -I../../external/esmf_time_f90
endif

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) $(CPPINCLUDES) -c $<
