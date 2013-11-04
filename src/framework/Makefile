.SUFFIXES: .F .o

DEPS := $(wildcard ../core_$(CORE)/Registry.xml)

OBJS = mpas_kind_types.o \
       mpas_framework.o \
       mpas_timer.o \
       mpas_timekeeping.o \
       mpas_configure.o \
       mpas_packages.o \
       mpas_constants.o \
       mpas_dmpar_types.o \
       mpas_attlist.o \
       mpas_grid_types.o \
       mpas_hash.o \
       mpas_sort.o \
       mpas_block_decomp.o \
       mpas_block_creator.o \
       mpas_dmpar.o \
       mpas_io.o \
       mpas_io_streams.o \
       mpas_io_input.o \
       mpas_io_output.o \
       mpas_io_units.o \
       streams.o

all: framework $(DEPS)

framework: $(OBJS)
	ar -ru libframework.a $(OBJS)

mpas_configure.o: mpas_dmpar.o mpas_io_units.o $(DEPS)

mpas_packages.o: $(DEPS)

mpas_framework.o: mpas_dmpar.o \
                  mpas_io_input.o \
                  mpas_io_output.o \
                  mpas_io.o \
                  mpas_grid_types.o \
                  mpas_configure.o \
                  mpas_timer.o \
                  mpas_sort.o \
                  mpas_io_units.o \
                  mpas_packages.o

mpas_constants.o: mpas_kind_types.o mpas_io_units.o

mpas_dmpar_types.o : mpas_kind_types.o mpas_io_units.o

mpas_attlist.o: mpas_kind_types.o mpas_io_units.o

mpas_grid_types.o: mpas_kind_types.o mpas_dmpar_types.o mpas_attlist.o mpas_io_units.o mpas_packages.o $(DEPS)

mpas_dmpar.o: mpas_sort.o streams.o mpas_kind_types.o mpas_grid_types.o mpas_hash.o mpas_io_units.o

mpas_sort.o: mpas_kind_types.o mpas_io_units.o

mpas_timekeeping.o: mpas_kind_types.o mpas_io_units.o

mpas_timer.o: mpas_kind_types.o mpas_io_units.o

mpas_block_decomp.o: mpas_grid_types.o mpas_hash.o mpas_configure.o mpas_io_units.o

mpas_block_creator.o: mpas_dmpar.o mpas_hash.o mpas_sort.o mpas_configure.o mpas_io_units.o $(DEPS)

mpas_io.o: mpas_dmpar_types.o mpas_io_units.o

mpas_io_streams.o: mpas_attlist.o mpas_grid_types.o mpas_timekeeping.o mpas_io.o mpas_io_units.o $(DEPS)

mpas_io_input.o: mpas_grid_types.o mpas_dmpar.o mpas_block_decomp.o mpas_block_creator.o mpas_sort.o mpas_configure.o mpas_timekeeping.o mpas_io_streams.o mpas_io_units.o $(DEPS)

mpas_io_output.o: mpas_grid_types.o mpas_dmpar.o mpas_sort.o mpas_configure.o mpas_io_streams.o mpas_io_units.o $(DEPS)

mpas_io_units.o:

clean:
	$(RM) *.o *.mod *.f90 libframework.a

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) -I../inc $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../external/esmf_time_f90
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../inc -I../external/esmf_time_f90
endif

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) $(CPPINCLUDES) -c $<
