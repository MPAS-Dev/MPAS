.SUFFIXES: .F .c .o

FILTER_OBJS = mpas_filter_test_compute_interval.o

all: filters

filters: mpas_filter_driver.o $(FILTER_OBJS)
	ar -ru libfilters.a mpas_filter_driver.o $(FILTER_OBJS)

clean:
	$(RM) *.o *.mod *.f90 libfilters.a
	@# Certain systems with intel compilers generate *.i files
	@# This removes them during the clean process
	$(RM) *.i

.F.o:
	$(RM) $@ $*.mod
ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../framework -I../operators
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) -I../framework -I../operators
endif

.c.o:
	$(CC) $(CFLAGS) $(CPPFLAGS) $(CPPINCLUDES) -I../framework -I../operators -c $<
