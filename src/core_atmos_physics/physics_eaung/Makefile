.SUFFIXES: .F .o

dummy:
	@echo "****** compile physics_eaung ******"

OBJS = \
	kinds.o               \
	module_eaung_params.o \
	module_eaung_util.o

all: physics_eaung

physics_eaung: $(OBJS)
	ar -ru ./../libphys.a $(OBJS)

# DEPENDENCIES:
module_eaung_params.o:    \
	kinds.o

module_eaung_util.o:      \
	kinds.o               \
	module_eaung_params.o

clean:
	$(RM) *.f90 *.o *.mod

.F.o:
	$(RM) $@ $*.mod
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES) -I../../framework -I../../operators
#	$(RM) $*.f90
