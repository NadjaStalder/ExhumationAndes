.SUFFIXES:.out .o .s .c .F .f .f90 .e .r .y .yr .ye .l .p .sh .csh .h

include Makefile.build
include Makefile.macbook

OBJECTS = \
module_definitions.o \
initialize_parameters.o 

.f90.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f90

.f.o:
	$(F90) $(FLAGS) $(INCLUDE) $*.f

.c.o:
	$(CC) $(FLAGS) $(INCLUDE) $*.c

SITE:   $(OBJECTS)
	$(F90) $(OPTIONS) $(OBJECTS) $(INCLUDE) $(LIBS) -o ../creategrid

clean:
	rm *.o *.mod

