# Object files:
# find all files ending in 90, and change .f90 and .F90 to .o, 
OBJECTS = $(patsubst %.F90,%.o,$(patsubst %.f90,%.o,$(wildcard *90)))

LIBDIR = ../lib

# include the dependencies
include makefile.dep

FFLAGS = $(EXTRAFLAGS) $(FDEBUGFLAGS) 

# "make" builds all
all: libchord.a

# Rule for building polychord library
libchord.a: $(OBJECTS)
	$(AR) $@ $^
	-cp $@ $(LIBDIR) 

# General rule for building object file (.o)  from fortran files (.f90/.F90)
%.o: %.f90
	$(FC) $(FFLAGS) -c $<
%.o: %.F90
	$(FC) $(FFLAGS) -c $<

# Utility targets
.PHONY: clean

clean:
	rm -f *.o *.mod *.MOD libchord.a $(LIBDIR)/libchord.a
