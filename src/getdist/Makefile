PROGS = getdist

FFLAGS = -Vaxlib -W0 -WB -vec_report0
BINDIR = ../../bin

LIBS = $(LAPACKLIB)
F90FLAGS = -DMATRIX_SINGLE $(FFLAGS)

DISTFILES = inifile.o utils.o ParamNames.o Matrix_utils.o settings.o GetDist.o

all : getdist

.f.o:
	f77 $(F90FLAGS) -c $<

%.o: %.f90
	$(FC) $(F90FLAGS) -c $*.f90

%.o: %.F90
	$(FC) $(F90FLAGS) -c $*.F90

clean:
	rm -f *.o *.mod *.d *.pc *.obj $(BINDIR)/getdist

getdist: $(DISTFILES)
	$(FC) -o $(BINDIR)/getdist $(DISTFILES) $(F90FLAGS) $(LIBS)
