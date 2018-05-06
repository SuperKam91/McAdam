PROGS = genCMBCM 
 
FFLAGS = $(EXTRAFLAGS) $(FDEBUGFLAGS) 
CFLAGS = $(EXTRAFLAGS) $(CDEBUGFLAGS) 
CXXFLAGS = $(EXTRAFLAGS) $(CDEBUGFLAGS)
 
AR = ar r  
LINKLIB = ld -shared 
BINDIR = ../../bin
 
OBJECTS = params.o telescope.o cbessi0.o utilities.o fits_handling.o globals.o confusion_noise.o inputs.o numerics.o integrals.o cmbcm.o
 
%.o: %.f90 
	$(FC) $(FFLAGS) -c -o $@ $^ 
	
%.o: %.F90 
	$(FC) $(FFLAGS) -c -o $@ $^ 

%.o: %.c 
	$(CC) $(CFLAGS) -c -o $@ $^ 
 
%.o: %.f 
	$(FC) $(FFLAGS) -c -o $@ $^ 
 
 
all:  genCMBCM 
 
genCMBCM: $(OBJECTS)
	$(FC) $(FFLAGS) -o genCMBCM $(OBJECTS) $(LAPACKLIB) -lcfitsio
	-mv $@ $(BINDIR) 
 
clean: 
	-rm *.o *.mod $(BINDIR)/genCMBCM
