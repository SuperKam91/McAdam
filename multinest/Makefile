PROGS = libnest3.a 

FFLAGS =  $(EXTRAFLAGS) $(FDEBUGFLAGS)

  
LINKLIB = ld -shared  
LIBDIR = ../lib
 
NSOBJECTS = utils.o utils1.o priors.o kmeans_clstr.o xmeans_clstr.o posterior.o nested.o

%.o: %.f90
	$(FC) $(FFLAGS) -c -o $@ $^ 

%.o: %.F90
	$(FC) $(FFLAGS) -c -o $@ $^ 

 
all: $(PROGS)
 
libnest3.so: $(NSOBJECTS) 
	$(LINKLIB) -o $(LIBS) $@ $^ 
 
libnest3.a: $(NSOBJECTS) 
	$(AR) $@ $^ 
	-mv $@ $(LIBDIR) 

clean: 
	-rm $(LIBDIR)/libnest3.*  *.o *.mod

