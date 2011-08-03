FC=gfortran
CERN_ROOT=/nfs/g/ilc/soft/utils64/cernlib_64/2005

FFLAGS=-fPIC -I$(KFLIBROOT)/gen_lib/sm/include 
LIBSRC=isr_function.o isr_remnant.o whizard_spectrum.o \
	isr_simple_spectrum.o isr_whizard_spectrum.o 

LIB=libwhizfunc.a

OBJS=testisr.o isr_function.o
OBJS2=testisr2.o isr_function.o
CXXFLAGS=$(shell root-config --cflags)
LDFLAGS=$(shell root-config --libs) -lgfortran

OBJGEN=genisr.o isr_remnant.o isr_function.o
OBJPLOT=funcplot.o isr_function.o isr_remnant.o

# all: lib funcplot genisr test test2
all: lib

lib: $(LIBSRC)
	ar crv $(LIB) $(LIBSRC)

funcplot:$(OBJPLOT)
	$(CXX) -o funcplot $(OBJPLOT) $(LDFLAGS)

genisr: $(OBJGEN)
	$(FC) -o genisr $(OBJGEN) 

test: $(OBJS)
	$(FC) -o test1 $(OBJS) -L$(WHIZDIR)/lib -lwhizard -lpdfdummy -lproc \
	$(WHIZDIR)/lib/mpi90.o $(WHIZDIR)/lib/circe.o \
	$(WHIZDIR)/lib/circe2.o \
	-L$(WHIZDIR)/a6f/lib -linclude \
	$(TAUOLADIR)/libtauola.a $(PHOTOSDIR)/libphotos.a \
	-L$(STDHEPDIR)/lib -lstdhep -lFmcfio \
	-L$(CERN_ROOT)/lib -lpythia6422 -lpdflib804 -lmathlib -lkernlib

test2: $(OBJS2)
	$(FC) -o test2 $(OBJS2) -L$(WHIZDIR)/lib -lwhizard -lpdfdummy -lproc \
	$(WHIZDIR)/lib/mpi90.o $(WHIZDIR)/lib/circe.o \
	$(WHIZDIR)/lib/circe2.o \
	-L$(WHIZDIR)/a6f/lib -linclude \
	$(TAUOLADIR)/libtauola.a $(PHOTOSDIR)/libphotos.a \
	-L$(STDHEPDIR)/lib -lstdhep -lFmcfio \
	-L$(CERN_ROOT)/lib -lpythia6422 -lpdflib804 -lmathlib -lkernlib

clean: 
	rm -f *.o *~

cleanall: clean
	rm -f test1 test2 genisr funcplot a.out *.root $(LIB)

funcplot.o: funcplot.cxx
	g++ -c funcplot.cxx $(CXXFLAGS)

genisr.o: genisr.f90
	gfortran -c -I$(WHIZDIR)/whizard-src genisr.f90

testisr.o: testisr.f90
	gfortran -c -I$(WHIZDIR)/whizard-src testisr.f90
testisr2.o: testisr2.f90
	gfortran -c -I$(WHIZDIR)/whizard-src testisr2.f90

whizard_spectrum.o: whizard_spectrum.f90
	gfortran -c -fPIC -I$(WHIZDIR)/whizard-src whizard_spectrum.f90

