PROGRAM1 = pltgrq
PROGRAM2 = pltsg0
OBJS     = sgtt0.o gmthet.o gqtabl.o

FC = gfortran -m32
FFLAGS = -fno-automatic

all: $(PROGRAM1) $(PROGRAM2)

$(PROGRAM1): $(PROGRAM1).F $(OBJS)
	$(FC) $(PROGRAM1).F -o $(PROGRAM1) $(OBJS) -L.. -lGENLIBTT -L../../TTPhys/tt_thresh/lib/ -lPOTQCD -lPHYLIBTT

$(PROGRAM2): $(PROGRAM2).F $(OBJS)
	$(FC) $(PROGRAM2).F -o $(PROGRAM2) $(OBJS) -L.. -lGENLIBTT -L../../TTPhys/tt_thresh/lib/ -lPOTQCD -lPHYLIBTT

.SUFFIXES: .F
.F.o:
	$(FC) -c $(FFLAGS) $<

clean:
	rm -f *.o
distclean:
	rm -f $(PROGRAM1) $(PROGRAM2) fort.*
