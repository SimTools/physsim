##############################################################################
## (Description)
##    Make file for precision EW tests.
##	run      : to run.
##	clean    : to clean the directory.
##	cleanall : to remove everything reproducible.
## (Update Record) 
##    94/11/25  K.Fujii 	Original.
##
##
##############################################################################
#
# Define the name of the executable. 
#
PROGRAM       = test

#
# Directory to store executables.
#
GODIR	      = .

#
# Routines to complie.
#
OBJS	      = cupltst.o \
		gamffs.o \
		gamffv.o \
		gamsff.o

#
# Libraries to link.
#
LIBS	    = ../../glib.a
HLSV1DIR    = $(HOME)/work/helas_v1
HLSV1LIB    = $(HLSV1DIR)/lib/glib.a
TOPAZLIB    = /proj/phys/TOPAZ/94a/lib
HBOOK_LIB   = /kek/cern/94a/lib
# TPPHYLIB      = /home/fujiik/work/korale_exp/tp_physics/glib.a
ANALLIB     = -Wl,-L$(TOPAZLIB) -lPRSTER  -lGEN \
	      -lHELASV204 -lBASES25 -lSIMJLC -lFORTLIB  -lCOM -lCLI \
	      -lMANAGER -lTBS -lFULL -lPRSTER -lTBS -lMANAGER \
	      -lHPKTOHBK -lFORTLIB -lCOM -lCLI \
	      -Wl,-L$(HBOOK_LIB) -lpacklib -lvec -lU77

#
# Specify compiler and linker.
#
F77	      = f77
LD	      = $(F77)

#
# Build, compiler, and linker options.
#
BLDOPT 	      = "R(UNIX,KEK,JLC)"
FFLAGS	      = -K -O +ppu
LDFLAGS	      = $(FFLAGS)

#
# Build and Compile the source codes.
#
.f.o :
	build $*.f $(BLDOPT)
	mv SYSPRINT SYSPRINT.$*.f
	$(F77) SYSPRINT.$*.f -c $(FFLAGS) -o $*.o 
	rm SYSPRINT.$*.f

#
# Make options.
#

all:		$(PROGRAM)

$(PROGRAM):     $(OBJS)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) $(ANALLIB) $(HLSV1LIB)\
		-o $(GODIR)/$(PROGRAM).exe
		@echo "done"

clean:;		@rm -f $(OBJS) core

cleanall:;	@rm -f $(OBJS) $(GODIR)/*.exe core
