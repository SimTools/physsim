########################################################################	
##
## --------------------------------------------------
##   Makefile for K.Fujii's physics study libraries.
## --------------------------------------------------
##
## (Version Information)
##   Version: 98a
##   Release: 1
##
## (Update Record)
##    1999/05/24  K.Fujii	A temporary JSF version (98a-1)
##				with C++ wrappers to be replaced
##				in future by a fully C++ version.
##
## (Description)
##   In order to use this package you should first set some
##   environmental variables:
##
##	$ cd <the directory containing this file>
##	$ export KFLIBROOT=`pwd`
##	Unix.*.Root.DynamicPath: .:$(ROOTSYS)/lib:$(JSFROOT)/lib
##      Unix.*.Root.MacroPath:   .:$(ROOTSYS)/macros:$(JSFROOT)/example/guiexam1
##
##   where the JSF standard environmental variables such as
##   ROOTSYS, JSFROOT, LCLIBROOT, CERN_ROOT, .... are assumed
##   to be set beforehand (see HowToInstall of JSF for details).
##   Once the environment is set up, do:
##
##	$ make
##
##   This will create libraries necessary to build various MC
##   event generators based on BASES/SPRING stored in
##	higgs	: zh, ...
##	susy	: sfsf, xcxc, ...
##	top	: eett, nntt, tth, tt, ttz, ...
##	twoph	: eeff (f=e,mu,tau,q)
##	wz	: eeww, eez, enw, nnww, nnz, ww, wwz, zz, ...
##   These subdirectories contain subsubdirectories named
##   <PROCESSNAME>Study. The top levelMakefile does NOT build 
##   individual Monte Carlo event generators: it only creates
##   necessary libraries for them.
##   You should "cd" to one of these subsubdirectories and do
##
##    	$ xmkmf -a
##      $ make
##   
##   For more details, see README.
##
## (Targets)
##	all       	: creates libraries in lib.
##	Makefiles 	: creates Makefiles.
##	clean     	: deletes *.o ...
##	distclean 	: deletes even lib/*.
## 
########################################################################	

SHELL	= /bin/bash
MFLAGS	=
CURRDIR	= .

ifeq ($(strip $(LCBASEDIR)),)
SUBDIR1	= dgen_lib gen_lib
else
SUBDIR1	= dgen_lib gen_lib
endif
SUBDIR2	= anl_lib top/TTLib top/TTPhys/tt_thresh/lib \
	  top/TTPhys/tt_thresh_new/lib Genlib Hellib
SUBDIRS = $(SUBDIR1) $(SUBDIR2)

all:
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIR1); do \
	(cd $$i; echo ``making'' all ``in $(CURRDIR)/$$i...''; \
	$(MAKE) $(MFLAGS) all); \
	done
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIR2); do \
	(cd $$i; echo ``making'' all ``in $(CURRDIR)/$$i...''; \
	xmkmf -a; \
	$(MAKE) $(MFLAGS) install); \
	done
	
Makefiles:
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIR1); do \
	(cd $$i; echo ``making'' Makefiles ``in $(CURRDIR)/$$i...''; \
	$(MAKE) $(MFLAGS) Makefiles); \
	done
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIR2); do \
	(cd $$i; echo ``making'' Makefiles ``in $(CURRDIR)/$$i...''; \
	xmkmf); \
	done
	
clean:
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS); do \
	(cd $$i; echo ``making'' clean ``in $(CURRDIR)/$$i...''; \
	$(MAKE) $(MFLAGS) clean); \
	done

distclean:
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS); do \
	(cd $$i; echo ``making'' distclean ``in $(CURRDIR)/$$i...''; \
	$(MAKE) $(MFLAGS) distclean); \
	done
	rm -f lib/*.a lib/*.so lib/*.so.*
