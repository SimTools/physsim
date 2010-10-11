all:
	(cd lib;    xmkmf; make)
	(cd plot;   xmkmf; make)

clean:
	(cd lib;      rm -f *.o *.bak *~)
	(cd plot;     make clean)

cleanall:
	(cd lib;      make clean)
	(cd plot;     make cleanall distclean)

