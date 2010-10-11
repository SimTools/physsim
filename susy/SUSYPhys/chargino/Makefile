all:
	(cd lib;    xmkmf; make)
	(cd plot;   xmkmf; make)

clean:
	(cd lib;      rm -f *.o *.bak *~)
	(cd plot;     make clean)
	rm -f *.o; rm -f *~ 

cleanall:
	(cd lib;      make clean distclean)
	(cd plot;     make cleanall distclean)
	rm -f *.o; rm -f *.a; rm -f *.exe; rm -f *~
