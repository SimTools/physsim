all:
	(cd alfs; xmkmf; make)
	(cd gut;  xmkmf; make)
	(cd sfm;  xmkmf; make)

clean:
	(cd alfs; make clean)
	(cd gut; make clean)
	(cd sfm; make clean)
	rm -f *.o; rm -f *~ 

cleanall:
	(cd alfs; make cleanall distclean)
	(cd gut; make cleanall distclean)
	(cd sfm; make cleanall distclean)
	rm -f *.o; rm -f *.a; rm -f *.exe; rm -f *~ 
