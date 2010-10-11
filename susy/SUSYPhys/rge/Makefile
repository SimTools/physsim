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
	(cd alfs; make cleanall)
	(cd gut; make cleanall)
	(cd sfm; make cleanall)
	rm -f *.o; rm -f *.a; rm -f *.exe; rm -f *~ 
