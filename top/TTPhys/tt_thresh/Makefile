all:
	(cd lib;       xmkmf; make; make install)
	(cd plots;     xmkmf; make)
	(cd contours;  xmkmf; make)

clean:
	(cd lib;       make clean)
	(cd plots;     make clean)
	(cd contours;  make clean)

cleanall:
	(cd lib;       make distclean)
	(cd plots;     make cleanall)
	(cd contours;  make cleanall)

tar:
	(cd lib;       make clean)
	(cd plots;     make cleanall)
	(cd contours;  make cleanall)
	rm -f *.o; rm -f *.a; rm -f *.exe; rm -f *~ 
	(cd ../;   tar cvf   tt_thresh.tar tt_thresh)
	(cd ../;   compress  tt_thresh.tar)
	(cd ../;   chmod -R 700 *)
