all:
	(cd lib;   xmkmf; make; make install)
	(cd plot;  make)

clean:
	(cd lib;      make clean)
	(cd plot;     make clean)

cleanall: clean
	(cd lib;      make distclean)
	(cd plot;     make cleanall)

tar:
	(cd lib;      make clean)
	(cd plot;     make cleanall)
	rm -f *.o; rm -f *.a; rm -f *.exe; rm -f *~ 
	(cd ../;   tar cvf   tt_thresh_new.tar tt_thresh_new)
	(cd ../;   compress  tt_thresh_new.tar)
	(cd ../;   chmod -R 700 *)
