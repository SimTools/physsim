all:
	(cd tt_thresh;      make)
	(cd tt_thresh_new;  make)

clean:
	(cd tt_thresh;      make clean)
	(cd tt_thresh_new;  make clean)
	rm -f *.o; rm -f *~ 

distclean:
	(cd tt_thresh;      make cleanall)
	(cd tt_thresh_new;  make cleanall)
	rm -f *.o; rm -f *.a; rm -f *.exe; rm -f *~ 

tar:
	(cd tt_thresh;      make cleanall)
	(cd tt_thresh_new;  make cleanall)
	rm -f *.o; rm -f *.a; rm -f *.exe; rm -f *~ 
	(cd ../;   tar cvf   TTPhys.tar TTPhys)
	(cd ../;   compress  TTPhys.tar)
	(cd ../;   chmod -R 700 *)

