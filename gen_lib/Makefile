all: Makefiles
	(cd bm;  make install)
	(cd kn;  make install)
	(cd hg;  make install)
	(cd sm;  make install)
	(cd ss;  make install)

Makefiles:
	(cd bm;  xmkmf)
	(cd kn;  xmkmf)
	(cd hg;  xmkmf)
	(cd sm;  xmkmf)
	(cd ss;  xmkmf)

clean:
	(cd bm;  make clean)
	(cd kn;  make clean)
	(cd hg;  make clean)
	(cd sm;  make clean)
	(cd ss;  make clean)

distclean:: clean
	(cd bm;  make distclean)
	(cd kn;  make distclean)
	(cd hg;  make distclean)
	(cd sm;  make distclean)
	(cd ss;  make distclean)
