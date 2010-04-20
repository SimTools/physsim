all: dir
	(cd pltppk;    	xmkmf; make)
	(cd pltptop;    xmkmf; make)
	(cd pltsga;    	xmkmf; make)
	(cd pltfbx;    	xmkmf; make)
	(cd pltfb1se;   xmkmf; make)
	(cd pltepk;   	xmkmf; make)

.PHONY.: dir
dir:
	(cd pltppk;     mkdir -p prod)	
	(cd pltptop;    mkdir -p prod)
	(cd pltsga;    	mkdir -p prod)
	(cd pltfbx;    	mkdir -p prod)
	(cd pltfb1se;   mkdir -p prod)
	(cd pltepk;   	mkdir -p prod)
clean:
	(cd pltppk;    	make clean)
	(cd pltptop;    make clean)
	(cd pltsga;    	make clean)
	(cd pltfbx;    	make clean)
	(cd pltfb1se;   make clean)
	(cd pltepk;   	make clean)
	rm -f *.o; rm -f *~ 

cleanall: clean
	(cd pltppk;    	make cleanall)
	(cd pltptop;    make cleanall)
	(cd pltsga;    	make cleanall)
	(cd pltfbx;    	make cleanall)
	(cd pltfb1se;   make cleanall)
	(cd pltepk;   	make cleanall)
	rm -f *.o; rm -f *.a; rm -f *.exe; rm -f *~ 
