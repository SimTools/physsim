all:
	(cd lib; xmkmf; make)
	(cd chargino; make)
	(cd neutralino; make)
	(cd slepton; make)
	(cd rge; make)

clean:
	(cd lib; rm -f *.o *.bak *~)
	(cd chargino; make clean)
	(cd neutralino; make clean)
	(cd slepton; make clean)
	(cd rge; make clean)

cleanall:
	(cd lib;      make clean)
	(cd chargino; make cleanall)
	(cd neutralino; make cleanall)
	(cd slepton; make cleanall)
	(cd rge; make cleanall)
	rm -f *.o; rm -f *.a; rm -f *.exe; rm -f *~ 

tar:
	(cd lib;      make cleanall)
	(cd chargino; make cleanall)
	(cd neutralino; make cleanall)
	(cd slepton; make cleanall)
	(cd rge; make cleanall)
	rm -f *.o; rm -f *.a; rm -f *.exe; rm -f *~ 
	(cd ../;   tar cvf   susy_phys.tar susy_phys)
	(cd ../;   compress  susy_phys.tar)
	(cd ../;   chmod -R 700 *)
