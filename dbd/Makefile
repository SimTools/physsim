
SUBDIR = TTHStudy TTZStudy TTBBStudy TTStudy

all:
	( cd whizfunc && make )
	for i in $(SUBDIR) ; do \
	( cd $$i && xmkmf -a && make ) ; \
	done

