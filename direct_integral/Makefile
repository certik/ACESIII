include ../GNUmakefile.chssi

DEFINES= 

FSOURCE=$(wildcard *.F)
F_OBJS=$(FSOURCE:.F=.o)
fSOURCE=$(wildcard *.f)
f_OBJS:=$(fSOURCE:.f=.o)
cSOURCE=$(wildcard *.c)
c_OBJS:=$(cSOURCE:.c=.o)

all rebuild relink: % : tabula_rasa $(F_OBJS) $(f_OBJS) $(c_OBJS)

tabula_rasa: ;
	rm -f make.out make.err

ppclean clean distclean: % : tabula_rasa
	rm -f *.o

libclean binclean archive: % : ;

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.F.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

