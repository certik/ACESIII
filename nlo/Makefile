include ../GNUmakefile.chssi

FSOURCE=$(wildcard *.F)
F_OBJS=$(FSOURCE:.F=.o)
fSOURCE=$(wildcard *.f)
f_OBJS:=$(fSOURCE:.f=.o)
cSOURCE=$(wildcard *.c)
c_OBJS:=$(cSOURCE:.c=.o)
CPP_SOURCE=$(wildcard *.cpp)
CPP_OBJS:=$(CPP_SOURCE:.cpp=.o)

NLOLIB=../lib/libnlo.a

all: tabula_rasa $(NLOLIB)

tabula_rasa: ;
	rm -f make.out make.err

$(NLOLIB): libnlo.a
	cp libnlo.a $(NLOLIB)

libnlo.a: $(F_OBJS) $(f_OBJS) $(c_OBJS)
	ar $(ARFLAGS) libnlo.a $(F_OBJS) $(f_OBJS) $(c_OBJS)

libclean: ;
	rm -f libnlo.a

ppclean clean distclean: % : tabula_rasa libclean
	rm -f *.o

rebuild: libclean $(NLOLIB)

binclean relink archive: % : ;

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.F.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

