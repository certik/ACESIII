include ../GNUmakefile.chssi

FSOURCE=$(wildcard *.F)
F_OBJS=$(FSOURCE:.F=.o)
fSOURCE=$(wildcard *.f)
f_OBJS:=$(fSOURCE:.f=.o)
cSOURCE=$(wildcard *.c)
c_OBJS:=$(cSOURCE:.c=.o)
CPP_SOURCE=$(wildcard *.cpp)
CPP_OBJS:=$(CPP_SOURCE:.cpp=.o)

OEDLIB=../lib/liboed.a

all: tabula_rasa $(OEDLIB)

tabula_rasa: ;
	rm -f make.out make.err

$(OEDLIB): liboedlib.a
	cp liboedlib.a $(OEDLIB)

liboedlib.a: $(F_OBJS) $(f_OBJS) $(c_OBJS)
	ar $(ARFLAGS) liboedlib.a $(F_OBJS) $(f_OBJS) $(c_OBJS)

libclean: ;
	rm -f liboedlib.a

ppclean clean distclean: % : tabula_rasa libclean
	rm -f *.o

rebuild: libclean $(FRAMELIB)

binclean relink archive: % : ;

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.F.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

