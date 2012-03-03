include ../GNUmakefile.chssi

ACES2_LIB=../lib/libaces2.a
EXE=../bin/xtest_compare

FSOURCE=$(wildcard *.F)
F_OBJS=$(FSOURCE:.F=.o)
fSOURCE=$(wildcard *.f)
f_OBJS:=$(fSOURCE:.f=.o)
cSOURCE=$(wildcard *.c)
c_OBJS:=$(cSOURCE:.c=.o)
CPP_SOURCE=$(wildcard *.cpp)
CPP_OBJS:=$(CPP_SOURCE:.cpp=.o)

MAIN_OBJS=$(F_OBJS) $(f_OBJS) $(c_OBJS) $(CPP_OBJS)

all: tabula_rasa $(EXE)

tabula_rasa: ;
	rm -f make.out make.err

$(EXE): xtest_compare
	cp xtest_compare $(EXE)

xtest_compare: $(ACES2_LIB) $(MAIN_OBJS)
	$(FC) $(FFLAGS) $(LFLAGS) -o xtest_compare \
	       $(MAIN_OBJS) \
	       $(LIB_DIRS) $(LIBS) && test -x xtest_compare  || rm -f xtest_compare

binclean: ;
	rm -f xtest_compare

ppclean clean distclean: % : tabula_rasa binclean
	rm -f *.o

rebuild relink: % : binclean $(EXE)

libclean archive: % : ;

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.F.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.cpp.o: 
	$(CPP) -c $(CPPFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

