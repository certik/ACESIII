include ../GNUmakefile.chssi

SIAL_INCLUDES=$(INCLUDE_DIRS) -I../sial
SIAL_LIB=../lib/libsial.a
SIAL_COMPILER=../bin/sial
SIP_SHARED_LIB=../lib/libsip_shared.a

FSOURCE=$(wildcard *.F)
F_OBJS=$(FSOURCE:.F=.o)
fSOURCE=$(wildcard *.f)
f_OBJS:=$(fSOURCE:.f=.o)
cSOURCE=$(wildcard *.c)
c_OBJS:=$(cSOURCE:.c=.o)
CPP_SOURCE=$(wildcard *.cpp)
CPP_OBJS:=$(CPP_SOURCE:.cpp=.o)

all: tabula_rasa $(SIAL_COMPILER)

tabula_rasa: ;
	rm -f make.out make.err

$(SIAL_COMPILER): sial
	cp sial $(SIAL_COMPILER)

sial: $(F_OBJS) $(f_OBJS) $(c_OBJS) $(CPP_OBJS) $(SIAL_LIB) $(SIP_SHARED_LIB) 
	$(SERIAL_CPP) $(CPPFLAGS) $(SIAL_INCLUDES) $(LIB_DIRS) -o sial \
	$(F_OBJS) $(f_OBJS) $(c_OBJS) $(CPP_OBJS) $(SIAL_COMPILER_LFLAGS) \
	$(SIAL_COMPILER_LIBS) \
	&& test -x sial || rm -f sial

binclean: ;
	rm -f sial

ppclean clean distclean: % : tabula_rasa binclean
	rm -f *.o

rebuild relink: % : binclean $(SIAL_COMPILER)

libclean archive: % : ;

.f.o:
	$(FC) -c $(FFLAGS) $(SIAL_INCLUDES) $< 1>>make.out 2>>make.err

.F.o:
	$(FC) -c $(FFLAGS) $(SIAL_INCLUDES) $< 1>>make.out 2>>make.err

.c.o:
	$(CC) -c $(CFLAGS) $(SIAL_INCLUDES) $< 1>>make.out 2>>make.err

.cpp.o: 
	$(CPP) -c $(CPPFLAGS) $(SIAL_INCLUDES) $< 1>>make.out 2>>make.err

