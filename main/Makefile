include ../GNUmakefile.chssi

SIP_DIR=../sip

FRAMELIB=../lib/libframelib.a
ERD_LIB=../lib/liberd.a
SIP_LIB=../lib/libsip.a
SIAL_LIB=../lib/libsial.a
SIP_SHARED_LIB=../lib/libsip_shared.a
ACES2_LIB=../lib/libaces2.a
GEOPT_LIB=../lib/libgeopt.a
SYMCOR_LIB=../lib/libsymcor.a
OED_LIB=../lib/liboed.a
DUP_LIB=../lib/libdup.a
INTEGRAL_OBJS=$(shell ls ../direct_integral/*.o 2>/dev/null)
WORKER_OBJS=$(shell ls ../worker/*.o 2>/dev/null)
MANAGER_OBJS=$(shell ls ../manager/*.o 2>/dev/null)

ACES3_EXE=../bin/xaces3

FSOURCE=$(wildcard *.F)
F_OBJS=$(FSOURCE:.F=.o)
fSOURCE=$(wildcard *.f)
f_OBJS:=$(fSOURCE:.f=.o)
cSOURCE=$(wildcard *.c)
c_OBJS:=$(cSOURCE:.c=.o)
CPP_SOURCE=$(wildcard *.cpp)
CPP_OBJS:=$(CPP_SOURCE:.cpp=.o)

MAIN_OBJS=$(F_OBJS) $(f_OBJS) $(c_OBJS) $(CPP_OBJS)

all: tabula_rasa $(ACES3_EXE)

tabula_rasa: ;
	rm -f make.out make.err

$(ACES3_EXE): xaces3
	cp xaces3 $(ACES3_EXE)

xaces3: $(FRAMELIB) $(ERD_LIB) $(SIP_LIB) $(SIAL_LIB) $(SIP_SHARED_LIB) \
        $(ACES2_LIB) $(GEOPT_LIB) $(SYMCOR_LIB) $(OED_LIB) $(DUP_LIB) $(INTEGRAL_OBJS) \
        $(MAIN_OBJS) $(WORKER_OBJS) $(MANAGER_OBJS)
	$(FC) $(FFLAGS) $(LFLAGS) -o xaces3 \
	       $(MAIN_OBJS) $(MANAGER_OBJS) $(WORKER_OBJS) $(INTEGRAL_OBJS) \
	       $(LIB_DIRS) $(LIBS) && test -x xaces3  || rm -f xaces3

binclean: ;
	rm -f xaces3

ppclean clean distclean: % : tabula_rasa binclean
	rm -f *.o

rebuild relink: % : binclean $(ACES3_EXE)

libclean archive: % : ;

.f.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.F.o:
	$(FC) -c $(FFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

.cpp.o: 
	$(CPP) -c $(CPPFLAGS) $(INCLUDE_DIRS) $< 1>>make.out 2>>make.err

