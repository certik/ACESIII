# ACES Top-Level makefile. Version 1.6

# Note: At this time...
#
#              TOP-LEVEL GOALS CANNOT BE CHAINED TOGETHER!
#
#       I.E., `gmake clean all` will not do what you think. Instead,
#       try `gmake clean && gmake`.

# Note: "sm_" means submake while "tl_" means top-level. This
#       partitions the flags used by top-level make and its children.

################################################################################

# ARGUMENTS FOR TOP-LEVEL MAKE

# TREES is a list of directories that require top-level making.
TREES ?= $(wildcard *.src)

# ARGUMENTS FOR SUB-MAKE

# If this is a fresh build, include "install" in sm_goals so linking
# the binaries does not crash.
sm_goals := all install
sm_flags := # Do not set sm_flags here if other users are using this file,
            # e.g., at qtp. If this is a direct port, i.e., single-user
            # mode, then it is safe to put flags for sub-make here.
ifeq (${64BIT},1)
   sm_flags += 64BIT=1
endif

# MASTER MAKEFILES

tl_makefile_root ?= $(shell pwd)/Makefiles
tl_makefile      := ${tl_makefile_root}/GNUmakefile
tl_makefile_src  := ${tl_makefile_root}/GNUmakefile.src
tl_makefile_tl   := ${tl_makefile_root}/GNUmakefile.tl

# DIRECTORIES

# directories that do not get made
tl_tmp :=
tl_tmp += \
          lib bin include \
          module% \
          scr tmp work test configure \
          CVS Basis Info Makefiles Scripts Test%
tl_tmp += \
          molcas.% \
          %ource \
          %rchive backup% old \
          %.skip skip.%
tl_skip := $(strip ${tl_tmp})

# MRCC
#    The build must run in two steps: 1) compile all source without
# archiving then 2) create a monolithic libmrcclibs.a archive from
# the file sorted.objs made by
# > lorder [ace]*/*.o mb*/*.o mrcc/*.o | tsort > sorted.objs
# To this end, linking mrcc requires only two libraries:
# libmrcclibs.a and libmrcca3.a.
tl_mrcc_libs := apg apg2 cse eomprocs mbdirect mbgrad mbdiag mbtools mrcca3
tl_mrcc_bins := mrcc vibron
#tl_tmp := $(foreach tl_dir, ${tl_mrcc_libs}, \
#            $(shell test -d $(tl_dir) && echo $(tl_dir)) \
#           )
#tl_mrcc_libs := $(strip ${tl_tmp})
#tl_mrcc_objs := $(addsuffix /*.o, ${tl_mrcc_libs})
tl_lorder := lorder
tl_ar     := ar -rcs
ifeq (${64BIT},1)
   ifeq ($(shell uname -s),AIX)
      tl_lorder += -X64
      tl_ar     += -X64
   endif
endif

# libs and bins are (currently) built exactly the same way. This allows
# us to list certains bins in tl_libs after the particular
# libs they depend on (e.g., joda:symcor:ecp:libr -> "libr ecp symcor joda").
# Use this if one AME depends on another AME (e.g., vscf_ks:intgrt:libra3).
tl_tmp :=
tl_tmp += \
          mpistubs mpw p_mpw blas lapack eispack linpack \
          lb tools \
          acescore p_acescore acesutil libr ecp symcor joda \
          crust libr2 librt3 trp vcc
tl_tmp += \
          oed erd gamess p_gamess \
          sbtools sbcore libra3 intgrt oeptools \
          molcas util1 util2 util3 aixio parautil sewutil getmem_ma molpro \
          ${tl_mrcc_libs} \
          diagtools symtools omm imm
tl_tmp += \
          vmol vmol2ja vscf vtran intprc lambda dens anti bcktrn vdint \
          vcceh vee pdens
tl_libs := $(strip ${tl_tmp})

# Mature subtrees have their own build environment and should minimally
# depend on any other library (numerical libraries the exception). They
# most likely have their own makefile, but if they don't, we'll link ours in.
tl_trees := $(foreach tl_dir, ${TREES}, \
              $(shell test -d $(tl_dir) && echo $(tl_dir)) \
             )

# ANALYSIS

tl_dirs := $(foreach tl_dir, $(wildcard *), \
             $(shell test -d $(tl_dir) && echo $(tl_dir)) \
            )
tl_dirs := $(filter-out ${tl_skip} $(wildcard *.src), ${tl_dirs})
tl_libs := $(filter     ${tl_dirs},                   ${tl_libs})
tl_bins := $(filter-out ${tl_libs},                   ${tl_dirs})

# MRCC PASS 2 (linking)

tl_tmp := $(filter ${tl_mrcc_bins}, ${tl_dirs})
tl_mrcc_pass2 := $(addsuffix -p2, ${tl_tmp})

################################################################################

.PHONY: all
all: tl_goal=all
all: sm_goals_tmp=$(strip ${sm_goals})
all: sm_flags_tmp=$(strip ${sm_flags})
all: ${tl_libs} ${tl_bins} ${tl_trees} ${tl_mrcc_pass2}
	@echo
	@echo Successfully made all.

.PHONY: pponly
pponly: tl_goal=pponly
pponly: sm_goals_tmp=pponly
pponly: sm_flags_tmp=$(strip ${sm_flags})
pponly: ${tl_libs} ${tl_bins} $(filter-out chssi.src, ${tl_trees})
	@echo
	@echo Successfully made pponly.

.PHONY: clean
clean: tl_goal=clean
clean: sm_goals_tmp=clean
clean: sm_flags_tmp=$(strip ${sm_flags}) "FAST=0"
clean: ${tl_libs} ${tl_bins} ${tl_trees}
	@echo; rm -f libmrcclibs.a
	@echo Successfully made clean.

.PHONY: ppclean
ppclean: tl_goal=ppclean
ppclean: sm_goals_tmp=ppclean
ppclean: sm_flags_tmp=$(strip ${sm_flags}) "FAST=0"
ppclean: ${tl_libs} ${tl_bins} ${tl_trees}
	@echo; rm -f libmrcclibs.a
	@echo Successfully made ppclean.

.PHONY: libclean
libclean: tl_goal=libclean
libclean: sm_goals_tmp=libclean
libclean: sm_flags_tmp=$(strip ${sm_flags})
libclean: ${tl_libs} ${tl_bins} ${tl_trees}
	@echo; rm -f libmrcclibs.a
	@echo Successfully made libclean.

.PHONY: binclean
binclean: tl_goal=binclean
binclean: sm_goals_tmp=binclean
binclean: sm_flags_tmp=$(strip ${sm_flags})
binclean: ${tl_libs} ${tl_bins} ${tl_trees}
	@echo; rm -f libmrcclibs.a
	@echo Successfully made binclean.

.PHONY: distclean
distclean: tl_goal=distclean
distclean: sm_goals_tmp=distclean
distclean: sm_flags_tmp=$(strip ${sm_flags}) "FAST=0"
distclean: ${tl_libs} ${tl_bins} ${tl_trees}
	@echo; rm -f libmrcclibs.a
	@echo Successfully made distclean.

.PHONY: archive
archive: tl_goal=archive
archive: sm_goals_tmp=archive
archive: sm_flags_tmp=$(strip ${sm_flags})
archive: ${tl_libs} ${tl_bins} ${tl_trees}
	@if [ -d ./include ] ; then \
	    if [ -d ./archive ] ; then \
	       tar -z -c -f ./archive/include.tgz ./include ; \
	    else \
	       tar -z -c -f ./include.tgz ./include ; \
	    fi ; \
	 fi
	@echo
	@echo Successfully made archive.

# run regression tests
.PHONY: test
test: ;
	@if test -d test; \
	 then $(MAKE) -C test | tee test/gmake.out; \
	 else echo; echo "ERROR: missing test directory"; echo; exit 1; \
	 fi

# These are for handling GNUmakefile goals.
.PHONY: install backup debug
install:;@echo;echo There is no top-level \`make install\` at this time.;echo
backup: ;@echo;echo There is no top-level \`make backup\`  at this time.;echo
debug:  ;@echo;echo There is no top-level \`make debug\`   at this time.;echo

# top-level macros
.PHONY: rebuild
rebuild: tl_goal=rebuild
rebuild: sm_goals_tmp=libclean all install
rebuild: sm_flags_tmp=$(strip ${sm_flags})
rebuild: ${tl_libs} ${tl_bins} ${tl_trees} ${tl_mrcc_pass2}
	@echo
	@echo Successfully made rebuild \(${sm_goals_tmp}\).

.PHONY: relink
relink: tl_goal=relink
relink: sm_goals_tmp=binclean all install
relink: sm_flags_tmp=$(strip ${sm_flags})
relink: ${tl_libs} ${tl_bins} ${tl_trees} ${tl_mrcc_pass2}
	@echo
	@echo Successfully made relink \(${sm_goals_tmp}\).

################################################################################

# these are used if we make a directory from the tl command line
sm_goals_tmp=$(strip ${sm_goals})
sm_flags_tmp=$(strip ${sm_flags})

# STUPID TEST! This is the flag to test for a symlink. Either -h or -L.
hL:=-h

.PHONY: ${tl_libs} ${tl_bins} ${tl_trees} FORCE

FORCE: ;

${tl_mrcc_libs} ${tl_mrcc_bins} : nolib=1
${tl_libs} ${tl_bins} : FORCE
#	@if [ -f $@/GNUmakefile ] && [ ! ${hL} $@/GNUmakefile ] ; \
#	 then echo Moving $@/GNUmakefile to $@/GNUmakefile.$$$$.old ; \
#	      mv $@/GNUmakefile $@/GNUmakefile.$$$$.old ; \
#	 fi
#	@if [ -f $@/GNUmakefile.src ] && [ ! ${hL} $@/GNUmakefile.src ] ; \
#	 then echo Moving $@/GNUmakefile.src to $@/GNUmakefile.src.$$$$.old ; \
#	      mv $@/GNUmakefile.src $@/GNUmakefile.src.$$$$.old ; \
#	 fi
	@if [ ! -f $@/GNUmakefile ] && [ ! ${hL} $@/GNUmakefile ] && \
	    [ ! -f $@/Makefile    ] && [ ! ${hL} $@/Makefile    ] && \
	    [ ! -f $@/makefile    ] && [ ! ${hL} $@/makefile    ] ; \
	 then echo "SymLink $@/GNUmakefile -> ${tl_makefile}" ; \
	      ln -s ${tl_makefile} $@/GNUmakefile ; \
	 fi
	@if [ ! -f $@/GNUmakefile.src ] && [ ! ${hL} $@/GNUmakefile.src ] && \
	    [ ! -f $@/Makefile.src    ] && [ ! ${hL} $@/Makefile.src    ] ; \
	 then echo "SymLink $@/GNUmakefile.src -> ${tl_makefile_src}" ; \
	      ln -s ${tl_makefile_src} $@/GNUmakefile.src ; \
	 fi
	@test "x${nolib}" = "x1" -a \
	      $(patsubst %clean, clean, ${tl_goal}) != clean \
	  && $(MAKE) -C $@ ${sm_goals_tmp} ${sm_flags_tmp} BUILDLIB=0 \
	  || $(MAKE) -C $@ ${sm_goals_tmp} ${sm_flags_tmp}
	@if [ "${tl_goal}" = "distclean" ] ; \
	 then printf "\nRemoving links to makefiles in $@ . . . " ; \
	      for link in GNUmakefile GNUmakefile.src \
	                  Makefile    Makefile.src \
	                  makefile    makefile.src ; \
	      do test ${hL} $@/$$link && rm -f $@/$$link ; \
	      done ; \
	      printf "done!\n\n" ; \
	 fi

${tl_trees} : FORCE
#	@if [ -f $@/makefile ] && [ ! ${hL} $@/makefile ] ; \
#	 then echo Moving $@/makefile to $@/makefile.$$$$.old ; \
#	      mv $@/makefile $@/makefile.$$$$.old ; \
#	 fi
	@if [ ! -f $@/Makefile ] && [ ! ${hL} $@/Makefile ] && \
	    [ ! -f $@/makefile ] && [ ! ${hL} $@/makefile ] ; \
	    then echo "SymLink $@/makefile -> ${tl_makefile_tl}" ; \
	         ln -s ${tl_makefile_tl} $@/makefile ; \
	 fi
#	@if [ "${tl_goal}" = "all" ] ; \
#	 then cd $@ && $(MAKE) TREES=1 tl_makefile_root=${tl_makefile_root} ; \
#	 else cd $@ && $(MAKE) TREES=1 tl_makefile_root=${tl_makefile_root} \
#	                       ${sm_goals_tmp} ; \
#	 fi
ifeq (${64BIT},1)
	@cd $@ && $(MAKE) ${tl_goal} 64BIT=1 \
	                  tl_makefile_root=${tl_makefile_root}
else
	@cd $@ && $(MAKE) ${tl_goal} \
	                  tl_makefile_root=${tl_makefile_root}
endif
	@if [ "${tl_goal}" = "distclean" ] ; \
	 then printf "\nRemoving build links in $@/ . . . " ; \
	      for link in GNUmakefile Makefile makefile Makefiles bin lib ; \
	      do test ${hL} $@/$$link && rm -f $@/$$link ; \
	      done ; \
	      printf "done!\n\n" ; \
	 fi

################################################################################

libmrcclibs.a : ${tl_libs} ${tl_bins} FORCE
	@printf "Recreating monolithic MRCC library . . . " ; \
	 rm -f libmrcclibs.a ; \
	 find . -name "*.o" | xargs ${tl_ar} libmrcclibs.a ; \
	 printf "done!\n\n"

.PHONY: ${tl_mrcc_pass2}
${tl_mrcc_pass2} : libmrcclibs.a
	@bindir="$(patsubst %-p2, %, $@)"; \
	 test ${hL} $$bindir/libmrcclibs.a || ln -s ../libmrcclibs.a $$bindir; \
	 $(MAKE) -C $$bindir ${sm_goals_tmp} ${sm_flags_tmp}

