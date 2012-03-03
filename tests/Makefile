include ../GNUmakefile.chssi

# Default variable settings
MPIRUN="mpirun ./xaces3 >./job.out"
RUN_COMPARE="mpirun .xtest_compare "
ENV_VARS=" "

# Machine-specific variable settings

ifeq ($(MACHNAME),submit1)
WORKDIR=/scratch/ufhpc/${USER}
endif

ifeq ($(MACHNAME),submit2)
WORKDIR=/scratch/ufhpc/${USER}
endif

ifeq ($(MACHNAME),chugach)
PBS_LINES="\#PBS -j oe"
PBS_LINES += "\#PBS -N ACES3_TEST"
PBS_LINES += "\#PBS -l mppwidth=64"
PBS_LINES += "\#PBS -l walltime=12:00:00"
PBS_LINES += "\#PBS -q standard"
PBS_LINES += "\#PBS -e ACESIII_tests.err"
PBS_LINES += "\#PBS -o ACESIII_tests.out"
PBS_LINES += "\#PBS -A ARONC24913112"
ENV_VARS="setenv MPICH_GNI_NDREG_LAZYMEM disabled"
SET_ACES_EXE_PATH=setenv ACES_EXE_PATH
MPIRUN="aprun -n 64 ./xaces3 >./job.out"
RUN_COMPARE="aprun -n 64 ./xtest_compare "
endif

ifeq ($(MACHNAME),pingo)
PBS_LINES="\#PBS -j oe"
PBS_LINES += "\#PBS -N ACES3_TEST"
PBS_LINES += "\#PBS -l mppwidth=64"
PBS_LINES += "\#PBS -l walltime=12:00:00"
PBS_LINES += "\#PBS -q standard" 
ENV_VARS="setenv MPICH_PTL_OTHER_EVENTS 20000" 
ENV_VARS +="setenv MPICH_PTL_UNEX_EVENTS  100000" 
ENV_VARS +="setenv MPICH_MSGS_PER_PROC 40000" 
ENV_VARS +="setenv MPICH_UNEX_BUFFER_SIZE 128000000" 
ENV_VARS +="setenv MPICH_MAX_SHORT_MSG_SIZE 64000" 
ENV_VARS +="lfs setstripe ./BLOCKDATA 0 -1 -1" 
ENV_VARS +="lfs setstripe ./BLOCKDATA2 0 -1 -1" 
SET_ACES_EXE_PATH=setenv ACES_EXE_PATH 
MPIRUN="aprun -n 64 ./xaces3 >./job.out"
RUN_COMPARE="aprun -n 64 ./xtest_compare "
endif

ifeq ($(MACHNAME),mjm)
PBS_LINES="\#BSUB -m \"mjm\""
PBS_LINES+="\#BSUB -a \"openmpi\""
PBS_LINES+="\#BSUB -J ACES3_TEST"
PBS_LINES+="\#BSUB -q standard"
PBS_LINES+="\#BSUB -n 64"
PBS_LINES+="\#BSUB -W 12:00"
PBS_LINES+="\#BSUB -P ARLAP02642520"
PBS_LINES+="\#BSUB -e job.err"
ENV_VARS="module load intel"
SET_ACES_EXE_PATH=setenv ACES_EXE_PATH
MPIRUN="mpirun.lsf ./xaces3 >./job.out"
RUN_COMPARE="mpirun.lsf ./xtest_compare "
endif

ifeq ($(MACHNAME),tbird)
WORKDIR=/gscratch1/$(USER)
PBS_LINES="\#PBS -l nodes=64:ppn=1:gscratch1"
PBS_LINES+="\#PBS -l walltime=12:00:00"
PBS_LINES+="\#PBS -A FY093344"
PBS_LINES+="\#PBS -N ACES3_TEST"
SET_ACES_EXE_PATH=setenv ACES_EXE_PATH
MPIRUN="mpiexec ./xaces3 >./job.out"
RUN_COMPARE="mpiexec ./xtest_compare "
endif

ifeq ($(MACHNAME),babbage)
PBS_LINES="\#PBS -q standard"
PBS_LINES+="\#PBS -l walltime=24:00:00"
PBS_LINES+="\#PBS -N ACES3_TEST"
PBS_LINES+="\#PBS -A ARLAP520"
PBS_LINES+="\#PBS -l select=4:ncpus=16:mpiprocs=16"
PBS_LINES+="\#PBS -l place=scatter:excl"
SET_ACES_EXE_PATH=setenv ACES_EXE_PATH
MPIRUN="poe ./xaces3 >./job.out"
RUN_COMPARE="poe ./xtest_compare "
endif

ifeq ($(MACHNAME),einstein)
PBS_LINES="\#PBS -N ACES3_TEST"
PBS_LINES+="\#PBS -o job.out"
PBS_LINES+="\#PBS -e job.err"
PBS_LINES+="\#PBS -A ARLAP02643520"
PBS_LINES+="\#PBS -l walltime=12:00:00"
PBS_LINES+="\#PBS -l mppwidth=64"
PBS_LINES+="\#PBS -l mppnppn=4"
PBS_LINES+="\#PBS -q standard"
SET_ACES_EXE_PATH=setenv ACES_EXE_PATH
MPIRUN="aprun -N 4 -n 64 ./xaces3 >./job.out"
RUN_COMPARE="aprun -N 4 -n 64 ./xtest_compare "
endif

ifeq ($(MACHNAME),davinci)
PBS_LINES="\#PBS -N ACES3_TEST"
PBS_LINES+="\#PBS -o job.out"
PBS_LINES+="\#PBS -e job.err"
PBS_LINES+="\#PBS -A ARLAP02643520"
PBS_LINES+="\#PBS -l walltime=24:00:00"
PBS_LINES+="\#PBS -l select=2:ncpus=32:mpiprocs=32"
PBS_LINES+="\#PBS -l place=scatter:excl"
PBS_LINES+="\#PBS -q standard"
SET_ACES_EXE_PATH=setenv ACES_EXE_PATH
MPIRUN="poe ./xaces3 >./job.out"
RUN_COMPARE="poe ./xtest_compare "
endif

ifeq ($(MACHNAME),hawk)
PBS_LINES="\#PBS -q standard"
PBS_LINES+="\#PBS -l select=16:ncpus=4"
PBS_LINES+="\#PBS -l walltime=24:00:00"
PBS_LINES+="\#PBS -N ACESIII_TEST"
PBS_LINES+="\#PBS -A WPCARLAP02642520"
ENV_VARS="setenv MPI_DEFAULT_SINGLE_COPY_OFF 1"
ENV_VARS+="unsetenv MPI_BUFFER_MAX"
ENV_VARS+="unsetenv MPI_MAPPED_HEAP_SIZE"
ENV_VARS+="unsetenv MPI_MAPPED_STACK_SIZE"
ENV_VARS+="setenv MPI_REQUEST_MAX 100000"
SET_ACES_EXE_PATH=setenv ACES_EXE_PATH
MPIRUN="mpirun -np 64 dplace -s1 -c0-63 ./xaces3 >job.out"
RUN_COMPARE="mpirun -np 64 dplace -s1 -c0-63 ./xtest_compare "
endif

ifeq ($(MACHNAME),midnight)
PBS_LINES="\#PBS -q standard"
PBS_LINES+="\#PBS -l select=16:ncpus=4:node_type=4way"
PBS_LINES+="\#PBS -l walltime=16:00:00.0"
PBS_LINES+="\#PBS -N ACESIII_TEST"
ENV_VARS="module switch PrgEnv PrgEnv.pgi"
SET_ACES_EXE_PATH=setenv ACES_EXE_PATH
MPIRUN="mpirun -np 64 ACES_EXE_PATH=$(ACES3_EXE) ./xaces3 >job.out"
RUN_COMPARE="mpirun -np 64 ACES_EXE_PATH=$(ACES3_EXE) ./xtest_compare "
endif

INSTALL_DIR=$(WORKDIR)/ACESIII_TESTS
ACES3_EXE=$(shell pwd)/../bin
GENBAS=./GENBAS
TESTS=1.1.1.1 1.1.2.1 1.1.2.2 1.1.3.1 1.1.3.2 1.1.4.1 1.1.4.2 \
      1.2.1.1 1.2.2.1 1.2.3.1 1.2.3.2 \
	1.3.1.1 1.3.2.1 \
	1.4.3.1 1.4.3.2 \
	2.1.1.1 2.1.2.1 2.1.2.2 2.1.3.1 2.1.3.2 2.1.4.1 2.1.4.2 \
	2.2.1.1 2.2.2.1 2.2.3.1 2.2.3.2 \
	2.3.1.1 2.3.2.1 \
	2.4.3.1 2.4.3.2 \
	3.1.1.1 3.1.2.1 3.1.3.1 3.1.3.2 3.1.4.1 3.1.4.2 \
	3.2.1.1 3.2.1.2 3.2.2.1 3.2.2.2 3.2.3.1

all: tests scripts

# Build directories of all test jobs.

tests: % : ;
	@cp test_results $(ACES3_EXE) ;
	@rm -rf $(INSTALL_DIR) ; 
	@mkdir $(INSTALL_DIR) ; 
	@for mytest in $(TESTS) ; \
	 do mkdir $(INSTALL_DIR)/$$mytest || exit 1; \
	    cp $(GENBAS) $(ACES3_EXE)/xaces3 $(ACES3_EXE)/xtest_compare $(INSTALL_DIR)/$$mytest || exit 1; \
	    cp zmat_$$mytest $(INSTALL_DIR)/$$mytest/ZMAT || exit 1; \
	 done
	
scripts: % : ;
	@echo \#!/bin/csh >runscript 
	@for line in $(PBS_LINES) ; \
	   do echo $$line >>runscript ; \
	   done
	@for line in $(ENV_VARS) ; \
	   do echo $$line >>runscript ; \
	   done
	@echo $(SET_ACES_EXE_PATH) $(ACES3_EXE) >>runscript ;
	@echo cd $(INSTALL_DIR) >>runscript ;
	@echo "echo Test of ACESIII system >monitor_tests" >>runscript ; 
	@echo "date >>monitor_tests" >>runscript  ; 
	@for test in $(TESTS) ; \
	   do echo >>runscript  ; \
	      echo \#\#\# Test $$test \#\#\# >>runscript ; \
	      echo cd $(INSTALL_DIR)/$$test >>runscript ; \
              echo lfs setstripe BLOCKDATA  0 -1 -1 >>runscript ; \
              echo lfs setstripe BLOCKDATA2 0 -1 -1 >>runscript ; \
	      echo "echo Running $$test >>../monitor_tests" >>runscript ; \
              echo "date >>../monitor_tests" >>runscript ; \
	      echo $(MPIRUN) >>runscript ; \
	      echo "echo Completed $$test >>../monitor_tests" >>runscript ; \
	      echo "$(RUN_COMPARE) $$test >>../monitor_tests" >>runscript ; \
	      echo sleep 10 >>runscript ; \
	   done  

clean: % : ;
	@rm -rf $(INSTALL_DIR) 
	@rm -rf runscript
	@rm -rf ACES3_TEST.o*

libclean binclean archive: % : ;
	

