#!/bin/bash
#################################################################################
#
# MAKEFILE WITH AUTOMATIC DEPENDENCIES for program prova.f90
#
#################################################################################
#
# Date of last update: 08/05/2013
#
#
# Author: Marco Giometto
#
# Environment: Ubuntu 12.04
#
#
# Description:
# 
# Some precompiler flags are there and used in the OUTPUT_SPC.f90, this is done
# so that anyone who uses the code can implement his output in specific standalone 
# modules without affecting the memory allocation at compile time of the code.
#
#
#################################################################################

#preprocessor stuff --------------------- please modify
BDG = yes
MPI = yes
SCALARS = yes

#define shell
SHELL = /bin/sh

#architecture ---------------------------- please modify
#ARCH = PC_GNU
ARCH = PC_INTEL
#ARCH = BELLATRIX_INTEL
#ARCH = CSCS_INTEL
#ARCH = DENEB_INTEL
#ARCH = MAMMOUTH

#mpi library path
MPIMOD = /usr/local/lib/mpi.mod

#executable name
EXE = program.exe

#mode
MODE = TRACEBACK
#MODE = AGGRESSIVE

#dependencies creator
MAKEDEP = ./src/makedepf90

#directory for the .o files
OBJDIR = ./src/.obj
MODDIR = ./src/.mod

#directories for sources
SRCS = ./src/*.f90
SRCS += ./src/src_sgs/*.f90
SRCS += ./src/src_sc/*.f90
SRCS += ./src/src_bdg/*.f90
SRCS += ./src/src_mpi/*.f90

FPPF = -DLES			#this does nothing at the moment
#FPPF += -DUCL
FPPF += -DMPI
FPPF += -DBDG


#################################################################################
#  setting compilers and libraries depending on the pc
#################################################################################


#settings for PC with GNU compiler (f90)
ifeq ($(ARCH),PC_GNU)
  ifeq ($(MPI),yes)
    FC = mpif90
  else
    FC = gfortran
  endif
  FFLAGS = -O3 -cpp -m64 -g -fopenmp -fbounds-check -I$(OBJDIR) -J$(MODDIR) 
  # -Wall -fbounds-check -fbacktrace 
  LIBPATH = -L/home/mg/Phd/libraries/bin_gfortran
  LIBS = $(LIBPATH) -lfftw3 
endif

#settings for PC with INTEL compiler (ifx)
ifeq ($(ARCH),PC_INTEL)
   ifeq ($(MPI),yes)
      FC = mpiifx -DMPI
   else
      FC = ifx
   endif
   FFLAGS = -fpp -O3 -qmkl -assume byterecl -qopenmp -no-wrap-margin $(FPPF) -module $(MODDIR) -I$(OBJDIR)
   #FFLAGS = -fpp -g -O0 -traceback -check bounds -assume byterecl -fpe0 -debug all -no-wrap-margin -qmkl -openmp -module $(MODDIR) -I$(OBJDIR)
   LDFLAGS = #-mcmodel medium
   LIBPATH = #-L/home/mg/Phd/libraries/bin_gfortran
   LIBS = $(LIBPATH) #-qmkl -lmkl_lapack95_lp64
endif

#settings for ROSA with INTEL compiler (ifx)
ifeq ($(ARCH),DENEB_INTEL)
  ifeq ($(MPI),yes)
    FC = mpiifx
  else
    FC = ifx
  endif
  FFLAGS = -fpp $(FPPF) -qmkl -assume byterecl -O3 -qopenmp -module $(MODDIR) -I$(OBJDIR) -axSSE4.2 
  #LDFLAGS = -mcmodel=medium
  LIBPATH = -L/ssoft/fftw/3.3.4/RH6/intel-15.0.0/x86_E5v2/intelmpi/lib
  LIBS = $(LIBPATH) -lfftw3 
endif


#settings for ROSA with INTEL compiler (ifx)
ifeq ($(ARCH),BELLATRIX_INTEL)
  ifeq ($(MPI),yes)
    FC = mpiifx
  else
    FC = ifx
  endif
  FFLAGS = -fpp $(FPPF) -qmkl -check bounds -traceback -assume byterecl -O3 -openmp -module $(MODDIR) -I$(OBJDIR) -axSSE4.2 
  #LDFLAGS = -mcmodel=medium
  LIBPATH = -L/ssoft/fftw/3.3.4/RH6/intel-15.0.0/x86_E5v2/intelmpi/lib
  LIBS = $(LIBPATH) -lfftw3 
endif

#settings for DORA with INTEL compiler (ifx)
ifeq ($(ARCH),CSCS_INTEL)
  FC = ftn
  #FC = scorep --user --static ftn
  FFLAGS = -fpp -O3 -assume byterecl -shared-intel -dynamic -qmkl -qopenmp -no-wrap-margin $(FPPF) -module $(MODDIR) -I$(OBJDIR)
  #FFLAGS = -fpp -g -O0 -traceback -assume byterecl -dynamic -check bounds -fpe0 -no-wrap-margin -qmkl -qopenmp $(FPPF) -module $(MODDIR) -I$(OBJDIR) 
  #FFLAGS = -fpp $(FPPF) -g -debug all -qmkl -traceback -check bounds -assume byterecl -no-wrap-margin -qopenmp -dynamic -shared-intel -mcmodel=large -qopt-streaming-stores always -module $(MODDIR) -I$(OBJDIR)
  LDFLAGS = 
  #LIBPATH = -L/opt/software/fftw/3.3.3/intel-intelmpi/13.0.1-4.1.0/lib
  LIBS = $(LIBPATH) #-lfftw3
endif

#settings for ROSA with INTEL compiler (ifx)
ifeq ($(ARCH),MAMMOUTH)
  FC = mpif90
  FFLAGS = -fpp $(FPPF) -qmkl -traceback -check bounds -assume byterecl -O3 -openmp -module $(MODDIR) -I$(OBJDIR)
  #LDFLAGS = -mcmodel=medium
  #LIBPATH = -L/opt/software/fftw/3.3.3/intel-intelmpi/13.0.1-4.1.0/lib
  LIBS = $(LIBPATH) -lfftw3
endif

#################################################################################
#  various checks and launch the makefile for compilation
#################################################################################

# preprocessor flags settings
COMPSTR = '$$(FC) -c -o $$@ $$(FFLAGS) $$<'

#include the dependency list created by makedepf90 below
include .depend

.depend: $(SRC)
	#cp $(MPIMOD) ./src/.mod/
	$(MAKEDEP) -r $(COMPSTR) -b $(OBJDIR) -o $(EXE) $(SRCS) > .depend


#################################################################################
# cleaning and extra useful commands
#################################################################################

.PHONY : clean

#cleaning out everything
clean:
	rm .depend
	rm ./src/.obj/*
	rm ./src/.mod/*
	rm program.exe
	clear
	@echo "All unnecessary files wiped out!"


