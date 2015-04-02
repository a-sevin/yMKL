# these values filled in by    yorick -batch make.i
Y_MAKEDIR=/home/sevin/yorick.git/relocate
Y_EXE=/home/sevin/yorick.git/relocate/bin/yorick
Y_EXE_PKGS=

Y_EXE_HOME=/home/sevin/yorick.git/relocate
Y_EXE_SITE=/home/sevin/yorick.git/relocate
Y_HOME_PKG=

# ----------------------------------------------------- optimization flags

# options for make command line, e.g.-   make COPT=-g TGT=exe
COPT=$(COPT_DEFAULT)
TGT=$(DEFAULT_TGT)
# ------------------------------------------------ macros for this package

PKG_NAME = libmkl
PKG_I=libmkl.i

USE_OPENMP=1

PKG_CFLAGS=-Wall -fPIC -O3 -fno-strict-aliasing
PKG_CFLAGS+=-I$(MKLROOT)/include -fopenmp  -DMKL_ILP64
PKG_DEPLIBS=-fPIC -fopenmp
PKG_DEPLIBS+= -lmkl_intel_thread -lmkl_core -lpthread 
PKG_DEPLIBS+= -lstdc++ -lm -liomp5 -lgfortran -lmkl_intel_ilp64 
PKG_DEPLIBS+= -L$(MKLROOT)/lib/intel64

OBJS=libmkl.o

# -------------------------------- standard targets and rules (in Makepkg)

#--compiler-options -fpermissive 

# set macros Makepkg uses in target and dependency names
# DLL_TARGETS, LIB_TARGETS, EXE_TARGETS
# are any additional targets (defined below) prerequisite to
# the plugin library, archive library, and executable, respectively
PKG_I_DEPS=$(PKG_I)
Y_DISTMAKE=distmake

include $(Y_MAKEDIR)/Make.cfg
include $(Y_MAKEDIR)/Makepkg
include $(Y_MAKEDIR)/Make$(TGT)

ifeq ($(USE_OPENMP),1)
	CC=icc
endif

# override macros Makepkg sets for rules and other macros
# Y_HOME and Y_SITE in Make.cfg may not be correct (e.g.- relocatable)
Y_HOME=$(Y_EXE_HOME)
Y_SITE=$(Y_EXE_SITE)

# reduce chance of yorick-1.5 corrupting this Makefile
MAKE_TEMPLATE = protect-against-1.5

# ------------------------------------- targets and rules for this package
