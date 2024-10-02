
# # This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it
# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create
# the file make.inc, which overrides the defaults below (which are
# for ubunutu linux/gcc system).


# compiler, and linking from C, fortran
CC=gcc
CXX=g++
FC=gfortran


# set compiler flags for c and fortran
FFLAGS= -fPIC -O3 -march=native -funroll-loops -std=legacy -w
FFLAGS_DYN= -shared -fPIC
CFLAGS= -fPIC -O3 -march=native -funroll-loops -std=gnu17

# set linking libraries
CLIBS = -lgfortran -lm -ldl
LIBS = -lm

# flags for MATLAB MEX compilation..
MFLAGS=-compatibleArrayDims -DMWF77_UNDERSCORE1 "CFLAGS=-std=gnu17"
MWFLAGS=-c99complex

# location of MATLAB's mex compiler
MEX=mex

# For experts, location of Mwrap executable
MWRAP=../../mwrap/mwrap
MEXLIBS=-lm -lstdc++ -ldl -lgfortran

SOM_INSTALL_DIR=$(PREFIX)
ifeq ($(PREFIX),)
	SOM_INSTALL_DIR = ${HOME}/lib
endif

DYLIBS = $(LIBS)

LIBNAME=libsommerfeld2d
DYNAMICLIB = $(LIBNAME).so
STATICLIB = $(LIBNAME).a
LIMPLIB = $(DYNAMICLIB)

LLINKLIB = -lsommerfeld2d



# For your OS, override the above by placing make variables in make.inc
-include make.inc

LIBDIR = lib-static

# objects to compile
#
# Common objects
COM = src/common
COMOBJS = $(COM)/alpertsqrt.o $(COM)/discretizer.o \
	$(COM)/legeexps.o $(COM)/prini.o 

# Helmholtz objects
HELM = src/helmholtz2d
HOBJS = $(HELM)/sommerfeld2drouts_alpert.o 

# Test objects
TOBJS = $(COM)/hkrand.o $(COM)/dlaran.o


OBJS = $(COMOBJS) $(HOBJS) $(LOBJS) $(STOBJS) $(EMOBJS)

.PHONY: usage install lib test all matlab test-dyn matlab-dyn 

default: usage

usage:
	@echo "------------------------------------------------------------------------"
	@echo "Makefile for sommerfeld. Specify what to make:"
	@echo ""
	@echo "  make install        compile and install the main library"
	@echo "  make install PREFIX=(INSTALL_DIR)  "
	@echo "                      compile and install the main library at custom"
	@echo "                      location given by PREFIX"
	@echo "  make lib            compile the main library (in lib/ and lib-static/)"
	@echo "  make test           compile and run validation tests"
	@echo "                      (will take a few seconds)"
	@echo "  make matlab         compile matlab interfaces"
	@echo "  make test-dyn       test successful installation by validation tests"
	@echo "                      linked to dynamic library (will take a couple of mins)"
	@echo "  make matlab-dyn     compile matlab interfaces with dynamic library linking"
	@echo "  make objclean       removal all object files, preserving lib & MEX"
	@echo "  make clean          also remove lib, MEX, py, and demo executables"
	@echo "  make mex            generate matlab interfaces"
	@echo "                      (for expert users only, requires mwrap)"
	@echo ""
	@echo "For faster (multicore) making, append the flag -j"
	@echo "------------------------------------------------------------------------"


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.cpp %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@
%.o: %.f 
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90
	$(FC) -c $(FFLAGS) $< -o $@

# build the library...
lib: $(STATICLIB) $(DYNAMICLIB)
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, single-threaded versions"


install: $(STATICLIB) $(DYNAMICLIB)
	echo $(SOM_INSTALL_DIR)
	mkdir -p $(SOM_INSTALL_DIR)
	cp -f lib/$(DYNAMICLIB) $(SOM_INSTALL_DIR)/
	cp -f lib-static/$(STATICLIB) $(SOM_INSTALL_DIR)/
	[ ! -f lib/$(LIMPLIB) ] || cp -f lib/$(LIMPLIB) $(SOM_INSTALL_DIR)/
	@echo "Make sure to include " $(SOM_INSTALL_DIR) " in the appropriate path variable"
	@echo "    LD_LIBRARY_PATH on Linux"
	@echo "    PATH on windows"
	@echo "    DYLD_LIBRARY_PATH on Mac OSX (not needed if default installation directory is used"
	@echo " "
	@echo "In order to link against the dynamic library, use -L"$(SOM_INSTALL_DIR) " -lsommerfeld2d"


$(STATICLIB): $(OBJS)
	ar rcs $(STATICLIB) $(OBJS)
	mv $(STATICLIB) lib-static/
$(DYNAMICLIB): $(OBJS)
	$(FC) $(FFLAGS_DYN) $(OBJS) -o $(DYNAMICLIB) $(DYLIBS)
	mv $(DYNAMICLIB) lib/
	[ ! -f $(LIMPLIB) ] || mv $(LIMPLIB) lib/


# matlab..
MWRAPFILE = somm2drouts
GATEWAY = $(MWRAPFILE)

matlab:	$(STATICLIB) matlab/$(GATEWAY).c 
	$(MEX) matlab/$(GATEWAY).c lib-static/$(STATICLIB) $(MFLAGS) \
	-output matlab/somm2drouts $(MEXLIBS)

matlab-dyn:	$(DYNAMICLIB) matlab/$(GATEWAY).c matlab/$(GATEWAY2).c
	$(MEX) matlab/$(GATEWAY).c $(MFLAGS) \
	-output matlab/somm2drouts $(MEXLIBS) -L$(SOM_INSTALL_DIR) $(LLINKLIB)

mex:  $(STATICLIB)
	cd matlab; $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY) -mb $(MWRAPFILE).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY) -c $(GATEWAY).c $(MWRAPFILE).mw;\
	$(MEX) $(GATEWAY).c ../lib-static/$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE) \
	$(MEXLIBS); \

## todo: continue from here

# testing routines
#
test: $(STATICLIB) $(TOBJS) test/somm test/dielec 
	(cd test; ./run_tests.sh)
	cat print_testres.txt
	rm print_testres.txt

test-dyn: $(DYNAMICLIB) $(TOBJS) test/somm-dyn test/dielec-dyn 
	(cd test/Helmholtz; ./run_tests.sh)
	cat print_testres.txt
	rm print_testres.txt

test/somm:
	$(FC) $(FFLAGS) test/sommerfeldsqrt5_dr.f $(TOBJS) -o test/int2-somm lib-static/$(STATICLIB) $(LIBS) -L$(SOM_INSTALL_DIR) -lfmm2d

test/dielec:
	$(FC) $(FFLAGS) test/dielectric_dr.f $(TOBJS) -o test/int2-dielec lib-static/$(STATICLIB) $(LIBS) -L$(SOM_INSTALL_DIR) -lfmm2d


## Linking against dynamic libraries
#
#
test/somm-dyn:
	$(FC) $(FFLAGS) test/sommerfeldsqrt5_dr.f $(TOBJS) -o test/int2-somm -L$(SOM_INSTALL_DIR) $(LLINKLIB) -lfmm2d

test/dielec-dyn:
	$(FC) $(FFLAGS) test/dielectric_dr.f $(TOBJS) -o test/int2-dielec -L$(SOM_INSTALL_DIR) $(LLINKLIB) -lfmm2d

clean: objclean
	rm -f lib-static/*.a lib/*.so lib/*.dll lib/*.lib
	rm -rf matlab/*.mex*
	rm -f test/int2-*

objclean:
	rm -f $(OBJS) $(COBJS) $(TOBJS)
	rm -f test/*.o 
