# makefile overrides
# OS:       macOS
# Compiler: gfortran X.X/Clang
# OpenMP:   enabled
#

CC=gcc
CXX=g++
FC=gfortran

ifeq ($(PREFIX),)
    SOM_INSTALL_DIR=/usr/local/lib
endif


# MATLAB interface:
FDIR=$$(dirname `gfortran --print-file-name libgfortran.dylib`)
MFLAGS +=-L${FDIR}
MEX = $(shell ls -d /Applications/MATLAB_R* | sort | tail -1)/bin/mex


