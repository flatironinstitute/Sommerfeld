#Compile the .f file:
gfortran -c -O3 -march=native -std=legacy helmholtz2d/sommerfeld2drouts_alpert.f -o helmholtz2d/sommerfeld2drouts_alpert.o
gfortran -c -O3 -march=native -std=legacy common/alpertsqrt.f -o common/alpertsqrt.o
gfortran -c -O3 -march=native -std=legacy common/discretizer.f -o common/discretizer.o
gfortran -c -O3 -march=native -std=legacy ../../FMM3D/src/Common/legeexps.f -o ../../FMM3D/src/Common/legeexps.o


#Generate the .m file from the .mw file
../../mwrap/mwrap -c99complex -list -mex somm2drouts -mb somm2drouts.mw
#Generate a .c gateway
../../mwrap/mwrap -c99complex -mex somm2drouts -c somm2drouts.c somm2drouts.mw
#Generate mexmaci file to run on Mac
/Applications/MATLAB_R2023a.app/bin/mex -v somm2drouts.c helmholtz2d/sommerfeld2drouts_alpert.o common/alpertsqrt.o common/discretizer.o /usr/local/lib/libfmm2d.a ../../FMM3D/src/Common/legeexps.o -compatibleArrayDims -DMWF77_UNDERSCORE1 -output somm2drouts -lgfortran -L/usr/local/Cellar/gcc/13.2.0/lib/gcc/current/


