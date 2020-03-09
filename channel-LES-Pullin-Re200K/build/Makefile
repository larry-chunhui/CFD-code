#Program : Makefile
# Purpose : Makefile for main
# B.C.: periodic in x- and z- direction

HOME = .

#FFTW3 = /home/fftw-3.3
#LAPACK	= /opt/intel/fce/10.1.018    
#
#/opt/intel/composer_xe_2013.0.079/mkl/include/fftw/fftw3.h


#FFTWINC = -I$(FFTW3)/include
#FFTWINC = -I/vol-th/home/zhangtq/software/fftw/fftw-3.3.5/ztqpathmpi/include
#FFTWINC = -I/vol-th/home/zhangtq/software/fftw/fftw-3.3.5/ztqpathmpi/include
FFTWINC = -I/export/home/lansl/liuchunhui/software/fftw-3.3.5/include


#FFTWLIB = -L$(FFTW3)/lib -lfftw3_mpi -lfftw3 -lm
FFTWLIB = -L/export/home/lansl/liuchunhui/software/fftw-3.3.5/lib -lfftw3_mpi -lfftw3 -lm

#LAPACKLIB = -L/opt/intel/fce/10.1.018/lib -llapack  -lblas -lm
#LAPACKLIB = -L/vol-th/intel_composer_xe_2013.0.079_lib -lmkl_intel_lp64 -lmkl_sequential  -lmkl_core -lpthread
LAPACKLIB = -L/usr/local/intel/composerxe-2011/composerxe-2011.5.220/mkl/lib/intel64/libmkl_blacs_ilp64.a -lmkl_intel_lp64 -lmkl_sequential  -lmkl_core -lpthread
#/opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64/libmkl_blacs_ilp64.a
# LAPACKLIB = -L/mkl/lib/intel64/ -lmkl_intel_lp64
    
INCS = $(FFTWINC)
LIBS = $(FFTWLIB) $(LAPACKLIB)




CFLAGS  =   -Wall -O3
#CFLAGS  =  -O3

#               -Mx,119,0x200000
#		-g -xildoff -Wall

CC  = /usr/local/mvapi2/bin/mpicc
F77 =/usr/local/mpi-intel2013/bin/mpif77 
#CC = /opt/openmpi/bin/mpif77


SRCDIR = ../src
VPATH = $(SRCDIR)

EXE = ../main

OBJS = main.o nrutil.o fft.o error_check.o derivatives__2.o derivatives__4.o convective.o convective_div__2.o convective_div__4.o convective_adv__2.o convective_adv__4.o share.o interpolation.o interpolation__les.o interpolation__2.o interpolation__4.o differentiate.o rhs.o time_march.o trisolver.o pentasolver.o poisson_solver_transpose.o correction.o statistics.o profiles.o velgrad_tensor.o bc__get.o les.o les_diagnostics.o spiral.o near_wall_ode.o fasttrig.o parameters.o

$(EXE): $(OBJS)
	$(CC) $(CFLAGS) -o $(EXE) $(OBJS) $(LIBS)

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@ $(INCS)

clean: 
	rm -f $(OBJS)

main.o:	main.c definitions.h parameters.h nrutil.h
parameters.o: parameters.c definitions.h nrutil.h parameters.h
nrutil.o: nrutil.c definitions.h nrutil.h
error_check.o: error_check.c definitions.h fft.h error_check.h
derivatives__2.o: derivatives__2.c definitions.h fft.h derivatives__2.h
derivatives__4.o: derivatives__4.c definitions.h fft.h derivatives__4.h
share.o: share.c definitions.h share.h
fft.o: fft.c definitions.h fft.h
interpolation.o: interpolation.c definitions.h interpolation.h
interpolation__2.o: interpolation__2.c definitions.h interpolation__2.h
interpolation__les.o: interpolation__les.c definitions.h interpolation__les.h
interpolation__4.o: interpolation__4.c definitions.h interpolation__4.h
convective.o: convective.c definitions.h fft.h convective.h
convective_div__2.o: convective_div__2.c definitions.h fft.h convective_div__2.h
convective_adv__2.o: convective_adv__2.c definitions.h fft.h convective_adv__2.h
convective_div__4.o: convective_div__4.c definitions.h fft.h convective_div__4.h
convective_adv__4.o: convective_adv__4.c definitions.h fft.h convective_adv__4.h
differentiate.o: differentiate.c definitions.h fft.h differentiate.h
rhs.o: rhs.c definitions.h rhs.h
time_march.o: time_march.c definitions.h time_march.h
trisolver.o: trisolver.c definitions.h trisolver.h
pentasolver.o: pentasolver.c definitions.h pentasolver.h
correction.o: correction.c definitions.h correction.h
statistics.o: statistics.c definitions.h statistics.h
profiles.o: profiles.c definitions.h profiles.h
velgrad_tensor.o: velgrad_tensor.c definitions.h velgrad_tensor.h
bc__get.o: bc__get.c definitions.h bc__get.h
les.o: les.c definitions.h les.h
les_diagnostics.o: les_diagnostics.c definitions.h les_diagnostics.h
spiral.o: spiral.c definitions.h spiral.h
near_wall_ode.o: near_wall_ode.c definitions.h near_wall_ode.h
poisson_solver_transpose.o: poisson_solver_transpose.c definitions.h poisson_solver_transpose.h
fasttrig.o: fasttrig.c fasttrig.h
