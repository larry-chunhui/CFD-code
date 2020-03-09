#f77 
#f77=ifort
f77=/scinet/niagara/software/2018a/opt/intel-2018.2/openmpi/3.0.1/bin/mpif90

#opt = -fast
opt=-O3 

srs=  channel.f90  correct.f90 \
deriv_3d.f90  Duvwp.f90 \
fftpack.f90  init.f90 \
lambda2.f90  plas_force.f90 \
setup_3d.f90  solve_3d.f90 \
solve_i_3d.f90 	statistic_kinetic_energy.f90 \
statistic_RS.f90  statistics_3d.f90 \
vbound.f90 \
xyzfft.f90  xyzfft2.f90 \

OBJS=$(srs:.f90=.o)

%.o:%.f90
	$(f77) $(opt) -c $<

default: $(OBJS)
	$(f77) -O3  -o channel.out $(OBJS)


clean:
	rm -f *.out *.o work.* *.pc *.mod

