

include ../make.inc

all_files = config.o main_module.o diagnose.o integrate.o timing_module.o 

all: with_mpi

with_mpi: main.o parallel_mpi.o ${all_files} 
	$(F90) main.o parallel_mpi.o ${all_files} $(F90FLAGS) $(CDFFLAGS) $(MPIFLAGS) -o ../bin/model.x

without_mpi: main.o parallel_none.o ${all_files} 
	$(F90) main.o parallel_none.o ${all_files} $(F90FLAGS) $(CDFFLAGS) -o ../bin/model.x


clean :
	 rm -f *.o *.mod 

timing_module.o    : timing_module.f90
	 	$(F90) $(F90FLAGS) -c timing_module.f90
main_module.o    : main_module.f90 timing_module.o 
	 	$(F90) $(F90FLAGS) -c main_module.f90
main.o         : main.f90 main_module.o timing_module.o 
	 	$(F90) $(F90FLAGS) -c main.f90
integrate.o    : integrate.f90 main_module.o timing_module.o
	 	$(F90) $(F90FLAGS) -c integrate.f90
diagnose.o  : diagnose.f90 main_module.o 
	 	$(F90) $(F90FLAGS) $(CDFFLAGS) -c diagnose.f90
config.o       : config.f90 main_module.o 
	 	$(F90) $(F90FLAGS)  $(CDFFLAGS)  -c config.f90
parallel_none.o : parallel_none.f90 main_module.o
	 	$(F90) $(F90FLAGS)  -c parallel_none.f90
parallel_mpi.o : parallel_mpi.f90 main_module.o
	 	$(F90) $(F90FLAGS) $(MPIFLAGS) -c parallel_mpi.f90

       
        
