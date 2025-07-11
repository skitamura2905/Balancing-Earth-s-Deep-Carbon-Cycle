COMPILER=mpif90
COMPILER_FLAGS=-O2 -fdefault-real-8
MODULES=phys_par_mantle.mod phys_par_core.mod
EXECUTABLES=earth_evo

all: $(EXECUTABLES)

phys_par_mantle.mod : phys_par_mantle_mod.f90
	$(COMPILER) $(COMPILER_FLAGS) -c phys_par_mantle_mod.f90
phys_par_core.mod : phys_par_core_mod.f90
	$(COMPILER) $(COMPILER_FLAGS) -c phys_par_core_mod.f90

earth_evo : $(OBJECTS) main.f90 phys_par_mantle.mod phys_par_core.mod
	$(COMPILER) $(COMPILER_FLAGS) -c main.f90
	$(COMPILER) -o earth_evo main.o $(OBJECTS)

clean:
	rm -rf earth_evo *.o *.mod 
