# makefile.sol

all : path

# Compilers
F90 = g95

# Compilation options
OPTF = -O3 -c
OPTO = -O3 -o

# Module
MOD = Defs.o

# Objects
OBJ = path.o HOMPACK90.o Defs.o ftir.o hybrd.o rkf45.o LAPACK.o user.o

# Dependencies
path : $(OBJ)
	$(F90) $(OBJ) $(OBJTIR) $(OPTO) path
	mv path ../bin

path.o: path.f90 $(MOD) HOMPACK90.o
	$(F90) $(OPTF) path.f90

Defs.o: Defs.f90
	$(F90) $(OPTF) Defs.f90

LAPACK.o : LAPACK.f
	$(F90) $(OPTF) LAPACK.f

HOMPACK90.o : HOMPACK90.f $(MOD)
	$(F90) $(OPTF) HOMPACK90.f

hybrd.o: hybrd.f
	$(F90) $(OPTF) hybrd.f

rkf45.o: rkf45.f
	$(F90) $(OPTF) rkf45.f

ftir.o: ftir.f90 $(MOD)
	$(F90) $(OPTF) ftir.f90

user.o: user.f90 $(MOD)
	$(F90) $(OPTF) user.f90

clean:	
	rm *.mod *.o
