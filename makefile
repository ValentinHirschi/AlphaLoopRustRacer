GFORTRAN = gfortran
all: racer

ifeq ($(OPTLEVEL),)
OPTLEVEL=3
endif

ifeq ($(OPTLEVEL),0)
FLAGS = -w -fbounds-check
else
FLAGS = -ffast-math -funroll-loops -march=native
endif

ifeq ($(MPFRLIBPATH),)
MPFRLIBPATH=/Users/vjhirsch/HEP_programs/mpfr-4.1.0_build/lib
endif

ifeq ($(GMPLIBPATH),)
GMPLIBPATH=/Users/vjhirsch/HEP_programs/gmp-6.2.1_build/lib
endif

ifeq ($(FORTRANMPFRPATH),)
FORTRANMPFRPATH=/Users/vjhirsch/HEP_programs/mpfr_fortran/mpfun20-mpfr-v24/fortran-var2
endif


$(info    Compiling with OPTLEVEL $(OPTLEVEL))

%.o : %.f
	$(GFORTRAN) -c $^ -O$(OPTLEVEL) -ffree-form -ffree-line-length-none $(FLAGS)

%.o : %.f90
	$(GFORTRAN) -c $^ -O$(OPTLEVEL) -ffree-form -ffree-line-length-none $(FLAGS) -I$(FORTRANMPFRPATH)

racer : racer.o
	$(GFORTRAN) -o racer racer.o -O$(OPTLEVEL) $(FLAGS)

racer_f128 : racer_f128.o
	$(GFORTRAN) -o racer_f128 racer_f128.o -O$(OPTLEVEL) $(FLAGS)

racer_mpfr : racer_mpfr.o
	$(GFORTRAN) -o racer_mpfr racer_mpfr.o $(FORTRANMPFRPATH)/mpmodule.o $(FORTRANMPFRPATH)/mpfuna.o $(FORTRANMPFRPATH)/mpfune.o $(FORTRANMPFRPATH)/mpfunf.o $(FORTRANMPFRPATH)/mpfung2.o $(FORTRANMPFRPATH)/mpfunh2.o $(FORTRANMPFRPATH)/mpinterface.o $(FORTRANMPFRPATH)/second.o -O$(OPTLEVEL) -fno-underscoring $(FLAGS) -L$(MPFRLIBPATH) -L$(GMPLIBPATH) -I$(FORTRANMPFRPATH) -lmpfr

clean:
	rm -f racer.o racer

clean_mpfr:
	rm -f racer_mpfr.o racer_mpfr

clean_f128:
	rm -f racer_f128.o racer_f128