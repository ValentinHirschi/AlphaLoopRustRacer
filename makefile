GFORTRAN = gfortran
GCC = gcc
GPP = g++
all: racer

ifeq ($(OPTLEVEL),)
OPTLEVEL=3
endif

ifeq ($(OPTLEVEL),0)
FLAGS = -w -fbounds-check
else
FLAGS = -ffast-math -funroll-loops -march=native
endif

ifeq ($(OPTLEVEL),0)
CFLAGS = 
else
CFLAGS = -ffast-math -funroll-loops -march=native -fcx-fortran-rules -fcx-limited-range
endif

ifeq ($(OPTLEVEL),0)
CPPFLAGS = -std=c++20
else
CPPFLAGS = -std=c++20 -ffast-math -funroll-loops -march=native -fcx-fortran-rules -fcx-limited-range -fext-numeric-literals
endif

ifeq ($(MPFRPATH),)
MPFRPATH=/Users/vjhirsch/HEP_programs/mpfr-4.1.0_build
endif

ifeq ($(GMPPATH),)
GMPPATH=/Users/vjhirsch/HEP_programs/gmp-6.2.1_build
endif

ifeq ($(FORTRANMPFRPATH),)
FORTRANMPFRPATH=/Users/vjhirsch/HEP_programs/mpfr_fortran/mpfun20-mpfr-v24/fortran-var2
endif

ifeq ($(MPPPATH),)
MPPPATH=/Users/vjhirsch/MG5/3.0.2.py3/PLUGIN/tderivative_alphaloop/libraries/form/mppp-0.26
endif

$(info    Compiling with OPTLEVEL $(OPTLEVEL))

%.o : %.f
	$(GFORTRAN) -c $^ -O$(OPTLEVEL) -ffree-form -ffree-line-length-none $(FLAGS)

%.o : %.f90
	$(GFORTRAN) -c $^ -O$(OPTLEVEL) -ffree-form -ffree-line-length-none $(FLAGS) -I$(FORTRANMPFRPATH)

%.o : %.c
	$(GCC) -c $^ -O$(OPTLEVEL) $(CFLAGS)

%.o : %.cpp
	$(GPP) -c $^ -O$(OPTLEVEL) $(CPPFLAGS) -I$(MPPPATH)/include -I$(GMPPATH)/include -I$(MPFRPATH)/include

racer : racer.o
	$(GFORTRAN) -o racer racer.o -O$(OPTLEVEL) $(FLAGS)

racer_f128 : racer_f128.o
	$(GFORTRAN) -o racer_f128 racer_f128.o -O$(OPTLEVEL) $(FLAGS)

racer_mpfr : racer_mpfr.o
	$(GFORTRAN) -o racer_mpfr racer_mpfr.o $(FORTRANMPFRPATH)/mpmodule.o $(FORTRANMPFRPATH)/mpfuna.o $(FORTRANMPFRPATH)/mpfune.o $(FORTRANMPFRPATH)/mpfunf.o $(FORTRANMPFRPATH)/mpfung2.o $(FORTRANMPFRPATH)/mpfunh2.o $(FORTRANMPFRPATH)/mpinterface.o $(FORTRANMPFRPATH)/second.o -O$(OPTLEVEL) -fno-underscoring $(FLAGS) -L$(MPFRPATH)/lib -L$(GMPPATH)/lib -I$(FORTRANMPFRPATH) -lmpfr

racer_c : racer_c.o
	$(GCC) -o racer_c racer_c.o -O$(OPTLEVEL) $(CFLAGS)

racer_c_f128 : racer_c_f128.o
	$(GPP) -o racer_c_f128 racer_c_f128.o -O$(OPTLEVEL) -lmp++ -lquadmath -L$(MPPPATH) $(CPPFLAGS)

clean:
	rm -f racer.o racer

clean_mpfr:
	rm -f racer_mpfr.o racer_mpfr

clean_f128:
	rm -f racer_f128.o racer_f128

clean_c:
	rm -f racer_c.o racer_c

clean_c_f128:
	rm -f racer_c_f128.o racer_c_f128