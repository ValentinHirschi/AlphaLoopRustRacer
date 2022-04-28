GFORTRAN = gfortran
all: racer

ifeq ($(OPTLEVEL),)
OPTLEVEL=3
endif

ifeq ($(OPTLEVEL),0)
FLAGS = -w -fbounds-check
else
FLAGS = -ffast-math
endif

$(info    Compiling with OPTLEVEL $(OPTLEVEL))

%.o : %.f
	$(GFORTRAN) -c $^ -O$(OPTLEVEL) -ffree-form -ffree-line-length-none $(FLAGS)

racer : racer.o
	$(GFORTRAN) -o racer racer.o -O$(OPTLEVEL) $(FLAGS)

clean:
	rm -f racer.o racer