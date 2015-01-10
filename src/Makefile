# Makefile for the NRmin project

# fortran compiler
FC = gfortran

# flags
FCFLAGS = -g -fbounds-check -Wall -Wextra # debug (comment if not needed)

# object file names
OBJS = mod_chi2.o mod_random.o NRmin_mod.o NRmin.o

NRmin.bin: $(OBJS)
	$(FC) -o NRmin.bin $(FCFLAGS) $(OBJS)

%.o: %.f90
	$(FC) -c $(FCFLAGS) $<

clean:
	rm -f *.o *.mod *.bin

