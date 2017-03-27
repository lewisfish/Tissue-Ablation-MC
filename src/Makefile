FCOMP = mpifort
FCFLAGS = -O2 -freal-4-real-8 
#FCFLAGS = -O2 -freal-4-real-8 -g -fcheck=all -fbounds-check -fbacktrace 
#-ffpe-trap=zero,overflow,underflow 
#FCFLAGS = -Wall -Wextra -pedantic
#LDFLAGS = 
	#goes after $^ -pg 
#-pg goes between fcomp and fcflags
PROGRAM =  mcgrid

SRCS =      ran2.f \
            constants.f90 \
            photon_vars.f90 \
            iarray.f90 \
            opt_prop.f90 \
            subs.f90 \
            ch_opt.f90 \
            density.f90 \
            gridset.f90 \
            inttau2.f90 \
            stokes.f90 \
            sourceph.f90 \
            writer.f90 \
            mcpolar.f90 
            
OBJECTS = $(SRCS:.f90=.o)

all:	$(PROGRAM)

$(PROGRAM):	$(OBJECTS)
	$(FCOMP) $(FCFLAGS) -o $@ $^ 
	
%.o:  %.f90
	$(FCOMP)  $(FCFLAGS) -c $<

.PHONY:	clean

clean:
	rm -f *.o *.mod *.MOD mcgridf mcgrid
