#http://swcarpentry.github.io/make-novice/
ifeq ($(FCOMP),mpifort)
FCFLAGS =  -freal-4-real-8 -O2 -march=native -flto
FCDEBUG = -g -fbacktrace -fcheck=all -fbounds-check -ffpe-trap=invalid,overflow,underflow,denormal
FCBUILD = -Wall -Wextra -pedantic -std=f2008
else ifeq ($(FCOMP),mpiifort)
FCFLAGS = -O2 -r8 -xHOST -ipo -heap-arrays
FCDEBUG = -g -debug -traceback -check all -CB -CU -CA -fpe0
#-traceback -check all -CB -CU -CA -fpe0 -debug
FCBUILD = -g -debug -traceback -check all -CB -CA -CU -warn all -stand f08
endif

SRCS =      ran2.f \
	      utils.f90 \
            memory.f90 \
            constants.f90 \
            photon_vars.f90 \
            thermalConst_mod.f90 \
            iarray.f90 \
            opt_prop.f90 \
            subs.f90 \
            ch_opt.f90 \
            gridset.f90 \
            inttau2.f90 \
            sourceph.f90 \
            3dFD.f90 \
            writer.f90 \
            mcpolar.f90 
                        
OBJECTS = $(SRCS:.f90=.o)

all:	$(PROGRAM)
debug:	FCFLAGS += $(FCDEBUG)
debug:	$(PROGRAM)
build:	FCFLAGS += $(FCBUILD)
build:	$(PROGRAM)

$(PROGRAM):	$(OBJECTS)
	$(FCOMP) $(FCFLAGS) -o $@ $^

%.o:  %.f90
	$(FCOMP) $(FCFLAGS) -c $<

.PHONY:	clean

clean:
	rm -f *.o *.mod *.MOD mcgrid