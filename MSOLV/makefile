MA48_LIB = 1

utils_dir = ../UTILS

MODDIR = mod
OBJDIR = obj
SRCDIR = src
LIBDIR = lib
INCDIR = inc

GENLIB  =
GENLIB += $(LIBDIR)/librefblas.a
GENLIB += $(LIBDIR)/liblapack.a

MUMLIB  =
MUMLIB += $(LIBDIR)/libmpiseq.a
MUMLIB += $(LIBDIR)/libdmumps.a
MUMLIB += $(LIBDIR)/libmumps_common.a
MUMLIB += $(LIBDIR)/libmetis.a

MUMLIB += $(LIBDIR)/libscotch.a
MUMLIB += $(LIBDIR)/libpord.a
MUMLIB += $(LIBDIR)/libesmumps.a
MUMLIB += $(LIBDIR)/libscotcherr.a

UMFLIB  =
UMFLIB += $(LIBDIR)/libumfpack.a
UMFLIB += $(LIBDIR)/libamd.a
UMFLIB += $(LIBDIR)/libsuitesparseconfig.a
UMFLIB += $(LIBDIR)/libcholmod.a
UMFLIB += $(LIBDIR)/libcolamd.a
UMFLIB += $(LIBDIR)/libccolamd.a
UMFLIB += $(LIBDIR)/libcamd.a

SUPLIB  =
SUPLIB += $(LIBDIR)/libsuperlu.a

M48LIB  =
ifeq ($(MA48_LIB),1)
   M48LIB += $(LIBDIR)/libhsl_ma48.a
endif

SMALIB = -Wl,--start-group $(GENLIB) $(M48LIB) $(MUMLIB) $(SUPLIB) $(UMFLIB) $(FFTLIB) -Wl,--end-group

EXT = f90

VPATH = $(MODDIR):$(OBJDIR):$(SRCDIR)

EXE = prg

SRC = $(notdir $(wildcard $(SRCDIR)/*.$(EXT)))

OBJ      = $(SRC:.$(EXT)=.o)
ALL_OBJS = $(addprefix $(OBJDIR)/, *.o) $(addprefix $(utils_dir)/$(OBJDIR)/, *.o) 


CFLAGS  = -I$(MODDIR) -I$(INCDIR) -I$(utils_dir)/$(MODDIR)
CFLAGS += -cpp -DWITH_MA48=$(MA48_LIB) -ffree-form -ffree-line-length-none -march=native -fimplicit-none -fall-intrinsics -fmax-errors=1 -fopenmp

LFLAGS  = $(SMALIB)
LFLAGS += -lpthread -lm -lgomp


ifneq ('$(DEBUG)','')
	CFLAGS += -Og -g -Wall -Wextra -fbacktrace -pedantic -fbounds-check -Wuninitialized -ffpe-summary=invalid,zero,overflow,underflow -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan
else
	CFLAGS += -O3 -ffpe-summary=none
endif

ifneq ('$(GPROF)','')
	CFLAGS += -pg -g
	LFLAGS += -pg
endif

ifndef ($(FORT))
	FORT = gfortran-7
endif

%.o: %.$(EXT)
	$(FORT) $(CFLAGS) -c $< -o $(OBJDIR)/$@
	-mv -f *.mod $(MODDIR)

$(EXE): $(OBJ)
	$(FORT) $(ALL_OBJS) $(LFLAGS) -o $(EXE)
	rm $(OBJDIR)/prg.o

mod_gen_par.o :

umfpack.o :
superlu.o :
dmumps_struc.o :
hsl_common90.o :
hsl_ddeps90.o :
hsl_ma48d.o : hsl_common90.o hsl_ddeps90.o

mod_solver.o : umfpack.o hsl_ma48d.o superlu.o dmumps_struc.o mod_gen_par.o

prg.o : mod_solver.o

#--------------------------------------------------
utils:
	cd $(utils_dir) ; make all

#--------------------------------------------------
.PHONY: clean debug gprof all

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(MODDIR)/*.mod
	rm -f prg

all: clean $(EXE)


