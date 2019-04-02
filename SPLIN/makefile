include ../makefile.inc

MODDIR = mod
OBJDIR = obj
SRCDIR = src

EXT = f90

VPATH = $(MODDIR):$(OBJDIR):$(SRCDIR)

EXE = prg

SRC = $(notdir $(wildcard $(SRCDIR)/*.$(EXT)))

OBJ      = $(SRC:.$(EXT)=.o)
ALL_OBJS = $(addprefix $(OBJDIR)/, *.o)

CFLAGS  = -I$(MODDIR) -fopenmp
CFLAGS += -ffree-form -ffree-line-length-none -march=native -fimplicit-none -fall-intrinsics -fmax-errors=1 -finit-real=nan -ffpe-summary=none

LFLAGS  = -lpthread -lm -lgomp

ifneq ('$(DEBUG)','')
	CFLAGS += -Og -g -Wall -Wextra -fbacktrace -pedantic -fbounds-check -Wuninitialized -fimplicit-none
else
	CFLAGS += -O3
endif

ifneq ('$(GPROF)','')
	CFLAGS += -pg -g
	LFLAGS += -pg
endif

%.o: %.$(EXT)
	$(FORT) $(CFLAGS) -c $< -o $(OBJDIR)/$@
	-mv -f *.mod $(MODDIR)

$(EXE): $(OBJ)
	$(FORT) $(ALL_OBJS) $(LFLAGS) -o $(EXE)
	rm $(OBJDIR)/prg.o

mod_bspline.o :
prg.o : mod_bspline.o

#--------------------------------------------------
.PHONY: clean debug gprof all last

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(MODDIR)/*.mod
	rm -f prg

debug:
	make "DEBUG=TRUE"

gprof:
	make "GPROF=TRUE"

all:
	(make clean ; make)

