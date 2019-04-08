include ../makefile.inc

MODDIR = mod
OBJDIR = obj
SRCDIR = src
INCDIR = inc
LIBDIR = lib

ALL_LIBS  = $(addprefix $(msolv_dir)/$(LIBDIR)/, *.a)

SMALIB = -Wl,--start-group $(ALL_LIBS) -Wl,--end-group

EXT = f90

VPATH = $(MODDIR):$(OBJDIR):$(SRCDIR)

EXE = main

SRC = $(notdir $(wildcard $(SRCDIR)/*.$(EXT)))

LFLAGS  = $(SMALIB)
LFLAGS += -lpthread -lm -lgomp


OBJ       = $(SRC:.$(EXT)=.o)
ALL_OBJS  = $(addprefix $(OBJDIR)/, *.o)
ALL_OBJS += $(addprefix $(utils_dir)/$(OBJDIR)/, *.o) 
ALL_OBJS += $(addprefix $(digis_dir)/$(OBJDIR)/, *.o) 
ALL_OBJS += $(addprefix $(msolv_dir)/$(OBJDIR)/, *.o) 
ALL_OBJS += $(addprefix $(splin_dir)/$(OBJDIR)/, *.o) 

ALL_MODS  = -I$(MODDIR)
ALL_MODS += -I$(utils_dir)/$(MODDIR)
ALL_MODS += -I$(digis_dir)/$(MODDIR)
ALL_MODS += -I$(msolv_dir)/$(MODDIR)
ALL_MODS += -I$(splin_dir)/$(MODDIR)

CFLAGS  = $(ALL_MODS) -fopenmp
CFLAGS += -cpp -ffree-form -ffree-line-length-none -march=native -fimplicit-none -fall-intrinsics -fmax-errors=10

ifneq ('$(DEBUG)','')
	CFLAGS += -Og -g -Wall -Wextra -fbacktrace -pedantic -fbounds-check -Wuninitialized -ffpe-summary=invalid,zero,overflow,underflow -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan
else
	CFLAGS += -O3 -fno-backtrace -ffpe-summary=none
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

mod_fluid.o :
mod_elements.o :
mod_mesh.o :
mod_num_par.o :

mod_vtk.o : mod_ir_vtk.o
mod_data.o : mod_fluid.o

mod_film.o : mod_data.o mod_num_par.o mod_mesh.o mod_elements.o
mod_ms_film.o : mod_film.o

mod_inout_files.o : mod_vtk.o mod_ms_film.o
mod_test_musst.o : mod_inout_files.o

main.o : mod_test_musst.o

#--------------------------------------------------
utils:
	cd $(utils_dir) ; make all

digis:
	cd $(digis_dir) ; make all

msolv:
	cd $(msolv_dir) ; make all

splin:
	cd $(splin_dir) ; make all

all:
	make utils
	make digis
	make msolv
	make splin
	make clean
	make

debug:
	make "DEBUG=TRUE"

gprof:
	make "GPROF=TRUE"

#--------------------------------------------------
.PHONY: clean debug prod gprof all

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(MODDIR)/*.mod
	rm -f $(EXE)


